use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree, Path},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};
use ark_std::marker::PhantomData;

use crate::{
    domain::Domain,
    ldt::{
        stir::{
            config::STIRConfig,
            proof::{STIRProof, STIRRoundProof},
            state::STIRRoundState,
        },
        Prover,
    },
    poly_utils,
    utils::{dedup, proof_of_work, squeeze_integer},
    witness::Witness,
};

pub struct STIRProver<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    W: Witness<F, M, MerkleConfig = M>,
{
    config: STIRConfig<W::MerkleConfig, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<W::MerkleConfig>,
    _sponge_config: PhantomData<S>,
}

impl<F, M, S, W> Prover<F> for STIRProver<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M, Commitment = MerkleTree<M>, CommittedValues = Vec<Vec<F>>>
        + Clone,
    W::ChallengeAnswers: Clone,
{
    type Witness = W;
    type ProverConfig = STIRConfig<M, S>;
    type Proof = STIRProof<F, M>;

    fn new(config: STIRConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }

    fn prove(&self, witness: &W) -> Self::Proof {
        // TODO: (z-tech) what does this check mean?
        assert!(witness.coeff().degree() < self.config.starting_degree);

        // Step 1: get initial round state
        let mut round_state = STIRRoundState::new(
            witness.domain(),
            witness.coeff(),
            witness.commitment(),
            witness.committed_values(),
            self.config.sponge_config.clone(),
        );

        // Step 2: compute inner rounds
        let mut round_proofs = Vec::with_capacity(self.config.num_rounds);
        for _round in 0..self.config.num_rounds {
            let (new_round_state, round_proof) =
                Self::compute_inner_round(&self.config, round_state);
            round_state = new_round_state;
            round_proofs.push(round_proof);
        }

        // Step 3: compute final round (v similar but fewer things)
        let final_round_proof = Self::compute_final_round(&self.config, round_state);
        round_proofs.push(final_round_proof.clone());

        // Boom.
        STIRProof::<F, M> { round_proofs }
    }
}

impl<F, M, S, W> STIRProver<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    W: Witness<F, M, MerkleConfig = M>,
{
    fn compute_final_round(
        config: &STIRConfig<W::MerkleConfig, S>,
        mut round_state: STIRRoundState<F, W::MerkleConfig, S>,
    ) -> STIRRoundProof<F, W::MerkleConfig> {
        // Step 1: Perfom fold operation
        let coeff = poly_utils::folding::poly_fold(
            &round_state.coeff,
            config.folding_factor,
            round_state.folding_randomness,
        );

        // Step 2: Generate challenges and answers
        let (_, committed_values, challenge_answers) = Self::generate_sampling(
            &mut round_state.sponge,
            round_state.domain,
            round_state.commitment.clone(),
            round_state.committed_values,
            config.repetitions[config.num_rounds],
            config.folding_factor,
        );

        // Step 3: Proof of work
        let proof_of_work_nonce = proof_of_work(
            &mut round_state.sponge,
            config.proof_of_work_bits[config.num_rounds],
        );

        // Step 4: Return
        STIRRoundProof::new(
            coeff,
            challenge_answers,
            committed_values,
            true,
            vec![],
            round_state.commitment.root(),
            proof_of_work_nonce,
            DensePolynomial::<F>::from_coefficients_vec(Vec::new()),
        )
    }
    fn compute_inner_round(
        config: &STIRConfig<W::MerkleConfig, S>,
        mut round_state: STIRRoundState<F, W::MerkleConfig, S>,
    ) -> (
        STIRRoundState<F, W::MerkleConfig, S>,
        STIRRoundProof<F, W::MerkleConfig>,
    ) {
        let last_round_domain = round_state.domain.clone();
        let last_round_commitment = round_state.commitment.clone();
        let last_round_committed_values = round_state.committed_values.clone();

        // Step 1: Perform fold/scale operation
        round_state.fold(config.folding_factor);

        // Step 2: Generate commitment on the folded stuff
        let round_commitment = MerkleTree::<W::MerkleConfig>::new(
            &config.merkle_leaf_hash_param,
            &config.merkle_two_to_one_param,
            &round_state.committed_values,
        )
        .unwrap();
        // put it in the sponge
        let round_commitment_digest = round_commitment.root();
        round_state.sponge_absorb(&round_commitment_digest);

        // Step 3: Out of domain samples
        let (out_of_domain_samples, out_of_domain_evaluations) =
            Self::out_of_domain_sample(&mut round_state, config.num_out_of_domain_samples);
        // put it in the sponge
        round_state.sponge_absorb(&out_of_domain_evaluations);

        // Step 4: Squeeze some randomness
        let proximity_generator_randomness = round_state.sponge_squeeze();
        round_state.update_folding_randomness();

        // Step 5: Generate challenges and answers
        let (challenges, committed_values, challenge_answers) = Self::generate_sampling(
            &mut round_state.sponge,
            last_round_domain,
            last_round_commitment,
            last_round_committed_values,
            config.repetitions[round_state.round_num],
            config.folding_factor,
        );

        // Step 6: Proof of work
        let proof_of_work_nonce = proof_of_work(
            &mut round_state.sponge,
            config.proof_of_work_bits[round_state.round_num],
        );

        // Step 7: Squeeze more randomness (used by only verifier)
        let _shake_randomness: F = round_state.sponge_squeeze();

        // Step 6: Generate quotient set and answers
        let (quotient_set, quotient_answers) = Self::get_quotient_set_and_answers(
            round_state.domain.clone(),
            round_state.coeff.clone(),
            challenges,
            out_of_domain_samples,
            config.folding_factor,
        );

        // Step 7: Compute polynomials
        let (answer_polynomial, shake_polynomial, witness_polynomial) = Self::compute_polynomials(
            quotient_set,
            quotient_answers,
            round_state.coeff,
            proximity_generator_randomness,
        );

        // Step 8: Return
        (
            STIRRoundState {
                domain: round_state.domain,
                coeff: witness_polynomial,
                commitment: round_commitment,
                committed_values: round_state.committed_values,
                folding_randomness: round_state.folding_randomness,
                round_num: round_state.round_num + 1,
                sponge: round_state.sponge,
            },
            STIRRoundProof {
                commitment_digest: round_commitment_digest,
                out_of_domain_evaluations,
                committed_values,
                challenge_answers,
                coeff: answer_polynomial,
                is_final_round: false,
                shake_coeff: shake_polynomial,
                proof_of_work_nonce,
            },
        )
    }
    fn get_leaf_values_from_queries(queries: Vec<usize>, evaluations: Vec<Vec<F>>) -> Vec<Vec<F>> {
        queries
            .iter()
            .map(|index| evaluations[*index].clone())
            .collect()
    }
    fn squeeze_queries(
        sponge: &mut S,
        domain: Domain<F>,
        num_repetitions: usize,
        folding_factor: usize,
    ) -> Vec<usize> {
        let scaling_factor: usize = domain.size() / folding_factor;
        // TODO: how would you get a dupe? And would't you be short one query if you did get one?
        let random_queries: Vec<usize> =
            dedup((0..num_repetitions).map(|_| squeeze_integer(sponge, scaling_factor)));
        return random_queries;
    }
    fn get_inclusion_proofs(
        p_commitment: MerkleTree<W::MerkleConfig>,
        leaf_indices: Vec<usize>,
    ) -> Vec<Path<W::MerkleConfig>> {
        // TODO: change this back to multiproof API
        let mut inclusion_proofs: Vec<Path<W::MerkleConfig>> =
            Vec::with_capacity(leaf_indices.len());
        for leaf_index in leaf_indices {
            inclusion_proofs.push(p_commitment.generate_proof(leaf_index).unwrap());
        }
        inclusion_proofs
    }
    fn generate_sampling(
        sponge: &mut S,
        domain: Domain<F>,
        p_commitment: MerkleTree<W::MerkleConfig>,
        evaluations: Vec<Vec<F>>,
        num_repetitions: usize,
        folding_factor: usize,
    ) -> (Vec<usize>, Vec<Vec<F>>, Vec<Path<W::MerkleConfig>>) {
        let random_queries: Vec<usize> =
            Self::squeeze_queries(sponge, domain.clone(), num_repetitions, folding_factor);
        let leaf_values_of_queries: Vec<Vec<F>> =
            Self::get_leaf_values_from_queries(random_queries.clone(), evaluations);
        let inclusion_proofs_of_queries: Vec<Path<W::MerkleConfig>> =
            Self::get_inclusion_proofs(p_commitment, random_queries.clone());
        (
            random_queries,
            leaf_values_of_queries,
            inclusion_proofs_of_queries,
        )
    }
    fn out_of_domain_sample(
        round_state: &mut STIRRoundState<F, M, S>,
        num_samples: usize,
    ) -> (Vec<F>, Vec<F>) {
        let points: Vec<F> = round_state.sponge_squeeze_multiple(num_samples);
        let evals: Vec<F> = points
            .iter()
            .map(|point| round_state.coeff.evaluate(point))
            .collect();
        (points, evals)
    }
    fn get_quotient_set_and_answers(
        domain: Domain<F>,
        polynomial: DensePolynomial<F>,
        random_queries: Vec<usize>,
        out_of_domain_samples: Vec<F>,
        folding_factor: usize,
    ) -> (Vec<F>, Vec<(F, F)>) {
        let stir_randomness: Vec<F> = random_queries
            .iter()
            .map(|index| domain.scale(folding_factor).element(*index))
            .collect();

        // Then compute the set we are quotienting by
        let quotient_set: Vec<F> = out_of_domain_samples
            .into_iter()
            .chain(stir_randomness.iter().cloned())
            .collect();

        // TODO: We can probably reuse this in quotient
        let quotient_answers: Vec<(F, F)> = quotient_set
            .iter()
            .map(|x| (*x, polynomial.evaluate(x)))
            .collect::<Vec<_>>();
        (quotient_set, quotient_answers)
    }
    fn compute_polynomials(
        quotient_set: Vec<F>,
        quotient_answers: Vec<(F, F)>,
        polynomial: DensePolynomial<F>,
        comb_randomness: F,
    ) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
        // Perform naive interpolation to get the answer polynomial
        let answer_polynomial = poly_utils::interpolation::naive_interpolation(&quotient_answers);

        // Initialize shake_polynomial as an empty polynomial
        let mut shake_polynomial = DensePolynomial::from_coefficients_vec(vec![]);
        for (x, y) in &quotient_answers {
            let num_polynomial =
                &answer_polynomial - &DensePolynomial::from_coefficients_vec(vec![*y]);
            let den_polynomial = DensePolynomial::from_coefficients_vec(vec![-*x, F::ONE]);
            shake_polynomial = shake_polynomial + (&num_polynomial / &den_polynomial);
        }

        // Compute the quotient polynomial
        let quotient_polynomial = poly_utils::quotient::poly_quotient(&polynomial, &quotient_set);

        // Compute the scaling polynomial: 1 + r * x + r^2 * x^2 + ... + r^n * x^n
        let scaling_polynomial = DensePolynomial::from_coefficients_vec(
            (0..=quotient_set.len())
                .map(|i| comb_randomness.pow([i as u64]))
                .collect(),
        );

        // Compute the witness polynomial
        let witness_polynomial = &quotient_polynomial * &scaling_polynomial;

        (answer_polynomial, shake_polynomial, witness_polynomial)
    }
}
