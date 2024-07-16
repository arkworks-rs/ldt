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
        let tmp_round_num = round_state.round_num;
        let tmp_domain = round_state.domain.clone();
        let challenges = Self::challenges(
            &mut round_state,
            tmp_domain.size() / config.folding_factor,
            config.repetitions[tmp_round_num],
        );
        let (committed_values, challenge_answers) = Self::challenge_answers(
            challenges,
            round_state.commitment.clone(),
            round_state.committed_values,
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
            Self::out_of_domain_samples(&mut round_state, config.num_out_of_domain_samples);
        // put it in the sponge
        round_state.sponge_absorb(&out_of_domain_evaluations);

        // Step 4: Squeeze some randomness
        let proximity_generator_randomness = round_state.sponge_squeeze();
        round_state.update_folding_randomness();

        // Step 5: Generate challenges and answers
        let tmp_round_num = round_state.round_num;
        let challenges = Self::challenges(
            &mut round_state,
            last_round_domain.size() / config.folding_factor,
            config.repetitions[tmp_round_num],
        );
        let (challenge_values, challenge_answers) = Self::challenge_answers(
            challenges.clone(),
            last_round_commitment,
            last_round_committed_values,
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
        let (answer_coeff, shake_coeff, witness_coeff) = Self::compute_polynomials(
            quotient_set,
            quotient_answers,
            round_state.coeff,
            proximity_generator_randomness,
        );

        // Step 8: Return
        let new_round_state = STIRRoundState {
            answer_coeff: answer_coeff.clone(),
            domain: round_state.domain,
            challenge_answers: challenge_answers,
            coeff: witness_coeff, // witness_coeff
            commitment: round_commitment.clone(),
            committed_values: round_state.committed_values,
            folding_randomness: round_state.folding_randomness,
            out_of_domain_evaluations: out_of_domain_evaluations.clone(),
            proof_of_work_nonce: proof_of_work_nonce.clone(),
            round_num: round_state.round_num + 1,
            shake_coeff: shake_coeff.clone(),
            sponge: round_state.sponge,
        };
        let new_round_proof = STIRRoundProof {
            commitment_digest: new_round_state.commitment.root(),
            out_of_domain_evaluations: new_round_state.out_of_domain_evaluations.clone(),
            challenge_values,
            challenge_answers: new_round_state.challenge_answers.clone(),
            coeff: answer_coeff,
            is_final_round: false,
            shake_coeff: new_round_state.shake_coeff.clone(),
            proof_of_work_nonce: new_round_state.proof_of_work_nonce,
        };
        (
            new_round_state,
            new_round_proof,
        )
    }
    fn challenges(
        round_state: &mut STIRRoundState<F, M, S>,
        scaling_factor: usize,
        num_repetitions: usize,
    ) -> Vec<usize> {
        dedup(
            (0..num_repetitions).map(|_| squeeze_integer(&mut round_state.sponge, scaling_factor)),
        )
    }
    fn challenge_answers(
        challenges: Vec<usize>,
        last_round_commitment: MerkleTree<M>,
        last_round_committed_values: Vec<Vec<F>>,
    ) -> (Vec<Vec<F>>, Vec<Path<M>>) {
        let challenge_values: Vec<Vec<F>> = challenges
            .iter()
            .map(|index| last_round_committed_values[*index].clone())
            .collect();
        let mut challenge_answers: Vec<Path<M>> = Vec::with_capacity(challenge_values.len());
        for challenge in challenges {
            challenge_answers.push(last_round_commitment.generate_proof(challenge).unwrap());
        }
        (challenge_values, challenge_answers)
    }
    fn out_of_domain_samples(
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
        coeff: DensePolynomial<F>,
        proximity_generator_randomness: F,
    ) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
        // Perform naive interpolation to get the answer polynomial
        let answer_coeff = poly_utils::interpolation::naive_interpolation(&quotient_answers);

        // Initialize shake_polynomial as an empty polynomial
        let mut shake_coeff = DensePolynomial::from_coefficients_vec(vec![]);
        for (x, y) in &quotient_answers {
            let num_coeff = &answer_coeff - &DensePolynomial::from_coefficients_vec(vec![*y]);
            let den_coeff = DensePolynomial::from_coefficients_vec(vec![-*x, F::ONE]);
            shake_coeff = shake_coeff + (&num_coeff / &den_coeff);
        }

        // Compute the quotient polynomial
        let quotient_coeff = poly_utils::quotient::poly_quotient(&coeff, &quotient_set);

        // Compute the scaling polynomial: 1 + r * x + r^2 * x^2 + ... + r^n * x^n
        let scaling_polynomial = DensePolynomial::from_coefficients_vec(
            (0..=quotient_set.len())
                .map(|i| proximity_generator_randomness.pow([i as u64]))
                .collect(),
        );

        // Compute the witness polynomial
        let witness_coeff = &quotient_coeff * &scaling_polynomial;

        (answer_coeff, shake_coeff, witness_coeff)
    }
}
