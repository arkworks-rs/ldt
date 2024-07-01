use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree, Path},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};
use ark_std::marker::PhantomData;

use crate::{
    domain::Domain,
    ldt::Prover,
    poly_utils,
    stir::{
        config::STIRConfig,
        proof::{STIRFinalRoundProof, STIRInnerRoundProof, STIRProof},
    },
    utils::{dedup, proof_of_work, squeeze_integer, stack_evaluations},
    witness::Witness,
};

pub struct STIRRoundState<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    domain: Domain<F>,
    polynomial: DensePolynomial<F>,
    p_commitment: MerkleTree<M>,
    p_evaluations: Vec<Vec<F>>,
    folding_randomness: F,
    round_num: usize,
    sponge: S,
}

pub struct STIRProver<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    config: STIRConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}

impl<F: FftField + PrimeField + Absorb, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge>
    Prover<F> for STIRProver<F, M, S>
where
    M::InnerDigest: Absorb,
{
    type Config = STIRConfig<M, S>;
    type Proof = STIRProof<F, M>;

    fn new(config: STIRConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }

    fn prove(&self, witness: impl Witness<F>) -> Self::Proof {
        assert!(witness.polynomial().degree() < self.config.starting_degree);

        // get evaluations over a domain
        let domain: Domain<F> =
            Domain::<F>::new(self.config.starting_degree, self.config.starting_rate).unwrap();
        // let evals: Vec<F> = polynomials[0]
        //     .evaluate_over_domain_by_ref(domain.backing_domain)
        //     .evals;
        // let p_evaluations: Vec<Vec<F>> = utils::stack_evaluations(evals, folding_factor);
        let committed_values =
            witness.folded_evaluations_over_domain(domain.clone(), self.config.folding_factor);

        // generate the committment
        let p_commitment = MerkleTree::<M>::new(
            &self.config.merkle_leaf_hash_param,
            &self.config.merkle_two_to_one_param,
            &committed_values,
        )
        .unwrap();

        // Step 1: get initial state of the protocol
        let mut current_round_state: STIRRoundState<F, M, S> =
            Self::get_round_state_from_commitment(
                &p_commitment,
                committed_values,
                &domain,
                witness.polynomial(),
                &self.config.sponge_config,
            );

        // Step 2: compute inner rounds
        let mut inner_round_proofs: Vec<STIRInnerRoundProof<F, M>> =
            Vec::with_capacity(self.config.num_rounds);
        for _round in 0..self.config.num_rounds {
            let (round_state, round_proof) =
                Self::compute_inner_round(&self.config, current_round_state);
            current_round_state = round_state;
            inner_round_proofs.push(round_proof);
        }

        // Step 3: compute final round (v similar but fewer things)
        let final_round_proof = Self::compute_final_round(&self.config, current_round_state);

        // Boom.
        Self::Proof {
            initial_p_commitment_root: p_commitment.root(),
            inner_round_proofs,
            final_round_proof,
        }
    }
}

impl<F: FftField + PrimeField + Absorb, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge>
    STIRProver<F, M, S>
where
    M::InnerDigest: Absorb,
{
    fn get_round_state_from_commitment(
        commitment: &MerkleTree<M>,
        committed_values: Vec<Vec<F>>,
        domain: &Domain<F>,
        polynomial: DensePolynomial<F>,
        sponge_config: &S::Config,
    ) -> STIRRoundState<F, M, S> {
        let mut sponge = S::new(sponge_config);
        sponge.absorb(&commitment.root());
        let folding_randomness = sponge.squeeze_field_elements(1)[0];
        STIRRoundState {
            domain: domain.clone(),
            polynomial: polynomial.clone(),
            p_commitment: commitment.clone(),
            p_evaluations: committed_values.clone(),
            folding_randomness,
            round_num: 0,
            sponge,
        }
    }
    fn compute_final_round(
        config: &STIRConfig<M, S>,
        mut round_state: STIRRoundState<F, M, S>,
    ) -> STIRFinalRoundProof<F, M> {
        // Step 1: Perfom fold operation
        let polynomial = poly_utils::folding::poly_fold(
            &round_state.polynomial,
            config.folding_factor,
            round_state.folding_randomness,
        );

        // Step 2: Generate challenges and answers
        let (_, leaf_values_of_queries, inclusion_proofs_of_queries) = Self::generate_sampling(
            &mut round_state.sponge,
            round_state.domain,
            round_state.p_commitment,
            round_state.p_evaluations,
            config.repetitions[config.num_rounds],
            config.folding_factor,
        );

        // Step 3: Proof of work
        let proof_of_work_nonce = proof_of_work(
            &mut round_state.sponge,
            config.proof_of_work_bits[config.num_rounds],
        );

        // Boom.
        STIRFinalRoundProof {
            polynomial,
            leaf_values_of_queries,
            inclusion_proofs_of_queries,
            proof_of_work_nonce,
        }
    }
    fn compute_inner_round(
        config: &STIRConfig<M, S>,
        mut round_state: STIRRoundState<F, M, S>,
    ) -> (STIRRoundState<F, M, S>, STIRInnerRoundProof<F, M>) {
        // Step 1: Perform fold/scale operation
        let (folded_polynomial, mut scaled_domain, folded_evaluations) = Self::fold_polynomial(
            round_state.polynomial.clone(),
            round_state.domain.clone(),
            config.folding_factor,
            round_state.folding_randomness,
        );

        // Step 2: Generate commitment using a Merkle Tree
        let folded_p_commitment = MerkleTree::<M>::new(
            &config.merkle_leaf_hash_param,
            &config.merkle_two_to_one_param,
            &folded_evaluations,
        )
        .unwrap();
        let folded_p_commitment_root = folded_p_commitment.root();
        round_state.sponge.absorb(&folded_p_commitment_root);

        // Step 3: Out of domain sampling
        let (out_of_domain_samples, out_of_domain_evaluations) =
            Self::get_out_of_domain_evaluations(
                &mut round_state.sponge,
                folded_polynomial.clone(),
                config.num_out_of_domain_samples,
            );
        round_state.sponge.absorb(&out_of_domain_evaluations);

        // Step 4: Squeeze some randomness
        let proximity_generator_randomness: F = round_state.sponge.squeeze_field_elements(1)[0];
        let next_round_folding_randomness: F = round_state.sponge.squeeze_field_elements(1)[0];

        // Step 5: Generate challenges and answers
        let (random_queries, leaf_values_of_queries, inclusion_proofs_of_queries) =
            Self::generate_sampling(
                &mut round_state.sponge,
                round_state.domain,
                round_state.p_commitment,
                round_state.p_evaluations,
                config.repetitions[round_state.round_num],
                config.folding_factor,
            );

        // Step 6: Proof of work
        let proof_of_work_nonce = proof_of_work(
            &mut round_state.sponge,
            config.proof_of_work_bits[round_state.round_num],
        );

        // Step 7: Squeeze more randomness (used by only verifier)
        let _shake_randomness: F = round_state.sponge.squeeze_field_elements(1)[0];

        // Step 6: Generate quotient set and answers
        let (quotient_set, quotient_answers) = Self::get_quotient_set_and_answers(
            &mut scaled_domain,
            folded_polynomial.clone(),
            random_queries,
            out_of_domain_samples,
            config.folding_factor,
        );

        // Step 7: Compute polynomials
        let (answer_polynomial, shake_polynomial, witness_polynomial) = Self::compute_polynomials(
            quotient_set,
            quotient_answers,
            folded_polynomial,
            proximity_generator_randomness,
        );

        // Boom.
        (
            STIRRoundState {
                domain: scaled_domain,
                polynomial: witness_polynomial,
                p_commitment: folded_p_commitment,
                p_evaluations: folded_evaluations,
                folding_randomness: next_round_folding_randomness,
                round_num: round_state.round_num + 1,
                sponge: round_state.sponge,
            },
            STIRInnerRoundProof {
                p_commitment_root: folded_p_commitment_root,
                out_of_domain_evaluations,
                leaf_values_of_queries,
                inclusion_proofs_of_queries,
                answer_polynomial,
                shake_polynomial,
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
    fn get_inclusion_proofs(p_commitment: MerkleTree<M>, leaf_indices: Vec<usize>) -> Vec<Path<M>> {
        // TODO: change this back to multiproof API
        let mut inclusion_proofs: Vec<Path<M>> = Vec::with_capacity(leaf_indices.len());
        for leaf_index in leaf_indices {
            inclusion_proofs.push(p_commitment.generate_proof(leaf_index).unwrap());
        }
        inclusion_proofs
    }
    fn generate_sampling(
        sponge: &mut S,
        domain: Domain<F>,
        p_commitment: MerkleTree<M>,
        evaluations: Vec<Vec<F>>,
        num_repetitions: usize,
        folding_factor: usize,
    ) -> (Vec<usize>, Vec<Vec<F>>, Vec<Path<M>>) {
        let random_queries: Vec<usize> =
            Self::squeeze_queries(sponge, domain.clone(), num_repetitions, folding_factor);
        let leaf_values_of_queries: Vec<Vec<F>> =
            Self::get_leaf_values_from_queries(random_queries.clone(), evaluations);
        let inclusion_proofs_of_queries: Vec<Path<M>> =
            Self::get_inclusion_proofs(p_commitment, random_queries.clone());
        (
            random_queries,
            leaf_values_of_queries,
            inclusion_proofs_of_queries,
        )
    }
    fn fold_polynomial(
        polynomial: DensePolynomial<F>,
        domain: Domain<F>,
        folding_factor: usize,
        folding_randomness: F,
    ) -> (DensePolynomial<F>, Domain<F>, Vec<Vec<F>>) {
        let folded_p =
            poly_utils::folding::poly_fold(&polynomial, folding_factor, folding_randomness);
        // TODO: there is some other option than FFT?
        let scaled_domain = domain.scale_offset(2);
        let evaluations = folded_p
            .evaluate_over_domain_by_ref(scaled_domain.backing_domain)
            .evals;
        let folded_evaluations = stack_evaluations(evaluations, folding_factor);
        (folded_p, scaled_domain, folded_evaluations)
    }
    fn get_out_of_domain_evaluations(
        sponge: &mut S,
        polynomial: DensePolynomial<F>,
        num_samples: usize,
    ) -> (Vec<F>, Vec<F>) {
        let out_of_domain_samples: Vec<F> = sponge.squeeze_field_elements(num_samples);
        let evaluations: Vec<F> = out_of_domain_samples
            .iter()
            .map(|sample| polynomial.evaluate(sample))
            .collect();
        (out_of_domain_samples, evaluations)
    }
    fn get_quotient_set_and_answers(
        domain: &mut Domain<F>,
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
