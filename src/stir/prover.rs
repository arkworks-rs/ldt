use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree, Path},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};
use ark_std::marker::PhantomData;

use crate::{
    commitment::Commitment,
    domain::Domain,
    ldt::Prover,
    poly_utils,
    stir::{
        config::STIRConfig,
        proof::{STIRFinalRoundProof, STIRInnerRoundProof, STIRProof},
    },
    utils::{dedup, proof_of_work, squeeze_integer, stack_evaluations},
};

pub struct STIRRoundState<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    domain: Domain<F>,
    polynomial: DensePolynomial<F>,
    p_commitment: MerkleTree<M>,
    p_evaluations: Vec<Vec<F>>,
    folding_randomness: F,
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
    type Commitment = Commitment<F, M>;
    type Proof = STIRProof<F, M>;

    fn new(config: STIRConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }

    fn prove(&self, commitment: &Self::Commitment) -> Self::Proof {
        assert!(commitment.polynomials[0].degree() < self.config.starting_degree);

        let round_state = self.generate_round_state_from_commitment(commitment);
        let (mut final_round_state, inner_round_proofs): (
            STIRRoundState<F, M, S>,
            Vec<STIRInnerRoundProof<F, M>>,
        ) = self.compute_inner_rounds(round_state);

        // let final_round_proof = self.compute_final_round(domain, polynomial, p_commitment, p_evaluations, folding_randomness, sponge);

        let final_polynomial = poly_utils::folding::poly_fold(
            &final_round_state.polynomial,
            self.config.folding_factor,
            final_round_state.folding_randomness,
        );

        let (_, leaf_values_of_queries, inclusion_proofs_of_queries) = Self::generate_sampling(
            &mut final_round_state.sponge,
            final_round_state.domain,
            final_round_state.p_commitment,
            final_round_state.p_evaluations,
            self.config.repetitions[self.config.num_rounds],
            self.config.folding_factor,
        );

        let pow_nonce = proof_of_work(
            &mut final_round_state.sponge,
            self.config.proof_of_work_bits[self.config.num_rounds],
        );

        let final_round_proof = STIRFinalRoundProof {
            polynomial: final_polynomial,
            leaf_values_of_queries,
            inclusion_proofs_of_queries,
            proof_of_work_nonce: pow_nonce,
        };

        Self::Proof {
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
    fn generate_round_state_from_commitment(
        &self,
        commitment: &Commitment<F, M>,
    ) -> STIRRoundState<F, M, S> {
        let mut sponge = S::new(&self.config.sponge_config);
        sponge.absorb(&commitment.p_commitment.root());
        let folding_randomness = sponge.squeeze_field_elements(1)[0];
        STIRRoundState {
            domain: commitment.domain.clone(),
            polynomial: commitment.polynomials[0].clone(),
            p_commitment: commitment.p_commitment.clone(),
            p_evaluations: commitment.p_evaluations.clone(),
            folding_randomness,
            sponge,
        }
    }
    fn compute_inner_rounds(
        &self,
        mut round_state: STIRRoundState<F, M, S>,
    ) -> (STIRRoundState<F, M, S>, Vec<STIRInnerRoundProof<F, M>>) {
        // For each round
        let mut inner_round_proofs = vec![];
        for round_num in 0..self.config.num_rounds {
            // 1. perform fold / scale
            let (folded_polynomial, mut scaled_domain, folded_evaluations) = Self::fold_polynomial(
                round_state.polynomial.clone(),
                round_state.domain.clone(),
                self.config.folding_factor,
                round_state.folding_randomness,
            );

            // 2. generate commitment
            let folded_p_commitment = MerkleTree::<M>::new(
                &self.config.merkle_leaf_hash_param,
                &self.config.merkle_two_to_one_param,
                &folded_evaluations,
            )
            .unwrap();
            let folded_p_commitment_root = folded_p_commitment.root();
            round_state.sponge.absorb(&folded_p_commitment_root);

            // 3. out of domain sampling
            let (out_of_domain_samples, out_of_domain_evaluations) =
                Self::get_out_of_domain_evaluations(
                    &mut round_state.sponge,
                    folded_polynomial.clone(),
                    self.config.num_out_of_domain_samples,
                );
            round_state.sponge.absorb(&out_of_domain_evaluations);

            // TODO is there a reason these occur here rather than immediately before their usage?
            // Proximity generator
            let comb_randomness: F = round_state.sponge.squeeze_field_elements(1)[0];
            // Folding randomness for next round_num
            let new_folding_randomness: F = round_state.sponge.squeeze_field_elements(1)[0];

            // 4. generate challenges and answers
            // The verifier queries the previous oracle at the indexes of L^k (reading the corresponding evals)
            let (random_queries, leaf_values_of_queries, inclusion_proofs_of_queries) =
                Self::generate_sampling(
                    &mut round_state.sponge,
                    round_state.domain,
                    round_state.p_commitment,
                    round_state.p_evaluations,
                    self.config.repetitions[round_num],
                    self.config.folding_factor,
                ); // used by final round

            // 5. Proof of work
            let pow_nonce = proof_of_work(
                &mut round_state.sponge,
                self.config.proof_of_work_bits[round_num],
            ); // used by final round

            // Not used
            let _shake_randomness: F = round_state.sponge.squeeze_field_elements(1)[0];

            // 6. Generate quotient set and answers
            let (quotient_set, quotient_answers) = Self::get_quotient_set_and_answers(
                &mut scaled_domain,
                folded_polynomial.clone(),
                random_queries,
                out_of_domain_samples,
                self.config.folding_factor,
            );

            // 7. compute polynomials
            let (answer_polynomial, shake_polynomial, witness_polynomial) =
                Self::compute_polynomials(
                    quotient_set,
                    quotient_answers,
                    folded_polynomial,
                    comb_randomness,
                );
            round_state.domain = scaled_domain;
            round_state.polynomial = witness_polynomial;
            round_state.p_commitment = folded_p_commitment;
            round_state.p_evaluations = folded_evaluations;
            round_state.folding_randomness = new_folding_randomness;

            inner_round_proofs.push(STIRInnerRoundProof {
                p_commitment_root: folded_p_commitment_root,
                out_of_domain_evaluations,
                leaf_values_of_queries,
                inclusion_proofs_of_queries,
                answer_polynomial,
                shake_polynomial,
                proof_of_work_nonce: pow_nonce,
            });
        }
        return (round_state, inner_round_proofs);
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
