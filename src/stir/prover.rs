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
        proof::{STIRProof, STIRRoundProof},
    },
    utils::{dedup, proof_of_work, squeeze_integer, stack_evaluations},
};

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

        // Stuff we're gonna use this round
        let mut domain = commitment.domain.clone();
        let mut polynomial = commitment.polynomials[0].clone();
        let mut p_commitment = commitment.p_commitment.clone();
        let mut p_evaluations = commitment.p_evaluations.clone();
        let (mut folding_randomness, mut sponge) = self.sponge_up(p_commitment.root());

        // For each round
        let mut round_proofs = vec![];
        for round_num in 0..self.config.num_rounds {
            // 1. perform fold / scale
            let (folded_polynomial, mut scaled_domain, folded_evaluations) = Self::fold_polynomial(
                polynomial.clone(),
                domain.clone(),
                self.config.folding_factor,
                folding_randomness,
            );
            // used by final round except for domain

            // 2. generate commitment
            let folded_p_commitment = MerkleTree::<M>::new(
                &self.config.merkle_leaf_hash_param,
                &self.config.merkle_two_to_one_param,
                &folded_evaluations,
            )
            .unwrap();
            let folded_p_commitment_root = folded_p_commitment.root();
            sponge.absorb(&folded_p_commitment_root);

            // 3. out of domain sampling
            let (out_of_domain_samples, out_of_domain_evaluations) =
                Self::get_out_of_domain_evaluations(
                    &mut sponge,
                    folded_polynomial.clone(),
                    self.config.num_out_of_domain_samples,
                );
            sponge.absorb(&out_of_domain_evaluations);

            // TODO is there a reason these occur here rather than immediately before their usage?
            // Proximity generator
            let comb_randomness: F = sponge.squeeze_field_elements(1)[0];
            // Folding randomness for next round_num
            folding_randomness = sponge.squeeze_field_elements(1)[0];

            // 4. generate challenges and answers
            // The verifier queries the previous oracle at the indexes of L^k (reading the corresponding evals)
            let (random_queries, leaf_values_of_queries, inclusion_proofs_of_queries) =
                Self::generate_sampling(
                    &mut sponge,
                    domain,
                    p_commitment,
                    p_evaluations,
                    self.config.repetitions[round_num],
                    self.config.folding_factor,
                ); // used by final round
            let queries_to_prev = (leaf_values_of_queries, inclusion_proofs_of_queries);

            // 5. Proof of work
            let pow_nonce = proof_of_work(&mut sponge, self.config.proof_of_work_bits[round_num]); // used by final round

            // Not used
            let _shake_randomness: F = sponge.squeeze_field_elements(1)[0];

            // 6. Generate quotient set and answers
            let (quotient_set, quotient_answers) = Self::get_quotient_set_and_answers(
                &mut scaled_domain,
                folded_polynomial.clone(),
                random_queries,
                out_of_domain_samples,
                self.config.folding_factor,
            );

            // 7. compute polynomials
            let (answer_polynomial, shake_polynomial, witness_polynomial) = Self::compute_polynomials(
                quotient_set,
                quotient_answers,
                folded_polynomial,
                comb_randomness,
            );
            domain = scaled_domain;
            polynomial = witness_polynomial;
            p_commitment = folded_p_commitment;
            p_evaluations = folded_evaluations;
            folding_randomness = folding_randomness;

            round_proofs.push(STIRRoundProof {
                p_commitment_root: folded_p_commitment_root,
                out_of_domain_evaluations,
                queries_to_prev,
                answer_polynomial,
                shake_polynomial,
                proof_of_work_nonce: pow_nonce,
            });
        }

        let final_polynomial = poly_utils::folding::poly_fold(
            &polynomial,
            self.config.folding_factor,
            folding_randomness,
        );

        let (_, leaf_values_of_queries, inclusion_proofs_of_queries) = Self::generate_sampling(
            &mut sponge,
            domain,
            p_commitment,
            p_evaluations,
            self.config.repetitions[self.config.num_rounds],
            self.config.folding_factor,
        );

        let pow_nonce = proof_of_work(
            &mut sponge,
            self.config.proof_of_work_bits[self.config.num_rounds],
        );

        Self::Proof {
            round_proofs,
            polynomial: final_polynomial,
            queries_to_final: (leaf_values_of_queries, inclusion_proofs_of_queries),
            proof_of_work_nonce: pow_nonce,
        }
    }
}

impl<F: FftField + PrimeField + Absorb, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge>
    STIRProver<F, M, S>
where
    M::InnerDigest: Absorb,
{
    fn sponge_up(&self, digest: M::InnerDigest) -> (F, S) {
        let mut sponge = S::new(&self.config.sponge_config);
        sponge.absorb(&digest);
        (sponge.squeeze_field_elements(1)[0], sponge)
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
        let answer_polynomial = poly_utils::interpolation::naive_interpolation(&quotient_answers);

        let mut shake_polynomial = DensePolynomial::from_coefficients_vec(vec![]);
        for (x, y) in quotient_answers {
            let num_polynomial = &answer_polynomial - &DensePolynomial::from_coefficients_vec(vec![y]);
            let den_polynomial = DensePolynomial::from_coefficients_vec(vec![-x, F::ONE]);
            shake_polynomial = shake_polynomial + (&num_polynomial / &den_polynomial);
        }

        // The quotient polynomial is then computed
        let quotient_polynomial = poly_utils::quotient::poly_quotient(&polynomial, &quotient_set);

        // This is the polynomial 1 + r * x + r^2 * x^2 + ... + r^n * x^n where n = |quotient_set|
        let scaling_polynomial = DensePolynomial::from_coefficients_vec(
            (0..quotient_set.len() + 1)
                .map(|i| comb_randomness.pow([i as u64]))
                .collect(),
        );

        let witness_polynomial = &quotient_polynomial * &scaling_polynomial;
        (answer_polynomial, shake_polynomial, witness_polynomial)
    }
}
