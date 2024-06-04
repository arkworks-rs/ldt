use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, MerkleTree, Path, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};
use ark_std::marker::PhantomData;

use crate::{
    commitment::Commitment,
    domain::Domain,
    ldt::Prover,
    poly_utils::{self, folding},
    stir::{
        config::STIRConfig,
        proof::{STIRProof, STIRRoundProof},
    },
    utils::{dedup, proof_of_work, squeeze_integer, stack_evaluations},
};

pub struct WitnessExtended<F: FftField, M: MerkleConfig> {
    pub domain: Domain<F>,
    pub polynomial: DensePolynomial<F>,
    pub merkle_tree: MerkleTree<M>,
    pub folded_evals: Vec<Vec<F>>,
    pub num_round: usize,
    pub folding_randomness: F,
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

        let (folding_randomness, mut sponge) = self.sponge_up(commitment.p_commitment.root());

        let mut witness = WitnessExtended {
            domain: commitment.domain.clone(),
            polynomial: commitment.polynomials[0].clone(),
            merkle_tree: commitment.p_commitment.clone(),
            folded_evals: commitment.p_evaluations.clone(),
            num_round: 0,
            folding_randomness,
        };

        let mut round_proofs = vec![];
        for _ in 0..self.config.num_rounds {
            let (new_witness, round_proof) = self.compute_round(&mut sponge, &witness);
            witness = new_witness;
            round_proofs.push(round_proof);
        }

        let final_polynomial = poly_utils::folding::poly_fold(
            &witness.polynomial,
            self.config.folding_factor,
            witness.folding_randomness,
        );

        let (_, leaf_values_of_queries, inclusion_proofs_of_queries) = Self::generate_sampling(
            &mut sponge,
            witness.domain.clone(),
            witness.merkle_tree.clone(),
            witness.folded_evals.clone(),
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
    fn compute_round(
        &self,
        sponge: &mut S,
        witness: &WitnessExtended<F, M>,
    ) -> (WitnessExtended<F, M>, STIRRoundProof<F, M>) {
        // 1. perform fold / scale
        let (folded_polynomial, scaled_domain, folded_evaluations) = Self::fold_polynomial(
            witness.polynomial.clone(),
            witness.domain.clone(),
            self.config.folding_factor,
            witness.folding_randomness,
        );

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
                sponge,
                folded_polynomial.clone(),
                self.config.num_out_of_domain_samples,
            );
        sponge.absorb(&out_of_domain_evaluations);

        // TODO is there a reason these occur here rather than immediately before their usage?
        // Proximity generator
        let comb_randomness: F = sponge.squeeze_field_elements(1)[0];
        // Folding randomness for next round
        let folding_randomness = sponge.squeeze_field_elements(1)[0];

        // 4. generate challenges and answers
        // The verifier queries the previous oracle at the indexes of L^k (reading the corresponding evals)
        let (random_queries, leaf_values_of_queries, inclusion_proofs_of_queries) =
            Self::generate_sampling(
                sponge,
                witness.domain.clone(),
                witness.merkle_tree.clone(),
                witness.folded_evals.clone(),
                self.config.repetitions[witness.num_round],
                self.config.folding_factor,
            );
        let queries_to_prev = (leaf_values_of_queries, inclusion_proofs_of_queries);

        let pow_nonce = proof_of_work(sponge, self.config.proof_of_work_bits[witness.num_round]);

        // Not used
        let _shake_randomness: F = sponge.squeeze_field_elements(1)[0];

        // Here, we update the witness
        // First, compute the set of points we are actually going to query at
        let stir_randomness: Vec<_> = random_queries
            .iter()
            .map(|index| {
                witness
                    .domain
                    .scale(self.config.folding_factor)
                    .element(*index)
            })
            .collect();

        // Then compute the set we are quotienting by
        let quotient_set: Vec<_> = out_of_domain_samples
            .into_iter()
            .chain(stir_randomness.iter().cloned())
            .collect();

        // TODO: We can probably reuse this in quotient
        let quotient_answers = quotient_set
            .iter()
            .map(|x| (*x, folded_polynomial.evaluate(x)))
            .collect::<Vec<_>>();

        let ans_polynomial = poly_utils::interpolation::naive_interpolation(&quotient_answers);

        let mut shake_polynomial = DensePolynomial::from_coefficients_vec(vec![]);
        for (x, y) in quotient_answers {
            let num_polynomial = &ans_polynomial - &DensePolynomial::from_coefficients_vec(vec![y]);
            let den_polynomial = DensePolynomial::from_coefficients_vec(vec![-x, F::ONE]);
            shake_polynomial = shake_polynomial + (&num_polynomial / &den_polynomial);
        }

        // The quotient polynomial is then computed
        let quotient_polynomial =
            poly_utils::quotient::poly_quotient(&folded_polynomial, &quotient_set);

        // This is the polynomial 1 + r * x + r^2 * x^2 + ... + r^n * x^n where n = |quotient_set|
        let scaling_polynomial = DensePolynomial::from_coefficients_vec(
            (0..quotient_set.len() + 1)
                .map(|i| comb_randomness.pow([i as u64]))
                .collect(),
        );

        let witness_polynomial = &quotient_polynomial * &scaling_polynomial;

        (
            WitnessExtended {
                domain: scaled_domain,
                polynomial: witness_polynomial,
                merkle_tree: folded_p_commitment,
                folded_evals: folded_evaluations,
                num_round: witness.num_round + 1,
                folding_randomness,
            },
            STIRRoundProof {
                g_root: folded_p_commitment_root,
                betas: out_of_domain_evaluations,
                queries_to_prev,
                ans_polynomial,
                shake_polynomial,
                proof_of_work_nonce: pow_nonce,
            },
        )
    }

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
}
