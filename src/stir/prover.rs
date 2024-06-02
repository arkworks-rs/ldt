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

pub struct RoundWitness<F: FftField, M: MerkleConfig> {
    pub domain: Domain<F>,
    pub folding_randomness: F,
    pub polynomial: DensePolynomial<F>,
    pub p_evaluations: Vec<Vec<F>>,
    pub p_commitment: MerkleTree<M>,
    pub round_num: usize,
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
        // TODO: fix this
        // assert!(witness.polynomial.degree() < self.parameters.starting_degree);

        let mut sponge = S::new(&self.config.sponge_config);
        sponge.absorb(&commitment.p_commitment.root());
        let folding_randomness = sponge.squeeze_field_elements(1)[0];

        let mut round_witness: RoundWitness<F, M> = RoundWitness {
            domain: commitment.domain.clone(),
            polynomial: commitment.polynomials[0].clone(),
            p_commitment: commitment.p_commitment.clone(),
            p_evaluations: commitment.p_evaluations.clone(),
            round_num: 0,
            folding_randomness,
        };

        let mut round_proofs = vec![];
        for _ in 0..self.config.num_rounds {
            let (new_witness, round_proof) = self.compute_round(&mut sponge, &round_witness);
            round_witness = new_witness;
            round_proofs.push(round_proof);
        }

        let final_polynomial = poly_utils::folding::poly_fold(
            &round_witness.polynomial,
            self.config.folding_factor,
            round_witness.folding_randomness,
        );

        let final_repetitions = self.config.repetitions[self.config.num_rounds];
        let scaling_factor = round_witness.domain.size() / self.config.folding_factor;
        let final_randomness_indexes =
            dedup((0..final_repetitions).map(|_| squeeze_integer(&mut sponge, scaling_factor)));

        let queries_to_final_ans: Vec<_> = final_randomness_indexes
            .iter()
            .map(|index| round_witness.p_evaluations[*index].clone())
            .collect();

        // TODO `generate_multi_proof`` doesn't exist w/ my version of ark_crypto_primitives
        // let queries_to_final_proof = round_witness
        //     .p_commitment
        //     .generate_multi_proof(final_randomness_indexes)
        //     .unwrap();
        let mut queries_to_final_proof: Vec<Path<M>> =
            Vec::with_capacity(final_randomness_indexes.len());
        for query in final_randomness_indexes.clone() {
            queries_to_final_proof.push(round_witness.p_commitment.generate_proof(query).unwrap());
        }

        let queries_to_final: (Vec<Vec<F>>, Vec<Path<M>>) =
            (queries_to_final_ans, queries_to_final_proof);

        let pow_nonce = proof_of_work(
            &mut sponge,
            self.config.proof_of_work_bits[self.config.num_rounds],
        );

        // TODO: where is the commitment?
        Self::Proof {
            round_proofs,
            polynomial: final_polynomial,
            queries_to_final,
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
        sponge: &mut impl CryptographicSponge,
        witness: &RoundWitness<F, M>,
    ) -> (RoundWitness<F, M>, STIRRoundProof<F, M>) {
        let g_poly = poly_utils::folding::poly_fold(
            &witness.polynomial,
            self.config.folding_factor,
            witness.folding_randomness,
        );

        // TODO: For now, I am FFTing
        let g_domain = witness.domain.scale_offset(2);
        let g_evaluations = g_poly
            .evaluate_over_domain_by_ref(g_domain.backing_domain)
            .evals;

        let g_folded_evaluations = stack_evaluations(g_evaluations, self.config.folding_factor);
        let g_merkle = MerkleTree::<M>::new(
            &self.config.merkle_leaf_hash_param,
            &self.config.merkle_two_to_one_param,
            &g_folded_evaluations,
        )
        .unwrap();
        let g_root = g_merkle.root();
        sponge.absorb(&g_root);

        // Out of domain sample
        let ood_randomness = sponge.squeeze_field_elements(self.config.num_out_of_domain_samples);
        let betas = ood_randomness
            .iter()
            .map(|alpha| g_poly.evaluate(alpha))
            .collect();
        sponge.absorb(&betas);

        // Proximity generator
        let comb_randomness: F = sponge.squeeze_field_elements(1)[0];

        // Folding randomness for next round
        let folding_randomness = sponge.squeeze_field_elements(1)[0];

        // Sample the indexes of L^k that we are going to use for querying the previous Merkle tree
        let scaling_factor = witness.domain.size() / self.config.folding_factor;
        let num_repetitions = self.config.repetitions[witness.round_num];
        let stir_randomness_indexes =
            dedup((0..num_repetitions).map(|_| squeeze_integer(sponge, scaling_factor)));

        let pow_nonce = proof_of_work(sponge, self.config.proof_of_work_bits[witness.round_num]);

        // Not used
        let _shake_randomness: F = sponge.squeeze_field_elements(1)[0];

        // The verifier queries the previous oracle at the indexes of L^k (reading the
        // corresponding evals)
        let queries_to_prev_ans: Vec<_> = stir_randomness_indexes
            .iter()
            .map(|&index| witness.p_evaluations[index].clone())
            .collect();

        // TODO `generate_multi_proof`` doesn't exist w/ my version of ark_crypto_primitives
        // let queries_to_prev_proof = witness
        //     .p_commitment
        //     .generate_multi_proof(stir_randomness_indexes.clone())
        //     .unwrap();
        let mut queries_to_prev_proof: Vec<Path<M>> =
            Vec::with_capacity(stir_randomness_indexes.len());
        for query in stir_randomness_indexes.clone() {
            queries_to_prev_proof.push(witness.p_commitment.generate_proof(query).unwrap());
        }
        let queries_to_prev = (queries_to_prev_ans, queries_to_prev_proof);

        // Here, we update the witness
        // First, compute the set of points we are actually going to query at
        let stir_randomness: Vec<_> = stir_randomness_indexes
            .iter()
            .map(|index| {
                witness
                    .domain
                    .scale(self.config.folding_factor)
                    .element(*index)
            })
            .collect();

        // Then compute the set we are quotienting by
        let quotient_set: Vec<_> = ood_randomness
            .into_iter()
            .chain(stir_randomness.iter().cloned())
            .collect();

        // TODO: We can probably reuse this in quotient
        let quotient_answers = quotient_set
            .iter()
            .map(|x| (*x, g_poly.evaluate(x)))
            .collect::<Vec<_>>();

        let ans_polynomial = poly_utils::interpolation::naive_interpolation(&quotient_answers);

        let mut shake_polynomial = DensePolynomial::from_coefficients_vec(vec![]);
        for (x, y) in quotient_answers {
            let num_polynomial = &ans_polynomial - &DensePolynomial::from_coefficients_vec(vec![y]);
            let den_polynomial = DensePolynomial::from_coefficients_vec(vec![-x, F::ONE]);
            shake_polynomial = shake_polynomial + (&num_polynomial / &den_polynomial);
        }

        // The quotient polynomial is then computed
        let quotient_polynomial = poly_utils::quotient::poly_quotient(&g_poly, &quotient_set);

        // This is the polynomial 1 + r * x + r^2 * x^2 + ... + r^n * x^n where n = |quotient_set|
        let scaling_polynomial = DensePolynomial::from_coefficients_vec(
            (0..quotient_set.len() + 1)
                .map(|i| comb_randomness.pow([i as u64]))
                .collect(),
        );

        let witness_polynomial = &quotient_polynomial * &scaling_polynomial;

        (
            RoundWitness {
                domain: g_domain,
                polynomial: witness_polynomial,
                p_commitment: g_merkle,
                p_evaluations: g_folded_evaluations,
                round_num: witness.round_num + 1,
                folding_randomness,
            },
            STIRRoundProof {
                g_root,
                betas,
                queries_to_prev,
                ans_polynomial,
                shake_polynomial,
                proof_of_work_nonce: pow_nonce,
            },
        )
    }
}
