use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};

#[cfg(not(feature = "std"))]
use ark_std::{vec, vec::Vec};

use crate::{
    domain::Domain,
    fri::{config::FRIConfig, proof::FRIProof},
    poly_utils,
    utils::{dedup, proof_of_work, squeeze_integer, stack_evaluations},
};

pub struct FRIProverState<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    pub coeff: DensePolynomial<F>,
    pub domain: Domain<F>,
    commitment_digest: M::InnerDigest,
    committed_values: Vec<Vec<F>>,
    config: FRIConfig<M, S>,
    initial_domain_size: usize,
    round_num: usize,
    sponge: S,
}

impl<F, M, S> FRIProverState<F, M, S>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>> + Clone,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
{
    pub fn new(
        coeff: DensePolynomial<F>,
        commitment_digest: M::InnerDigest,
        committed_values: Vec<Vec<F>>,
        config: FRIConfig<M, S>,
        domain: Domain<F>,
    ) -> Self {
        let mut sponge = S::new(&config.sponge_config);
        sponge.absorb(&commitment_digest);
        Self {
            coeff,
            domain: domain.clone(),
            commitment_digest,
            committed_values: committed_values.clone(),
            config: config.clone(),
            initial_domain_size: domain.size(),
            round_num: 0,
            sponge,
        }
    }
    pub fn challenges(&mut self) -> Vec<usize> {
        let num_leaves = self.initial_domain_size / self.config.folding_factor;
        dedup((0..self.config.repetitions).map(|_| squeeze_integer(&mut self.sponge, num_leaves)))
    }
    pub fn folded_evaluations(
        &mut self,
        folding_randomness: F,
        prev_evals: Vec<Vec<F>>,
    ) -> Vec<Vec<F>> {
        // Fold the initial polynomial
        self.coeff = poly_utils::folding::poly_fold(
            &self.coeff,
            self.config.folding_factor,
            folding_randomness,
        );

        // let prev_evals = folded_evals.last().unwrap();

        // The following lines are just precomputations, to avoid having to do inversion
        // and exponentiations in the inner loop
        let domain_size = self.domain.size();
        let generator = self
            .domain
            .backing_domain
            .element(domain_size / self.config.folding_factor);
        let generator_inv = generator.inverse().unwrap();
        let size_inv = F::from(self.config.folding_factor as u64)
            .inverse()
            .unwrap();
        let coset_offsets: Vec<_> = self
            .domain
            .backing_domain
            .elements()
            .take(prev_evals.len())
            .collect();
        let mut counter = F::ONE;
        let scale = self.domain.backing_domain.element(1).inverse().unwrap();
        let mut coset_offsets_inv: Vec<_> = vec![];
        for _ in 0..prev_evals.len() {
            coset_offsets_inv.push(counter);
            counter *= scale;
        }

        // Compute the evalations of the folded polynomial
        let g_evaluations: Vec<_> = prev_evals
            .iter()
            .zip(coset_offsets.into_iter())
            .zip(coset_offsets_inv.into_iter())
            .map(|((e, c), ci)| (e, c, ci))
            .map(|(evals, coset_offset, coset_offset_inv)| {
                poly_utils::interpolation::fft_interpolate(
                    generator,
                    coset_offset,
                    generator_inv,
                    coset_offset_inv,
                    size_inv,
                    evals,
                )
                .evaluate(&folding_randomness)
            })
            .collect();
        stack_evaluations(g_evaluations, self.config.folding_factor)
    }
    pub fn proof_of_work(&mut self) -> Option<usize> {
        proof_of_work(&mut self.sponge, self.config.proof_of_work_bits)
    }
    pub fn sponge_absorb(&mut self, element: impl Absorb) {
        self.sponge.absorb(&element);
    }
    pub fn sponge_squeeze(&mut self) -> F {
        self.sponge.squeeze_field_elements(1)[0]
    }
    pub fn sponge_squeeze_multiple(&mut self, num_elements: usize) -> Vec<F> {
        self.sponge.squeeze_field_elements(num_elements)
    }
}
