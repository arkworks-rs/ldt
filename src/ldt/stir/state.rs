use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::DensePolynomial;

use crate::{domain::Domain, poly_utils, utils::stack_evaluations};

pub struct STIRRoundState<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    pub domain: Domain<F>,
    pub coeff: DensePolynomial<F>,
    pub commitment: MerkleTree<M>,
    pub committed_values: Vec<Vec<F>>,
    pub folding_randomness: F,
    pub round_num: usize,
    pub sponge: S,
}

impl<F, M, S> STIRRoundState<F, M, S>
where
    F: FftField + PrimeField,
    M: MerkleConfig,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
{
    pub fn new(
        domain: Domain<F>,
        coeff: DensePolynomial<F>,
        commitment: MerkleTree<M>,
        committed_values: Vec<Vec<F>>,
        sponge_config: S::Config,
    ) -> Self {
        let mut sponge = S::new(&sponge_config);
        sponge.absorb(&commitment.root());
        Self {
            domain,
            coeff,
            commitment,
            committed_values,
            folding_randomness: sponge.squeeze_field_elements(1)[0],
            round_num: 0, // TODO: is this needed?
            sponge,
        }
    }
    // pub fn coeff(&self) -> DensePolynomial<F> {
    //     self.coeff.clone()
    // }
    // pub fn commitment(&self) -> MerkleTree<M> {
    //     self.commitment.clone()
    // }
    // pub fn committed_values(&self) -> Vec<Vec<F>> {
    //     self.committed_values.clone()
    // }
    // pub fn domain(&self) -> Domain<F> {
    //     self.domain.clone()
    // }
    pub fn fold(&mut self, folding_factor: usize) {
        let folded_coeff = poly_utils::folding::poly_fold(
            &self.coeff.clone(),
            folding_factor,
            self.folding_randomness,
        );
        let scaled_domain = self.domain.clone().scale_offset(2);
        let evals = folded_coeff
            .evaluate_over_domain_by_ref(scaled_domain.backing_domain)
            .evals;
        let folded_committed_values = stack_evaluations(evals, folding_factor);
        self.coeff = folded_coeff;
        self.domain = scaled_domain;
        self.committed_values = folded_committed_values;
    }
    // pub fn folding_randomness(&self) -> F {
    //     self.folding_randomness
    // }
    // pub fn round_num(&self) -> usize {
    //     self.round_num
    // }
    // pub fn sponge(&self) -> S {
    //     self.sponge
    // }
    pub fn sponge_absorb(&mut self, element: impl Absorb) {
        self.sponge.absorb(&element);
    }
    pub fn sponge_squeeze(&mut self) -> F {
        self.sponge.squeeze_field_elements(1)[0]
    }
    pub fn sponge_squeeze_multiple(&mut self, num_elements: usize) -> Vec<F> {
        self.sponge.squeeze_field_elements(num_elements)
    }
    pub fn update_folding_randomness(&mut self) {
        self.folding_randomness = self.sponge_squeeze();
    }
}
