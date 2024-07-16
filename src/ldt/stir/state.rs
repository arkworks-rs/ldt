use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree, Path},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};

use crate::{domain::Domain, poly_utils, utils::stack_evaluations};

use super::proof::STIRRoundProof;

pub struct STIRRoundState<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    pub answer_coeff: DensePolynomial<F>,
    pub domain: Domain<F>,
    pub coeff: DensePolynomial<F>,
    pub challenge_answers: Vec<Path<M>>,
    pub challenge_values: Vec<Vec<F>>,
    pub commitment: MerkleTree<M>,
    pub committed_values: Vec<Vec<F>>,
    pub folding_randomness: F,
    pub out_of_domain_evaluations: Vec<F>,
    pub proof_of_work_nonce: Option<usize>,
    pub round_num: usize,
    pub shake_coeff: DensePolynomial<F>,
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
            answer_coeff: DensePolynomial::from_coefficients_vec(vec![]),
            domain,
            challenge_answers: vec![],
            challenge_values: vec![],
            coeff,
            commitment,
            committed_values,
            folding_randomness: sponge.squeeze_field_elements(1)[0],
            out_of_domain_evaluations: vec![],
            proof_of_work_nonce: None,
            round_num: 0, // TODO: is this needed?
            shake_coeff: DensePolynomial::from_coefficients_vec(vec![]),
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
    pub fn round_num(&self) -> usize {
        self.round_num
    }
    pub fn round_proof(&self) -> STIRRoundProof<F, M> {
        STIRRoundProof {
            commitment_digest: self.commitment.root(),
            out_of_domain_evaluations: self.out_of_domain_evaluations.clone(),
            challenge_values: self.challenge_values.clone(),
            challenge_answers: self.challenge_answers.clone(),
            coeff: self.answer_coeff.clone(),
            is_final_round: false,
            shake_coeff: self.shake_coeff.clone(),
            proof_of_work_nonce: self.proof_of_work_nonce,
        }
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
    pub fn update_folding_randomness(&mut self) {
        self.folding_randomness = self.sponge_squeeze();
    }
}
