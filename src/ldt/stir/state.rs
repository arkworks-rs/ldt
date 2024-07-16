use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, MerkleTree, Path, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};

use crate::{
    domain::Domain,
    poly_utils,
    utils::{dedup, proof_of_work, squeeze_integer, stack_evaluations},
};

use super::proof::STIRRoundProof;

pub struct STIRRoundState<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    pub answer_coeff: DensePolynomial<F>,
    pub domain: Domain<F>,
    pub challenge_answers: Vec<Path<M>>,
    pub challenge_values: Vec<Vec<F>>,
    pub challenges: Vec<usize>,
    pub commitment: MerkleTree<M>,
    pub committed_values: Vec<Vec<F>>,
    pub folding_factor: usize,
    pub folding_randomness: F,
    pub last_round_commitment: MerkleTree<M>,
    pub last_round_committed_values: Vec<Vec<F>>,
    pub last_round_domain_size: usize,
    pub merkle_leaf_hash_param: LeafParam<M>,
    pub merkle_two_to_one_param: TwoToOneParam<M>,
    pub num_out_of_domain_samples: usize,
    pub num_proof_of_work_bits: Vec<usize>,
    pub num_repetitions: Vec<usize>,
    pub out_of_domain_samples: Vec<F>,
    pub out_of_domain_evaluations: Vec<F>,
    pub proof_of_work_nonce: Option<usize>,
    pub proximity_generator_randomness: F,
    pub quotient_answers: Vec<F>,
    pub quotient_set: Vec<F>,
    pub round_num: usize,
    pub shake_coeff: DensePolynomial<F>,
    pub sponge: S,
    pub witness_coeff: DensePolynomial<F>,
}

impl<F, M, S> STIRRoundState<F, M, S>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
{
    pub fn new(
        domain: Domain<F>,
        commitment: MerkleTree<M>,
        committed_values: Vec<Vec<F>>,
        folding_factor: usize,
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        num_out_of_domain_samples: usize,
        num_proof_of_work_bits: Vec<usize>,
        num_repetitions: Vec<usize>,
        sponge_config: S::Config,
        witness_coeff: DensePolynomial<F>,
    ) -> Self {
        let mut sponge = S::new(&sponge_config);
        sponge.absorb(&commitment.root());
        Self {
            answer_coeff: DensePolynomial::from_coefficients_vec(vec![]),
            domain: domain.clone(),
            challenge_answers: vec![],
            challenge_values: vec![],
            challenges: vec![],
            commitment,
            committed_values: committed_values.clone(),
            folding_factor,
            folding_randomness: sponge.squeeze_field_elements(1)[0],
            last_round_domain_size: domain.size(),
            last_round_commitment: MerkleTree::<M>::new(
                &merkle_leaf_hash_param,
                &merkle_two_to_one_param,
                &committed_values,
            )
            .unwrap(),
            last_round_committed_values: committed_values,
            merkle_leaf_hash_param,
            merkle_two_to_one_param,
            num_out_of_domain_samples,
            num_proof_of_work_bits,
            num_repetitions,
            out_of_domain_samples: vec![],
            out_of_domain_evaluations: vec![],
            proof_of_work_nonce: None,
            proximity_generator_randomness: F::one(),
            quotient_answers: vec![],
            quotient_set: vec![],
            round_num: 0,
            shake_coeff: DensePolynomial::from_coefficients_vec(vec![]),
            sponge,
            witness_coeff,
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
        self.last_round_domain_size = self.domain.size();
        let folded_coeff = poly_utils::folding::poly_fold(
            &self.witness_coeff.clone(),
            folding_factor,
            self.folding_randomness,
        );
        let scaled_domain = self.domain.clone().scale_offset(2);
        let evals = folded_coeff
            .evaluate_over_domain_by_ref(scaled_domain.backing_domain)
            .evals;
        let folded_committed_values = stack_evaluations(evals, folding_factor);
        self.witness_coeff = folded_coeff;
        self.domain = scaled_domain;
        self.last_round_committed_values = self.committed_values.clone();
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
    pub fn update_challenges(&mut self) {
        self.challenges = dedup((0..self.num_repetitions[self.round_num]).map(|_| {
            squeeze_integer(
                &mut self.sponge,
                self.last_round_domain_size / self.folding_factor,
            )
        }));
        self.challenge_values = self
            .challenges
            .iter()
            .map(|index| self.last_round_committed_values[*index].clone())
            .collect();
        let mut challenge_answers: Vec<Path<M>> = Vec::with_capacity(self.challenge_values.len());
        for challenge in &self.challenges {
            challenge_answers.push(
                self.last_round_commitment
                    .generate_proof(*challenge)
                    .unwrap(),
            );
        }
        self.challenge_answers = challenge_answers;
    }
    pub fn update_commitment(&mut self) {
        self.last_round_commitment = self.commitment.clone();
        self.commitment = MerkleTree::<M>::new(
            &self.merkle_leaf_hash_param,
            &self.merkle_two_to_one_param,
            &self.committed_values,
        )
        .unwrap();
        // put it in the sponge
        self.sponge_absorb(&self.commitment.root());
    }
    pub fn update_folding_randomness(&mut self) {
        self.folding_randomness = self.sponge_squeeze();
    }
    pub fn update_out_of_domain_samples(&mut self) {
        self.out_of_domain_samples = self.sponge_squeeze_multiple(self.num_out_of_domain_samples);
        self.out_of_domain_evaluations = self
            .out_of_domain_samples
            .iter()
            .map(|point| self.witness_coeff.evaluate(point))
            .collect();
        self.sponge_absorb(&self.out_of_domain_evaluations.clone());
    }
    pub fn update_proof_of_work(&mut self) {
        self.proof_of_work_nonce = proof_of_work(
            &mut self.sponge,
            self.num_proof_of_work_bits[self.round_num],
        );
    }
    pub fn update_proximity_generator_randomness(&mut self) {
        self.proximity_generator_randomness = self.sponge_squeeze();
    }
    pub fn update_quotient_answers(&mut self) {
        let stir_randomness: Vec<F> = self
            .challenges
            .iter()
            .map(|index| self.domain.scale(self.folding_factor).element(*index))
            .collect();
        self.quotient_set = self
            .out_of_domain_samples
            .clone()
            .into_iter()
            .chain(stir_randomness.iter().cloned())
            .collect();
        self.quotient_answers = self
            .quotient_set
            .iter()
            .map(|x| self.witness_coeff.evaluate(x))
            .collect();
    }
}
