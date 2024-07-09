use ark_crypto_primitives::merkle_tree::{Config as MerkleConfig, Path};
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

pub struct STIRInnerRoundProof<F: Field, M: MerkleConfig> {
    pub answer_coeff: DensePolynomial<F>,
    pub challenge_answers: Vec<Path<M>>,
    pub committed_values: Vec<Vec<F>>,
    pub out_of_domain_evaluations: Vec<F>,
    pub commitment_digest: M::InnerDigest,
    pub proof_of_work_nonce: Option<usize>,
    pub shake_coeff: DensePolynomial<F>,
}

pub struct STIRFinalRoundProof<F: Field, M: MerkleConfig> {
    pub challenge_answers: Vec<Path<M>>,
    pub committed_values: Vec<Vec<F>>,
    pub coeff: DensePolynomial<F>,
    pub proof_of_work_nonce: Option<usize>,
}

pub struct STIRProof<F: Field, M: MerkleConfig> {
    pub final_round_proof: STIRFinalRoundProof<F, M>,
    pub initial_commitment_digest: M::InnerDigest,
    pub inner_round_proofs: Vec<STIRInnerRoundProof<F, M>>,
}
