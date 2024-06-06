use ark_crypto_primitives::merkle_tree::{Config as MerkleConfig, Path};
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

pub struct STIRInnerRoundProof<F: Field, M: MerkleConfig> {
    pub answer_polynomial: DensePolynomial<F>,
    pub inclusion_proofs_of_queries: Vec<Path<M>>,
    pub leaf_values_of_queries: Vec<Vec<F>>,
    pub out_of_domain_evaluations: Vec<F>,
    pub p_commitment_root: M::InnerDigest,
    pub proof_of_work_nonce: Option<usize>,
    pub shake_polynomial: DensePolynomial<F>,
}

pub struct STIRFinalRoundProof<F: Field, M: MerkleConfig> {
    pub inclusion_proofs_of_queries: Vec<Path<M>>,
    pub leaf_values_of_queries: Vec<Vec<F>>,
    pub polynomial: DensePolynomial<F>,
    pub proof_of_work_nonce: Option<usize>,
}

pub struct STIRProof<F: Field, M: MerkleConfig> {
    pub final_round_proof: STIRFinalRoundProof<F, M>,
    pub inner_round_proofs: Vec<STIRInnerRoundProof<F, M>>,
}
