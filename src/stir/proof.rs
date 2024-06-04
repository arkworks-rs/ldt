use ark_crypto_primitives::merkle_tree::{Config as MerkleConfig, Path};
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

pub struct STIRRoundProof<F: Field, M: MerkleConfig> {
    pub ans_polynomial: DensePolynomial<F>,
    pub betas: Vec<F>,
    pub g_root: M::InnerDigest,
    pub proof_of_work_nonce: Option<usize>,
    pub queries_to_prev: (Vec<Vec<F>>, Vec<Path<M>>),
    pub shake_polynomial: DensePolynomial<F>,
}

pub struct STIRProof<F: Field, M: MerkleConfig> {
    pub polynomial: DensePolynomial<F>,
    pub round_proofs: Vec<STIRRoundProof<F, M>>,
    pub proof_of_work_nonce: Option<usize>,
    pub queries_to_final: (Vec<Vec<F>>, Vec<Path<M>>),
}
