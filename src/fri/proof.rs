use ark_crypto_primitives::merkle_tree::{Config as MerkleConfig, Path};
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

pub struct FRIRoundProof<F: Field, M: MerkleConfig> {
    pub queries_to_prev: (Vec<Vec<F>>, Vec<Path<M>>),
}

pub struct FRIProof<F: Field, M: MerkleConfig> {
    pub commitment_digests: Vec<M::InnerDigest>,
    pub initial_commitment_digest: M::InnerDigest,
    pub coeff: DensePolynomial<F>,
    pub round_proofs: Vec<FRIRoundProof<F, M>>,
    pub proof_of_work_nonce: Option<usize>,
}
