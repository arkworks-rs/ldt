use ark_crypto_primitives::merkle_tree::{Config as MerkleConfig, Path};
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

pub struct FRIRoundProof<F: Field, M: MerkleConfig> {
    pub queries_to_prev: (Vec<Vec<F>>, Vec<Path<M>>),
}

pub struct FRIProof<F: Field, M: MerkleConfig> {
    pub commitments: Vec<<M>::InnerDigest>,
    pub polynomial: DensePolynomial<F>,
    pub round_proofs: Vec<FRIRoundProof<F, M>>,
    pub proof_of_work_nonce: Option<usize>,
}
