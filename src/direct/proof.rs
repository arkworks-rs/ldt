use ark_crypto_primitives::merkle_tree::{Config as MerkleConfig, Path};
use ark_ff::Field;

pub struct DirectProof<F: Field, M: MerkleConfig> {
    pub p_commitment_root: M::InnerDigest,
    pub p_evaluations: Vec<Vec<F>>,
    pub inclusion_proofs: Vec<Path<M>>,
}
