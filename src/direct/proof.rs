use ark_crypto_primitives::merkle_tree::{Config as MerkleConfig, Path};
use ark_ff::Field;

pub struct DirectProof<F: Field, M: MerkleConfig> {
    pub commitment_digest: M::InnerDigest,
    pub committed_values: Vec<Vec<F>>,
    pub challenge_answers: Vec<Path<M>>,
}
