use ark_crypto_primitives::merkle_tree::{Config as MerkleConfig, Path};

pub struct DirectProof<M: MerkleConfig> {
    pub inclusion_proofs: Vec<Path<M>>,
}
