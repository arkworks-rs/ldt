use ark_crypto_primitives::merkle_tree::Config as MerkleConfig;

pub mod single;

pub trait Statement<M: MerkleConfig> {
    type Argument;
    fn new(argument: Self::Argument) -> Self;
    fn commitment_digest(&self) -> M::InnerDigest;
}
