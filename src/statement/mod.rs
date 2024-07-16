use ark_crypto_primitives::{merkle_tree::Config as MerkleConfig, sponge::CryptographicSponge};
use ark_ff::FftField;

pub mod single;

pub trait Statement<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    type Argument;
    fn new(argument: Self::Argument) -> Self;
    fn commitment_digest(&self) -> M::InnerDigest;
}
