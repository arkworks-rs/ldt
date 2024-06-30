use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree, Path},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;

use crate::utils::squeeze_integer;

pub fn generate_challenges<F: FftField, S: CryptographicSponge>(
    commitment_digest: impl Absorb,
    num_challenges: usize,
    sponge_config: &S::Config,
) -> Vec<usize> {
    // absorb committment digest
    let mut sponge = S::new(&sponge_config);
    sponge.absorb(&commitment_digest);
    // squeeze out the challenges as indices
    let mut challenges: Vec<usize> = Vec::with_capacity(num_challenges);
    for _ in 0..num_challenges {
        challenges.push(squeeze_integer(&mut sponge, 32)); // TODO (z-tech): this range must be set properly
    }
    // return as vec of usizes
    challenges
}

pub trait Argument<F: FftField, M: MerkleConfig> {
    fn generate_challenges(&self) -> Vec<usize>;
    fn generate_challenge_answers(&self, challenges: Vec<usize>) -> Vec<Path<M>>;
}

pub struct SingleArgument<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    pub commitment: MerkleTree<M>,
    pub committed_values: Vec<Vec<F>>,
    pub num_challenges: usize,
    pub sponge_config: S::Config,
}

impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> Argument<F, M>
    for SingleArgument<F, M, S>
where
    M::InnerDigest: Absorb,
{
    fn generate_challenges(&self) -> Vec<usize> {
        generate_challenges::<F, S>(
            self.commitment.root(),
            self.num_challenges,
            &self.sponge_config,
        )
    }
    fn generate_challenge_answers(&self, challenges: Vec<usize>) -> Vec<Path<M>> {
        let mut challenge_answers: Vec<Path<M>> = Vec::with_capacity(challenges.len());
        for challenge in challenges {
            challenge_answers.push(self.commitment.generate_proof(challenge).unwrap());
        }
        challenge_answers
    }
}
