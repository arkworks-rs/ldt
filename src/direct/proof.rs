use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, Path, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;

use crate::utils::squeeze_integer;

pub struct DirectProof<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    challenge_answers: Vec<Path<M>>,
    committed_values: Vec<Vec<F>>,
    merkle_leaf_hash_param: LeafParam<M>,
    merkle_two_to_one_param: TwoToOneParam<M>,
    sponge_config: S::Config,
}

impl<F, M, S> DirectProof<F, M, S>
where
    F: FftField,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
{
    pub fn new(
        challenges: Vec<usize>,
        challenge_answers: Vec<Path<M>>,
        committed_values: Vec<Vec<F>>,
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        sponge_config: S::Config,
    ) -> Self {
        Self {
            challenge_answers,
            committed_values,
            merkle_leaf_hash_param,
            merkle_two_to_one_param,
            sponge_config,
        }
    }
    pub fn challenges(
        &self,
        commitment_digest: M::InnerDigest,
        num_challenges: usize,
    ) -> Vec<usize> {
        // absorb committment digest
        let mut sponge = S::new(&self.sponge_config);
        sponge.absorb(&commitment_digest);
        // squeeze out the challenges as indices
        let mut challenges = Vec::with_capacity(num_challenges);
        for _ in 0..num_challenges {
            challenges.push(squeeze_integer(&mut sponge, 32)); // TODO (z-tech): this range must be set properly
        }
        challenges
    }
    pub fn verify(&self, commitment_digest: M::InnerDigest, challenges: Vec<usize>) -> bool {
        for (&challenge, answer) in challenges.iter().zip(self.challenge_answers.clone()) {
            // the answer given should correspond to the correct challenge
            if !answer.leaf_index == challenge {
                return false;
            }

            // the proof should be valid with the given value against the digest
            if !answer
                .verify(
                    &self.merkle_leaf_hash_param,
                    &self.merkle_two_to_one_param,
                    &commitment_digest,
                    self.committed_values[challenge].clone(),
                )
                .unwrap()
            {
                return false;
            }
        }
        true
    }
}
