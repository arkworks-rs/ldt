use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, MultiPath, TwoToOneParam},
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
    challenge_answers: MultiPath<M>,
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
        challenge_answers: MultiPath<M>,
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
            challenges.push(squeeze_integer(&mut sponge, 32));
        }
        challenges
    }
    pub fn verify(&self, commitment_digest: M::InnerDigest, challenges: Vec<usize>) -> bool {
        if self.challenge_answers.leaf_indexes
            != challenges.iter().rev().cloned().collect::<Vec<usize>>()
        {
            // TODO: IDK why self.challenge_answers.leaf_indexes comes back in reverse
            return false;
        }

        let challenge_values: Vec<Vec<F>> = self
            .challenge_answers
            .leaf_indexes
            .iter()
            .map(|&i| self.committed_values.get(i).unwrap().clone())
            .collect();
        self.challenge_answers
            .verify(
                &self.merkle_leaf_hash_param,
                &self.merkle_two_to_one_param,
                &commitment_digest,
                challenge_values,
            )
            .unwrap()
    }
}
