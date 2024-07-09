use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, Path, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;

pub struct DirectProof<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    challenges: Vec<usize>,
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
            challenges,
            challenge_answers,
            committed_values,
            merkle_leaf_hash_param,
            merkle_two_to_one_param,
            sponge_config,
        }
    }
    fn verify(&self, commitment_digest: M::InnerDigest) -> bool {
        for (&challenge, answer) in self.challenges.iter().zip(self.challenge_answers.clone()) {
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
