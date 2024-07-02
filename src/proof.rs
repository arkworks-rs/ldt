use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;

use crate::commitment::Witness;

pub trait Proof<F: FftField, S: CryptographicSponge, W: Witness<F>>
where
    <W as Witness<F>>::MerkleConfig: ark_crypto_primitives::merkle_tree::Config,
{
    fn new(
        merkle_leaf_hash_param: LeafParam<W::MerkleConfig>,
        merkle_two_to_one_param: TwoToOneParam<W::MerkleConfig>,
        num_challenges: usize,
        sponge_config: <S as CryptographicSponge>::Config,
        starting_degree: usize,
        starting_rate: usize,
        witness: W,
    ) -> Self;
    fn verify(&self) -> bool;
}

pub struct SingleProof<F: FftField, S: CryptographicSponge, W: Witness<F>>
where
    <W as Witness<F>>::MerkleConfig: ark_crypto_primitives::merkle_tree::Config,
{
    pub commitment_digest: <<W as Witness<F>>::MerkleConfig as MerkleConfig>::InnerDigest,
    pub committed_values: W::CommittedValues,
    pub challenge_answers: W::ChallengeAnswers,
    pub merkle_leaf_hash_param: LeafParam<W::MerkleConfig>,
    pub merkle_two_to_one_param: TwoToOneParam<W::MerkleConfig>,
    pub num_challenges: usize,
    pub sponge_config: S::Config,
    pub witness: W,
}

impl<F: FftField, S: CryptographicSponge, W: Witness<F>> Proof<F, S, W> for SingleProof<F, S, W>
where
    <W as Witness<F>>::MerkleConfig: ark_crypto_primitives::merkle_tree::Config,
    <<W as Witness<F>>::MerkleConfig as ark_crypto_primitives::merkle_tree::Config>::InnerDigest:
        Absorb,
    <W as Witness<F>>::ChallengeAnswers: Clone,
{
    fn new(
        merkle_leaf_hash_param: LeafParam<W::MerkleConfig>,
        merkle_two_to_one_param: TwoToOneParam<W::MerkleConfig>,
        num_challenges: usize,
        sponge_config: <S as CryptographicSponge>::Config,
        starting_degree: usize,
        starting_rate: usize,
        witness: W,
    ) -> Self {
        let challenges = witness.challenges(num_challenges);
        let mut challenge_answers = witness.challenge_answers(challenges);

        Self {
            commitment_digest: witness.commitment_digest(),
            committed_values: witness.committed_values(),
            challenge_answers,
            merkle_leaf_hash_param,
            merkle_two_to_one_param,
            num_challenges,
            sponge_config,
            witness,
        }
    }
    fn verify(&self) -> bool {
        let challenges = self.witness.challenges(self.num_challenges);
        self.witness
            .verify(challenges, self.challenge_answers.clone())
    }
}
