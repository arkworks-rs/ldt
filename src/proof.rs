use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;

use crate::witness::Witness;

pub trait Proof<F: FftField, M: MerkleConfig, S: CryptographicSponge, W: Witness<F, M>> {
    fn new(
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        num_challenges: usize,
        sponge_config: <S as CryptographicSponge>::Config,
        starting_degree: usize,
        starting_rate: usize,
        witness: W,
    ) -> Self;
    fn verify(&self) -> bool;
}

pub struct SingleProof<F: FftField, M: MerkleConfig, S: CryptographicSponge, W: Witness<F, M>> {
    pub commitment_digest: M::InnerDigest,
    pub committed_values: W::CommittedValues,
    pub challenge_answers: W::ChallengeAnswers,
    pub merkle_leaf_hash_param: LeafParam<M>,
    pub merkle_two_to_one_param: TwoToOneParam<M>,
    pub num_challenges: usize,
    pub sponge_config: S::Config,
    pub witness: W,
}

impl<F: FftField, M: MerkleConfig, S: CryptographicSponge, W: Witness<F, M>> Proof<F, M, S, W>
    for SingleProof<F, M, S, W>
where
    M::InnerDigest: Absorb,
    <W as Witness<F, M>>::ChallengeAnswers: Clone,
{
    fn new(
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        num_challenges: usize,
        sponge_config: <S as CryptographicSponge>::Config,
        _starting_degree: usize,
        _starting_rate: usize,
        witness: W,
    ) -> Self {
        let challenges = witness.challenges(num_challenges);
        let challenge_answers = witness.challenge_answers(challenges);

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
