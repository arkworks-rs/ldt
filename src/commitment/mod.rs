use ark_crypto_primitives::merkle_tree::Config as MerkleConfig;
use ark_ff::FftField;
use ark_poly::univariate::DensePolynomial;

use crate::domain::Domain;

pub mod single;

pub trait Witness<F: FftField, M: MerkleConfig> {
    type Argument;
    type Commitment;
    type Challenges;
    type ChallengeAnswers;
    type CommittedValues;
    type MerkleConfig;

    fn new(argument: Self::Argument) -> Self;
    fn coeff(&self) -> DensePolynomial<F>;
    fn commitment(&self) -> Self::Commitment;
    fn commitment_digest(&self) -> M::InnerDigest;
    fn committed_values(&self) -> Self::CommittedValues;
    fn challenges(&self, num_challenges: usize) -> Self::Challenges;
    fn challenge_answers(&self, challenges: Self::Challenges) -> Self::ChallengeAnswers;
    fn domain(&self) -> Domain<F>;
    fn verify(
        &self,
        challenges: Self::Challenges,
        challenge_answers: Self::ChallengeAnswers,
    ) -> bool;
}
