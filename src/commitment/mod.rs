use ark_crypto_primitives::merkle_tree::Config as MerkleConfig;
use ark_ff::FftField;
use ark_poly::univariate::DensePolynomial;

pub mod single;

// pub trait Commitment<F: FftField> where Self::MerkleConfig: MerkleConfig
// {
//     type MerkleConfig;
//     fn new(argument: <Self::MerkleConfig as MerkleConfig>::InnerDigest) -> Self;
//     fn commitment_digest(
//         &self,
//     ) -> <<Self as Commitment<F>>::MerkleConfig as MerkleConfig>::InnerDigest
//     where
//         <Self as Commitment<F>>::MerkleConfig: ark_crypto_primitives::merkle_tree::Config;
// }

pub trait Witness<F: FftField> {
    type Argument;
    type Commitment;
    type Challenges;
    type ChallengeAnswers;
    type CommittedValues;
    type MerkleConfig;

    fn new(argument: Self::Argument) -> Self;
    fn coeff(&self) -> DensePolynomial<F>;
    fn commitment(&self) -> Self::Commitment;
    fn commitment_digest(
        &self,
    ) -> <<Self as Witness<F>>::MerkleConfig as MerkleConfig>::InnerDigest
    where
        <Self as Witness<F>>::MerkleConfig: ark_crypto_primitives::merkle_tree::Config;
    fn committed_values(&self) -> Self::CommittedValues;
    fn challenges(&self, num_challenges: usize) -> Self::Challenges;
    fn challenge_answers(&self, challenges: Self::Challenges) -> Self::ChallengeAnswers;
    fn verify(
        &self,
        challenges: Self::Challenges,
        challenge_answers: Self::ChallengeAnswers,
    ) -> bool;
}
