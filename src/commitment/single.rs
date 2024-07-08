use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, MerkleTree, Path, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;
use ark_poly::univariate::DensePolynomial;

use crate::{
    commitment::Witness,
    domain::Domain,
    utils::{squeeze_integer, stack_evaluations},
};

// pub struct SingleCommitment<M: MerkleConfig> {
//     commitment_digest: M::InnerDigest,
// }
// impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>> Commitment<F>
//     for SingleCommitment<M>
// where
//     M::InnerDigest: Absorb,
// {
//     type MerkleConfig = M;
//     fn commitment_digest(&self) -> <<Self as Commitment<F>>::MerkleConfig as MerkleConfig>::InnerDigest {
//         self.commitment_digest.clone()
//     }
//     fn new(argument: <Self::MerkleConfig as MerkleConfig>::InnerDigest) -> Self {
//         Self {
//             commitment_digest: argument,
//         }
//     }
// }

pub struct SingleWitness<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    argument: SingleWitnessArgument<F, M, S>,
    coeff: DensePolynomial<F>,
    commitment: MerkleTree<M>,
    committed_values: Vec<Vec<F>>,
}
impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> Witness<F>
    for SingleWitness<F, M, S>
where
    M::InnerDigest: Absorb,
    S::Config: Clone,
{
    type Argument = SingleWitnessArgument<F, M, S>;
    type Commitment = MerkleTree<M>;
    type CommittedValues = Vec<Vec<F>>;
    type Challenges = Vec<usize>;
    type ChallengeAnswers = Vec<Path<Self::MerkleConfig>>;
    type MerkleConfig = M;

    fn new(argument: Self::Argument) -> Self {
        // commit to the witness
        let evals: Vec<F> = argument
            .coeff
            .evaluate_over_domain_by_ref(argument.domain.backing_domain)
            .evals;
        let committed_values = stack_evaluations(evals, argument.folding_factor);
        let commitment = MerkleTree::<M>::new(
            &argument.merkle_leaf_hash_param,
            &argument.merkle_two_to_one_param,
            &committed_values,
        )
        .unwrap();
        Self {
            argument: argument.clone(),
            coeff: argument.coeff,
            commitment,
            committed_values,
        }
    }
    fn coeff(&self) -> DensePolynomial<F> {
        self.coeff.clone()
    }
    fn commitment_digest(
        &self,
    ) -> <<Self as Witness<F>>::MerkleConfig as MerkleConfig>::InnerDigest {
        self.commitment.root()
    }
    fn commitment(&self) -> MerkleTree<Self::MerkleConfig> {
        self.commitment.clone()
    }
    fn committed_values(&self) -> Self::CommittedValues {
        self.committed_values.clone()
    }
    fn challenges(&self, num_challenges: usize) -> Self::Challenges {
        // absorb committment digest
        let mut sponge = S::new(&self.argument.sponge_config);
        sponge.absorb(&self.commitment.root());
        // squeeze out the challenges as indices
        let mut challenges: Self::Challenges = Vec::with_capacity(num_challenges);
        for _ in 0..num_challenges {
            challenges.push(squeeze_integer(&mut sponge, 32)); // TODO (z-tech): this range must be set properly
        }
        challenges
    }
    fn challenge_answers(&self, challenges: Self::Challenges) -> Self::ChallengeAnswers {
        let mut challenge_answers: Vec<Path<M>> = Vec::with_capacity(challenges.len());
        for challenge in challenges {
            challenge_answers.push(self.commitment.generate_proof(challenge).unwrap());
        }
        challenge_answers
    }
    fn verify(
        &self,
        challenges: Self::Challenges,
        challenge_answers: Self::ChallengeAnswers,
    ) -> bool {
        for (&challenge, answer) in challenges.iter().zip(challenge_answers) {
            // the answer given should correspond to the correct challenge
            if !answer.leaf_index == challenge {
                return false;
            }

            // the proof should be valid with the given value against the digest
            if !answer
                .verify(
                    &self.argument.merkle_leaf_hash_param,
                    &self.argument.merkle_two_to_one_param,
                    &self.commitment.root(),
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
impl<F: FftField, M: MerkleConfig, S: CryptographicSponge> Clone for SingleWitness<F, M, S>
where
    S::Config: Clone,
{
    fn clone(&self) -> Self {
        Self {
            argument: self.argument.clone(),
            coeff: self.coeff.clone(),
            commitment: self.commitment.clone(),
            committed_values: self.committed_values.clone(),
        }
    }
}

// Use this to instantiate a SingleWitness
pub struct SingleWitnessArgument<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    pub coeff: DensePolynomial<F>,
    pub domain: Domain<F>,
    pub folding_factor: usize,
    pub merkle_leaf_hash_param: LeafParam<M>,
    pub merkle_two_to_one_param: TwoToOneParam<M>,
    pub sponge_config: S::Config,
}

impl<F: FftField, M: MerkleConfig, S: CryptographicSponge> Clone for SingleWitnessArgument<F, M, S>
where
    S::Config: Clone,
{
    fn clone(&self) -> Self {
        Self {
            coeff: self.coeff.clone(),
            domain: self.domain.clone(),
            folding_factor: self.folding_factor,
            merkle_leaf_hash_param: self.merkle_leaf_hash_param.clone(),
            merkle_two_to_one_param: self.merkle_two_to_one_param.clone(),
            sponge_config: self.sponge_config.clone(),
        }
    }
}
