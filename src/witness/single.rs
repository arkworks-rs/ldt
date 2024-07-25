use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, MerkleTree, MultiPath, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;
use ark_poly::{univariate::DensePolynomial, Polynomial};

#[cfg(not(feature = "std"))]
use ark_std::vec::Vec;

use crate::{
    domain::Domain,
    statement::single::SingleStatement,
    utils::{squeeze_integer, stack_evaluations},
    witness::Witness,
};

pub struct SingleWitness<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    argument: SingleWitnessArgument<F, M, S>,
    coeff: DensePolynomial<F>,
    domain: Domain<F>,
    commitment: MerkleTree<M>,
    committed_values: Vec<Vec<F>>,
}

impl<F, M, S> Witness<F, M> for SingleWitness<F, M, S>
where
    F: FftField,
    M: MerkleConfig<Leaf = Vec<F>> + Clone,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
{
    type Argument = SingleWitnessArgument<F, M, S>;
    type Commitment = MerkleTree<M>;
    type CommittedValues = Vec<Vec<F>>;
    type Challenges = Vec<usize>;
    type ChallengeAnswers = MultiPath<M>;
    type Statement = SingleStatement<M>;
    type MerkleConfig = M;

    fn new(argument: Self::Argument) -> Self {
        // 1) Generate a commitment for the argument
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

        // 2) Keep everything the prover will need
        Self {
            argument: argument.clone(),
            coeff: argument.coeff,
            domain: argument.domain,
            commitment,
            committed_values,
        }
    }
    fn coeff(&self) -> DensePolynomial<F> {
        self.coeff.clone()
    }
    fn coeff_degree(&self) -> usize {
        self.coeff.degree()
    }
    fn commitment_digest(&self) -> M::InnerDigest {
        self.commitment.root()
    }
    fn commitment(&self) -> MerkleTree<Self::MerkleConfig> {
        self.commitment.clone()
    }
    fn committed_values(&self) -> Self::CommittedValues {
        self.committed_values.clone()
    }
    fn challenges(&self, num_challenges: usize) -> Vec<usize> {
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
        self.commitment.generate_multi_proof(challenges).unwrap()
    }
    fn domain(&self) -> Domain<F> {
        self.domain.clone()
    }
    fn statement(&self) -> Self::Statement {
        SingleStatement::<M>::new(self.commitment_digest())
    }
    fn verify(
        &self,
        challenges: Self::Challenges,
        challenge_answers: Self::ChallengeAnswers,
    ) -> bool {
        if challenge_answers.leaf_indexes != challenges {
            return false;
        }
        challenge_answers
            .verify(
                &self.argument.merkle_leaf_hash_param,
                &self.argument.merkle_two_to_one_param,
                &self.commitment.root(),
                self.committed_values.clone(),
            )
            .unwrap()
    }
}
impl<F, M, S> Clone for SingleWitness<F, M, S>
where
    F: FftField,
    M: MerkleConfig + Clone,
    S: CryptographicSponge,
    S::Config: Clone,
{
    fn clone(&self) -> Self {
        Self {
            argument: self.argument.clone(),
            coeff: self.coeff.clone(),
            domain: self.domain.clone(),
            commitment: self.commitment.clone(),
            committed_values: self.committed_values.clone(),
        }
    }
}

// Use SingleWitnessArguement to instantiate a SingleWitness
#[derive(Clone)]
pub struct SingleWitnessArgument<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    pub coeff: DensePolynomial<F>,
    pub domain: Domain<F>,
    pub folding_factor: usize,
    pub merkle_leaf_hash_param: LeafParam<M>,
    pub merkle_two_to_one_param: TwoToOneParam<M>,
    pub sponge_config: S::Config,
}
