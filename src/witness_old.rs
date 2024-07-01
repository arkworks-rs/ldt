use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, MerkleTree, Path, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;

pub trait Witness<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    fn commitment(&self) -> MerkleTree<M>;
    fn commitment_digest(&self) -> M::InnerDigest;
    fn committed_values(&self) -> Vec<Vec<F>>;
    fn generate_challenges(&self) -> Vec<usize>;
    fn generate_challenge_answers(&self, challenges: Vec<usize>) -> Vec<Path<M>>;
}

pub struct SingleArgument<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    commitment: MerkleTree<M>,
    committed_values: Vec<Vec<F>>,
    num_challenges: usize,
    sponge_config: S::Config,
}

impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> Argument<F, M, S>
    for SingleArgument<F, M, S>
where
    M::InnerDigest: Absorb,
{
    fn new(
        committed_values: Vec<Vec<F>>,
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        num_challenges: usize,
        sponge_config: S::Config,
    ) -> Self {
        // commit to the values
        let commitment = MerkleTree::<M>::new(
            &merkle_leaf_hash_param,
            &merkle_two_to_one_param,
            &committed_values,
        )
        .unwrap();
        // return self
        Self {
            commitment,
            committed_values,
            num_challenges,
            sponge_config,
        }
    }
    fn commitment(&self) -> MerkleTree<M> {
        self.commitment.clone()
    }
    fn commitment_digest(&self) -> M::InnerDigest {
        self.commitment.root()
    }
    fn committed_values(&self) -> Vec<Vec<F>> {
        self.committed_values.clone()
    }
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
