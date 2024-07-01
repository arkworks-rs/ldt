use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, MerkleTree, Path, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;

use crate::{
    domain::Domain,
    utils::squeeze_integer,
    witness::Witness,
};

pub fn generate_challenges<F: FftField, S: CryptographicSponge>(
    commitment_digest: impl Absorb,
    num_challenges: usize,
    sponge_config: &S::Config,
) -> Vec<usize> {
    // absorb committment digest
    let mut sponge = S::new(&sponge_config);
    sponge.absorb(&commitment_digest);
    // squeeze out the challenges as indices
    let mut challenges: Vec<usize> = Vec::with_capacity(num_challenges);
    for _ in 0..num_challenges {
        challenges.push(squeeze_integer(&mut sponge, 32)); // TODO (z-tech): this range must be set properly
    }
    // return as vec of usizes
    challenges
}

pub fn generate_challenge_answers<F: FftField, M: MerkleConfig>(
    commitment: MerkleTree<M>,
    challenges: Vec<usize>,
) -> Vec<Path<M>> {
    let mut challenge_answers: Vec<Path<M>> = Vec::with_capacity(challenges.len());
    for challenge in challenges {
        challenge_answers.push(commitment.generate_proof(challenge).unwrap());
    }
    challenge_answers
}

pub trait Proof<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    fn new(
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        num_challenges: usize,
        sponge_config: <S as CryptographicSponge>::Config,
        starting_degree: usize,
        starting_rate: usize,
        witness: impl Witness<F>,
    ) -> Self;
    fn verify(&self) -> bool;
}

pub struct SingleProof<F: FftField, M: MerkleConfig, S: CryptographicSponge>
where
    M::InnerDigest: Absorb,
{
    commitment_digest: M::InnerDigest,
    committed_values: Vec<Vec<F>>,
    challenge_answers: Vec<Path<M>>,
    merkle_leaf_hash_param: LeafParam<M>,
    merkle_two_to_one_param: TwoToOneParam<M>,
    num_challenges: usize,
    sponge_config: S::Config,
}

impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> Proof<F, M, S>
    for SingleProof<F, M, S>
where
    M::InnerDigest: Absorb,
{
    fn new(
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        num_challenges: usize,
        sponge_config: <S as CryptographicSponge>::Config,
        starting_degree: usize,
        starting_rate: usize,
        witness: impl Witness<F>,
    ) -> Self {
        // commit to the witness
        let domain = Domain::<F>::new(starting_degree, starting_rate).unwrap();
        let committed_values = witness.folded_evaluations_over_domain(domain, 1);
        let commitment = MerkleTree::<M>::new(
            &merkle_leaf_hash_param,
            &merkle_two_to_one_param,
            &committed_values,
        )
        .unwrap();

        // generate challenges
        let challenges =
            generate_challenges::<F, S>(commitment.root(), num_challenges, &sponge_config);
        
        // generate challenge answers
        let challenge_answers = generate_challenge_answers::<F, M>(commitment.clone(), challenges);

        Self {
            commitment_digest: commitment.root(),
            committed_values,
            challenge_answers,
            merkle_leaf_hash_param,
            merkle_two_to_one_param,
            num_challenges,
            sponge_config,
        }
    }
    fn verify(&self) -> bool {
        // absorb the digest to derive the challenges
        let challenges = generate_challenges::<F, S>(
            self.commitment_digest.clone(),
            self.num_challenges,
            &self.sponge_config,
        );

        // then verify each challenge
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
                    &self.commitment_digest,
                    self.committed_values[challenge].clone(),
                )
                .unwrap()
            {
                return false;
            }
        }

        // verification is accepted
        true
    }
}
