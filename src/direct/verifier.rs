use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    argument::generate_challenges,
    direct::{config::DirectConfig, proof::DirectProof},
    ldt::Verifier,
};

pub struct DirectVerifier<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    config: DirectConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> Verifier<F>
    for DirectVerifier<F, M, S>
where
    M::InnerDigest: Absorb,
{
    type Config = DirectConfig<M, S>;
    type Proof = DirectProof<F, M>;

    fn new(config: DirectConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn verify(&self, proof: &Self::Proof) -> bool {
        // absorb the digest to derive the challenges
        let challenges = generate_challenges::<F, S>(
            proof.commitment_digest.clone(),
            self.config.num_challenges,
            &self.config.sponge_config,
        );

        // then verify each challenge
        for (&challenge, answer) in challenges.iter().zip(proof.challenge_answers.clone()) {
            // the answer given should correspond to the correct challenge
            if !answer.leaf_index == challenge {
                return false;
            }

            // the proof should be valid with the given value against the digest
            if !answer
                .verify(
                    &self.config.merkle_leaf_hash_param,
                    &self.config.merkle_two_to_one_param,
                    &proof.commitment_digest,
                    proof.committed_values[challenge].clone(),
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
