use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    direct::{config::DirectConfig, proof::DirectProof},
    ldt::Verifier,
    witness::Witness,
};

pub struct DirectVerifier<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
    W: Witness<F, M>,
{
    config: DirectConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
    _witness: PhantomData<W>,
}
impl<F, M, S, W> Verifier<F> for DirectVerifier<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    W: Witness<F, M>,
    W::ChallengeAnswers: Clone,
{
    type VerifierConfig = DirectConfig<M, S>;
    type Proof = DirectProof<F, M, S>;

    fn new(config: DirectConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
            _witness: PhantomData::<W>,
        }
    }
    fn verify(&self, proof: &Self::Proof) -> bool {
        // proof.verify() HERE pass the thing to verify the proof against
        true
    }
}
