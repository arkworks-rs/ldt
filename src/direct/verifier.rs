use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    commitment::Witness,
    direct::config::DirectConfig,
    ldt::Verifier,
    proof::{Proof, SingleProof},
};

pub struct DirectVerifier<F: FftField, M: MerkleConfig, S: CryptographicSponge, W: Witness<F, M>> {
    config: DirectConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<W::MerkleConfig>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField, M: MerkleConfig, S: CryptographicSponge, W: Witness<F, M>> Verifier<F>
    for DirectVerifier<F, M, S, W>
where
    M::InnerDigest: Absorb,
    <W as Witness<F, M>>::ChallengeAnswers: Clone,
{
    type Config = DirectConfig<M, S>;
    type Proof = SingleProof<F, M, S, W>;

    fn new(config: DirectConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<W::MerkleConfig>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn verify(&self, proof: &Self::Proof) -> bool {
        proof.verify()
    }
}
