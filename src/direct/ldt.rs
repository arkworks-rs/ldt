use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    commitment::Witness,
    direct::{config::DirectConfig, prover::DirectProver, verifier::DirectVerifier},
    ldt::{LowDegreeTest, Prover, Verifier},
    proof::SingleProof,
};

pub struct DirectLDT<F: FftField, M: MerkleConfig, S: CryptographicSponge, W: Witness<F, M>> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
    _witness: PhantomData<W>,
}
impl<F, M, S, W> LowDegreeTest<F> for DirectLDT<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M> + Clone,
    W::ChallengeAnswers: Clone,
{
    type LDTConfig = DirectConfig<M, S>;
    type Proof = SingleProof<F, M, S, W>;
    type Prover = DirectProver<F, M, S, W>;
    type Verifier = DirectVerifier<F, M, S, W>;

    fn new(config: Self::LDTConfig) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
