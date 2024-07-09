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
impl<F: FftField, S: CryptographicSponge, M: MerkleConfig, W: Witness<F, M, MerkleConfig = M>>
    LowDegreeTest<F> for DirectLDT<F, M, S, W>
where
    S::Config: Clone,
    W: Clone,
    W::ChallengeAnswers: Clone,
    W::MerkleConfig: MerkleConfig,
    M::InnerDigest: Absorb,
{
    type Config = DirectConfig<W::MerkleConfig, S>;
    type Proof = SingleProof<F, M, S, W>;
    type Prover = DirectProver<F, M, S, W>;
    type Verifier = DirectVerifier<F, M, S, W>;

    fn new(config: Self::Config) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
