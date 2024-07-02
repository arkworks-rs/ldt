use ark_crypto_primitives::{merkle_tree::Config as MerkleConfig, sponge::{Absorb, CryptographicSponge}};
use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    commitment::Witness,
    direct::{config::DirectConfig, prover::DirectProver, verifier::DirectVerifier},
    ldt::{LowDegreeTest, Prover, Verifier},
    proof::SingleProof,
};

pub struct DirectLDT<F: FftField, S: CryptographicSponge, W: Witness<F>> {
    _field: PhantomData<F>,
    _sponge_config: PhantomData<S>,
    _witness: PhantomData<W>,
}
impl<F: FftField, S: CryptographicSponge, W: Witness<F>> LowDegreeTest<F> for DirectLDT<F, S, W>
where
    S::Config: Clone,
    W: Clone,
    W::ChallengeAnswers: Clone,
    W::MerkleConfig: MerkleConfig,
    <W::MerkleConfig as MerkleConfig>::InnerDigest: Absorb,
{
    type Config = DirectConfig<W::MerkleConfig, S>;
    type Proof = SingleProof<F, S, W>;
    type Prover = DirectProver<F, S, W>;
    type Verifier = DirectVerifier<F, S, W>;

    fn new(config: Self::Config) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
