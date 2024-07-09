use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_std::marker::PhantomData;

use crate::{
    commitment::Witness,
    ldt::{LowDegreeTest, Prover, Verifier},
    stir::{config::STIRConfig, proof::STIRProof, prover::STIRProver, verifier::STIRVerifier},
};

pub struct STIR<F: FftField + PrimeField, M: MerkleConfig, S: CryptographicSponge, W: Witness<F, M>>
{
    _field: PhantomData<F>,
    _merkle_config: PhantomData<W::MerkleConfig>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField + PrimeField, M: MerkleConfig, S: CryptographicSponge, W: Witness<F, M>>
    LowDegreeTest<F> for STIR<F, M, S, W>
where
    F: Absorb,
    W: Clone,
    W::ChallengeAnswers: Clone,
    W::MerkleConfig: MerkleConfig<Leaf = Vec<F>>,
    <W::MerkleConfig as MerkleConfig>::InnerDigest: Absorb,
    S::Config: Clone,
{
    type Config = STIRConfig<W::MerkleConfig, S>;
    type Proof = STIRProof<F, W::MerkleConfig>;
    type Prover = STIRProver<F, M, S, W>;
    type Verifier = STIRVerifier<F, M, S, W>;

    fn new(config: Self::Config) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
