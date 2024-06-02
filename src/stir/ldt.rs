use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_std::marker::PhantomData;

use crate::{
    ldt::{LowDegreeTest, Prover, Verifier},
    stir::{config::STIRConfig, proof::STIRProof, prover::STIRProver, verifier::STIRVerifier},
};

pub struct STIR<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField + PrimeField + Absorb, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge>
    LowDegreeTest<F> for STIR<F, M, S>
where
    M::InnerDigest: Absorb,
    S::Config: Clone,
{
    type Config = STIRConfig<M, S>;
    type Proof = STIRProof<F, M>;
    type Prover = STIRProver<F, M, S>;
    type Verifier = STIRVerifier<F, M, S>;

    fn new(config: Self::Config) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
