use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};

use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    direct::{
        config::DirectConfig, proof::DirectProof, prover::DirectProver, verifier::DirectVerifier,
    },
    ldt::{LowDegreeTest, Prover, Verifier},
};

pub struct DirectLDT<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> LowDegreeTest<F>
    for DirectLDT<F, M, S>
where
    M::InnerDigest: Absorb,
    S::Config: Clone,
{
    type Config = DirectConfig<M, S>;
    type Proof = DirectProof<M>;
    type Prover = DirectProver<F, M, S>;
    type Verifier = DirectVerifier<F, M, S>;

    fn new(config: Self::Config) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
