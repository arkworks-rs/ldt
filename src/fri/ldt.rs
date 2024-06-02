use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};

use ark_ff::{FftField, PrimeField};
use ark_std::marker::PhantomData;

use crate::{
    fri::{config::FRIConfig, proof::FRIProof, prover::FRIProver, verifier::FRIVerifier},
    ldt::{LowDegreeTest, Prover, Verifier},
};

pub struct FRI<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField + PrimeField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge>
    LowDegreeTest<F> for FRI<F, M, S>
where
    M::InnerDigest: Absorb,
    S::Config: Clone,
{
    type Config = FRIConfig<M, S>;
    type Proof = FRIProof<F, M>;
    type Prover = FRIProver<F, M, S>;
    type Verifier = FRIVerifier<F, M, S>;

    fn new(config: Self::Config) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
