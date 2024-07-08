use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_std::marker::PhantomData;

use crate::{
    commitment::Witness, fri::{config::FRIConfig, proof::FRIProof, prover::FRIProver, verifier::FRIVerifier}, ldt::{LowDegreeTest, Prover, Verifier}
};

pub struct FRI<F: FftField, S: CryptographicSponge, W: Witness<F>> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<W::MerkleConfig>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField + PrimeField, S: CryptographicSponge, W: Witness<F>>
    LowDegreeTest<F> for FRI<F, S, W>
where
    W: Clone,
    W::ChallengeAnswers: Clone,
    W::MerkleConfig: MerkleConfig<Leaf = Vec<F>>,
    <W::MerkleConfig as MerkleConfig>::InnerDigest: Absorb,
    S::Config: Clone,
{
    type Config = FRIConfig<W::MerkleConfig, S>;
    type Proof = FRIProof<F, W::MerkleConfig>;
    type Prover = FRIProver<F, S, W>;
    type Verifier = FRIVerifier<F, S, W>;

    fn new(config: Self::Config) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
