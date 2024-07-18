use crate::{
    ldt::{LowDegreeTest, Prover, Verifier},
    stir::{config::STIRConfig, proof::STIRProof, prover::STIRProver, verifier::STIRVerifier},
    witness::Witness,
};
use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_std::marker::PhantomData;

pub struct STIR<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
    W: Witness<F, M>,
{
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge: PhantomData<S>,
    _witness: PhantomData<W>,
}

impl<F, M, S, W> LowDegreeTest<F> for STIR<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>> + Clone,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M, Commitment = MerkleTree<M>, CommittedValues = Vec<Vec<F>>>
        + Clone,
    W::ChallengeAnswers: Clone,
{
    type LDTConfig = STIRConfig<M, S>;
    type Proof = STIRProof<F, M, S>;
    type Prover = STIRProver<F, M, S, W>;
    type Verifier = STIRVerifier<F, M, S, W>;

    fn new(config: Self::LDTConfig) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
