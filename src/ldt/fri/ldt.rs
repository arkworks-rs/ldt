use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_std::marker::PhantomData;

use crate::{
    ldt::{
        fri::{config::FRIConfig, proof::FRIProof, prover::FRIProver, verifier::FRIVerifier},
        LowDegreeTest, Prover, Verifier,
    },
    witness::Witness,
};

pub struct FRI<F, M, S, W>
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

impl<F, M, S, W> LowDegreeTest<F> for FRI<F, M, S, W>
where
    F: FftField + PrimeField,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    W: Witness<F, M, Commitment = MerkleTree<M>, MerkleConfig = M, CommittedValues = Vec<Vec<F>>>,
    W: Clone,
    W::ChallengeAnswers: Clone,
    S: CryptographicSponge,
    S::Config: Clone,
{
    type LDTConfig = FRIConfig<M, S>;
    type Proof = FRIProof<F, M>;
    type Prover = FRIProver<F, M, S, W>;
    type Verifier = FRIVerifier<F, M, S, W>;

    fn new(config: Self::LDTConfig) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
