use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MultiPath},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    direct::{config::DirectConfig, prover::DirectProver, verifier::DirectVerifier},
    ldt::{LowDegreeTest, Prover, Verifier},
    witness::Witness,
};

use super::proof::DirectProof;

pub struct DirectLDT<F: FftField, M: MerkleConfig, S: CryptographicSponge, W: Witness<F, M>> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge: PhantomData<S>,
    _witness: PhantomData<W>,
}
impl<F, M, S, W> LowDegreeTest<F> for DirectLDT<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<
            F,
            M,
            MerkleConfig = M,
            ChallengeAnswers = MultiPath<M>,
            CommittedValues = Vec<Vec<F>>,
            Challenges = Vec<usize>,
        > + Clone,
    W::ChallengeAnswers: Clone,
{
    type LDTConfig = DirectConfig<M, S>;
    type Proof = DirectProof<F, M, S>;
    type Prover = DirectProver<F, M, S, W>;
    type Verifier = DirectVerifier<F, M, S, W>;

    fn new(config: Self::LDTConfig) -> (Self::Prover, Self::Verifier) {
        (
            Self::Prover::new(config.clone()),
            Self::Verifier::new(config),
        )
    }
}
