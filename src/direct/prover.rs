use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};

use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    argument::{Argument, SingleArgument},
    direct::{config::DirectConfig, proof::DirectProof},
    ldt::Prover,
};

pub struct DirectProver<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    config: DirectConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> Prover<F>
    for DirectProver<F, M, S>
where
    M::InnerDigest: Absorb,
{
    type Argument = SingleArgument<F, M, S>;
    type Config = DirectConfig<M, S>;
    type Proof = DirectProof<F, M>;

    fn new(config: DirectConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn prove(&self, argument: &Self::Argument) -> Self::Proof {
        let challenges = argument.generate_challenges();
        let challenge_answers = argument.generate_challenge_answers(challenges);
        Self::Proof {
            commitment_digest: argument.commitment.root(),
            committed_values: argument.committed_values.clone(),
            challenge_answers: challenge_answers,
        }
    }
}
