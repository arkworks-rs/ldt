use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, Path},
    sponge::{Absorb, CryptographicSponge},
};

use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    direct::{config::DirectConfig, proof::DirectProof},
    ldt::Prover,
    witness::Witness,
};

pub struct DirectProver<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
    W: Witness<F, M>,
{
    config: DirectConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge: PhantomData<S>,
    _witness: PhantomData<W>,
}

impl<F, M, S, W> Prover<F> for DirectProver<F, M, S, W>
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
            CommittedValues = Vec<Vec<F>>,
            ChallengeAnswers = Vec<Path<M>>,
            Challenges = Vec<usize>,
        > + Clone,
    W::ChallengeAnswers: Clone,
{
    type Witness = W;
    type ProverConfig = DirectConfig<M, S>;
    type Proof = DirectProof<F, M, S>;

    fn new(config: DirectConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge: PhantomData::<S>,
            _witness: PhantomData::<W>,
        }
    }

    fn prove(&self, witness: &W) -> Self::Proof {
        let challenges = witness.challenges(self.config.num_challenges);
        DirectProof::<F, M, S>::new(
            witness.challenge_answers(challenges),
            witness.committed_values(),
            self.config.merkle_leaf_hash_param.clone(),
            self.config.merkle_two_to_one_param.clone(),
            self.config.sponge_config.clone(),
        )
    }
}
