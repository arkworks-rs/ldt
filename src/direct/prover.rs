use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};

use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    commitment::Witness,
    direct::config::DirectConfig,
    ldt::Prover,
    proof::{Proof, SingleProof},
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
    _merkle_config: PhantomData<W::MerkleConfig>,
    _sponge_config: PhantomData<S>,
}

impl<F, M, S, W> Prover<F> for DirectProver<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M> + Clone,
    W::ChallengeAnswers: Clone,
{
    type Witness = W;
    type Config = DirectConfig<W::MerkleConfig, S>;
    type Proof = SingleProof<F, M, S, W>;

    fn new(config: DirectConfig<W::MerkleConfig, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<W::MerkleConfig>,
            _sponge_config: PhantomData::<S>,
        }
    }

    fn prove(&self, witness: &W) -> Self::Proof {
        SingleProof::<F, M, S, W>::new(
            self.config.merkle_leaf_hash_param.clone(),
            self.config.merkle_two_to_one_param.clone(),
            self.config.num_challenges,
            self.config.sponge_config.clone(),
            self.config.degree,
            1,
            witness.clone(),
        )
    }
}
