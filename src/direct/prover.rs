use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};

use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    direct::config::DirectConfig,
    ldt::Prover,
    proof::{Proof, SingleProof},
    witness::Witness,
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
    S::Config: Clone,
{
    type Config = DirectConfig<M, S>;
    type Proof = SingleProof<F, M, S>;

    fn new(config: DirectConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn prove(&self, witness: impl Witness<F>) -> Self::Proof {
        Self::Proof::new(
            self.config.merkle_leaf_hash_param.clone(),
            self.config.merkle_two_to_one_param.clone(),
            self.config.num_challenges,
            self.config.sponge_config.clone(),
            self.config.degree,
            1,
            witness,
        )
    }
}
