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

pub struct DirectProver<F: FftField, S: CryptographicSponge, W: Witness<F>>
where
    <W as Witness<F>>::MerkleConfig: ark_crypto_primitives::merkle_tree::Config,
{
    config: DirectConfig<W::MerkleConfig, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<W::MerkleConfig>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField, S: CryptographicSponge, W: Witness<F>> Prover<F> for DirectProver<F, S, W>
where
    S::Config: Clone,
    W: Clone,
    W::ChallengeAnswers: Clone,
    W::MerkleConfig: MerkleConfig,
    <W::MerkleConfig as MerkleConfig>::InnerDigest: Absorb,
{
    type Witness = W;
    type Config = DirectConfig<W::MerkleConfig, S>;
    type Proof = SingleProof<F, S, W>;

    fn new(config: DirectConfig<W::MerkleConfig, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<W::MerkleConfig>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn prove(&self, witness: &W) -> Self::Proof {
        <Self::Proof as Proof<F, S, W>>::new(
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
