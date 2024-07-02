use ark_crypto_primitives::sponge::{Absorb, CryptographicSponge};
use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    commitment::Witness,
    direct::config::DirectConfig,
    ldt::Verifier,
    proof::{Proof, SingleProof},
};

pub struct DirectVerifier<F: FftField, S: CryptographicSponge, W: Witness<F>>
where
    <W as Witness<F>>::MerkleConfig: ark_crypto_primitives::merkle_tree::Config,
{
    config: DirectConfig<W::MerkleConfig, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<W::MerkleConfig>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField, S: CryptographicSponge, W: Witness<F>> Verifier<F> for DirectVerifier<F, S, W>
where
    <W as Witness<F>>::MerkleConfig: ark_crypto_primitives::merkle_tree::Config,
    <<W as Witness<F>>::MerkleConfig as ark_crypto_primitives::merkle_tree::Config>::InnerDigest:
        Absorb,
    <W as Witness<F>>::ChallengeAnswers: Clone,
{
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
    fn verify(&self, proof: &Self::Proof) -> bool {
        proof.verify()
    }
}
