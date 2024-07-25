use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::Polynomial;
use ark_std::marker::PhantomData;

#[cfg(not(feature = "std"))]
use ark_std::vec::Vec;

use crate::{
    ldt::Verifier,
    statement::single::SingleStatement,
    stir::{config::STIRConfig, proof::STIRProof, verifier_state::STIRVerifierState},
    witness::Witness,
};

pub struct STIRVerifier<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
    W: Witness<F, M, MerkleConfig = M>,
{
    config: STIRConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge: PhantomData<S>,
    _witness: PhantomData<W>,
}
impl<F, M, S, W> Verifier<F> for STIRVerifier<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>> + Clone,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M> + Clone,
    W::ChallengeAnswers: Clone,
{
    type Statement = SingleStatement<M>;
    type VerifierConfig = STIRConfig<M, S>;
    type Proof = STIRProof<F, M, S>;

    fn new(config: STIRConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge: PhantomData::<S>,
            _witness: PhantomData::<W>,
        }
    }
    fn verify(&self, claim: &Self::Statement, proof: &Self::Proof) -> bool {
        if proof.rounds.last().unwrap().coeff.degree() + 1 > self.config.stopping_degree {
            return false;
        }

        STIRVerifierState::new(
            self.config.clone(),
            claim.commitment_digest(),
            proof.clone(),
        )
        .into_iter()
        .all(|state| state.is_verified() == true)
    }
}
