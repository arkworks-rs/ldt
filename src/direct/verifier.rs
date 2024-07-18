use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    direct::{config::DirectConfig, proof::DirectProof},
    ldt::Verifier,
    statement::single::SingleStatement,
    utils::squeeze_integer,
    witness::Witness,
};

pub struct DirectVerifier<F, M, S, W>
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
impl<F, M, S, W> Verifier<F> for DirectVerifier<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    W: Witness<F, M>,
    W::ChallengeAnswers: Clone,
{
    type Statement = SingleStatement<M>;
    type VerifierConfig = DirectConfig<M, S>;
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
    fn verify(&self, claim: &Self::Statement, proof: &Self::Proof) -> bool {
        // regenerate the challenges
        let mut sponge = S::new(&self.config.sponge_config);
        sponge.absorb(&claim.commitment_digest());
        // squeeze out the challenges as indices
        let mut challenges = Vec::with_capacity(self.config.num_challenges);
        for _ in 0..self.config.num_challenges {
            challenges.push(squeeze_integer(&mut sponge, 32));
        }
        // verifiy the proof against the claim
        proof.verify(claim.commitment_digest(), challenges)
    }
}
