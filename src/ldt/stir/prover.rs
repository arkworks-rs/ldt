use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::Polynomial;
use ark_std::marker::PhantomData;

use crate::{
    ldt::{
        stir::{config::STIRConfig, proof::STIRProof, round::STIRRound},
        Prover,
    },
    witness::Witness,
};

pub struct STIRProver<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    W: Witness<F, M, MerkleConfig = M>,
{
    config: STIRConfig<W::MerkleConfig, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<W::MerkleConfig>,
    _sponge_config: PhantomData<S>,
}

impl<F, M, S, W> Prover<F> for STIRProver<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M, Commitment = MerkleTree<M>, CommittedValues = Vec<Vec<F>>>
        + Clone,
    W::ChallengeAnswers: Clone,
{
    type Witness = W;
    type ProverConfig = STIRConfig<M, S>;
    type Proof = STIRProof<F, M>;

    fn new(config: STIRConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn prove(&self, witness: &W) -> Self::Proof {
        assert!(witness.coeff().degree() < self.config.starting_degree);

        let proofs = STIRRound::<F, M, S>::new(
            witness.domain(),
            witness.commitment(),
            witness.committed_values(),
            self.config.clone(),
            witness.coeff(),
        )
        .map(|round| round.proof())
        .collect();

        STIRProof::<F, M> {
            round_proofs: proofs,
        }
    }
}
