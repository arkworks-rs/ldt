use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, Path},
    sponge::{Absorb, CryptographicSponge},
};

use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    commitment::Commitment,
    direct::{config::DirectConfig, proof::DirectProof},
    ldt::Prover,
    utils::squeeze_integer,
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
    type Config = DirectConfig<M, S>;
    type Commitment = Commitment<F, M>;
    type Proof = DirectProof<F, M>;

    fn new(config: DirectConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn prove(&self, commitment: &Self::Commitment) -> Self::Proof {
        // absorb committment
        let p_commitment_root: M::InnerDigest = commitment.p_commitment.root();
        let mut sponge = S::new(&self.config.sponge_config);
        sponge.absorb(&p_commitment_root);
        // squeeze out queries
        let mut queries: Vec<usize> = Vec::with_capacity(self.config.num_queries);
        for _ in 0..self.config.num_queries {
            queries.push(squeeze_integer(&mut sponge, 32));
        }
        // get the openings
        let mut inclusion_proofs: Vec<Path<M>> = Vec::with_capacity(self.config.num_queries);
        for query in queries {
            inclusion_proofs.push(commitment.p_commitment.generate_proof(query).unwrap());
        }
        Self::Proof {
            p_commitment_root,
            p_evaluations: commitment.p_evaluations.clone(),
            inclusion_proofs,
        }
    }
}
