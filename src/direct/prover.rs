use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree, Path},
    sponge::{Absorb, CryptographicSponge},
};

use ark_ff::FftField;
use ark_poly::univariate::DensePolynomial;
use ark_std::marker::PhantomData;

use crate::{direct::config::DirectConfig, domain::Domain, ldt::Prover, utils::squeeze_integer};

pub struct DirectCommitment<F: FftField, M: MerkleConfig> {
    pub p_commitment: MerkleTree<M>,
    pub p_evaluations_over_domain: Vec<Vec<F>>,
}
pub struct DirectProof<F: FftField, M: MerkleConfig> {
    pub commitment: DirectCommitment<F, M>,
    pub query_inclusion_proofs: Vec<Path<M>>,
}

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
    type Commitment = DirectCommitment<F, M>;
    type Proof = DirectProof<F, M>;

    fn new(config: DirectConfig<M, S>) -> Self {
        DirectProver {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn commit(&self, polynomial: DensePolynomial<F>) -> DirectCommitment<F, M> {
        // get evaluations over a domain
        let domain = Domain::<F>::new(self.config.degree, 0).unwrap();
        let p_evaluations_over_domain: Vec<Vec<F>> = polynomial
            .evaluate_over_domain_by_ref(domain.backing_domain)
            .evals
            .iter()
            .map(|f| -> Vec<F> { vec![*f] })
            .collect();
        // generate the committment
        let p_commitment = MerkleTree::<M>::new(
            &self.config.merkle_leaf_hash_param,
            &self.config.merkle_two_to_one_param,
            p_evaluations_over_domain.clone(),
        )
        .unwrap();
        DirectCommitment {
            p_evaluations_over_domain,
            p_commitment,
        }
    }
    fn prove(&self, commitment: DirectCommitment<F, M>) -> DirectProof<F, M> {
        // absorb committment
        let root_hash: M::InnerDigest = commitment.p_commitment.root();
        let mut sponge = S::new(&self.config.sponge_config);
        sponge.absorb(&root_hash);
        // squeeze out queries
        let mut queries: Vec<usize> = Vec::with_capacity(self.config.num_queries);
        for _ in 0..self.config.num_queries {
            queries.push(squeeze_integer(&mut sponge, 32));
        }
        // get the openings
        let mut query_inclusion_proofs: Vec<Path<M>> = Vec::with_capacity(self.config.num_queries);
        for query in queries {
            query_inclusion_proofs.push(commitment.p_commitment.generate_proof(query).unwrap());
        }
        DirectProof {
            commitment,
            query_inclusion_proofs,
        }
    }
}
