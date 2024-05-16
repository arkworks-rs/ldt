use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};

use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{direct::config::DirectConfig, utils::squeeze_integer};

use super::prover::DirectProof;

pub struct DirectVerifier<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    config: DirectConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}

impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> DirectVerifier<F, M, S>
where
    M::InnerDigest: Absorb,
{
    pub fn new(config: DirectConfig<M, S>) -> Self {
        DirectVerifier {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }
    pub fn verify(&self, proof: &DirectProof<F, M>) -> bool {
        // absorb the committment to derive the queries
        let commitment: M::InnerDigest = proof.commitment.p_commitment.root();
        let mut sponge: S = S::new(&self.config.sponge_config);
        sponge.absorb(&commitment);

        // for all the queries
        for query in 0_usize..self.config.num_queries {
            // squeeze out the query
            let leaf_index_of_query: usize = squeeze_integer(
                &mut sponge,
                proof.commitment.p_evaluations_over_domain.len(),  // TODO: (z-tech) verify evals.len() is always power of 2
            );

            // verify query was derived properly
            let is_correct_leaf_index =
                proof.query_inclusion_proofs[query].leaf_index == leaf_index_of_query;
            if !is_correct_leaf_index {
                return false;
            }

            // is valid proof of inclusion at index
            let is_valid_proof = proof.query_inclusion_proofs[query]
                .verify(
                    &self.config.merkle_leaf_hash_param,
                    &self.config.merkle_two_to_one_param,
                    &commitment,
                    proof.commitment.p_evaluations_over_domain[leaf_index_of_query].clone(),
                )
                .unwrap();
            if !is_valid_proof {
                return false;
            }
        }
        true
    }
}
