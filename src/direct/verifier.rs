use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::FftField;
use ark_std::marker::PhantomData;

use crate::{
    commitment::Commitment, direct::config::DirectConfig, direct::proof::DirectProof,
    ldt::Verifier, utils::squeeze_integer,
};

pub struct DirectVerifier<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    config: DirectConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> Verifier<F>
    for DirectVerifier<F, M, S>
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
    fn verify(&self, proof: &Self::Proof) -> bool {
        // absorb the committment to derive the queries
        let root_hash: M::InnerDigest = proof.p_commitment_root.clone();
        let mut sponge: S = S::new(&self.config.sponge_config);
        sponge.absorb(&root_hash);

        // for all the queries
        for query in 0_usize..self.config.num_queries {
            // squeeze out the query
            let leaf_index_of_query: usize = squeeze_integer(
                &mut sponge,
                proof.p_evaluations.len(), // TODO: (z-tech) verify evals.len() is always power of 2
            );

            // verify query was derived properly
            let is_correct_leaf_index =
                proof.inclusion_proofs[query].leaf_index == leaf_index_of_query;
            if !is_correct_leaf_index {
                return false;
            }

            // is valid proof of inclusion at index
            let is_valid_proof = proof.inclusion_proofs[query]
                .verify(
                    &self.config.merkle_leaf_hash_param,
                    &self.config.merkle_two_to_one_param,
                    &root_hash,
                    proof.p_evaluations[leaf_index_of_query].clone(),
                )
                .unwrap();
            if !is_valid_proof {
                return false;
            }
        }
        true
    }
}
