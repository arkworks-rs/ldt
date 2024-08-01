use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MultiPath, Path},
    sponge::CryptographicSponge,
};
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

#[cfg(not(feature = "std"))]
use ark_std::vec::Vec;

use super::config::FRIConfig;

pub struct FRIProofRound<F: Field, M: MerkleConfig, S: CryptographicSponge> {
    pub queries_to_prev: (Vec<Vec<F>>, Vec<Path<M>>),
    pub challenge_answers: MultiPath<M>,
    pub challenge_values: Vec<Vec<F>>,
    pub config: FRIConfig<M, S>,
    pub last_round_commitment_digest: M::InnerDigest,
}

pub struct FRIProof<F: Field, M: MerkleConfig, S: CryptographicSponge> {
    pub config: FRIConfig<M, S>,
    pub commitment_digests: Vec<M::InnerDigest>,
    pub initial_commitment_digest: M::InnerDigest,
    pub coeff: DensePolynomial<F>,
    pub round_proofs: Vec<FRIProofRound<F, M, S>>,
    // pub round_proofs_2: Vec<FRIProofRound<F, M>>,
    pub proof_of_work_nonce: Option<usize>,
}

impl<F: Field, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> FRIProofRound<F, M, S> {
    pub fn verify_challenge_answers(&self) -> bool {
        self.challenge_answers
            .verify(
                &self.config.merkle_leaf_hash_param,
                &self.config.merkle_two_to_one_param,
                &self.last_round_commitment_digest,
                self.challenge_values.clone(),
            )
            .unwrap()
    }
}

// pub struct FRIProof<F: Field, M: MerkleConfig, S: CryptographicSponge> {
//     pub coeff: DensePolynomial<F>,
//     pub commitment_digest: M::InnerDigest,
//     pub config: FRIConfig<M, S>,
//     pub proof_of_work_nonce: Option<usize>,
// }
