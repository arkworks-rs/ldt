use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, Path},
    sponge::CryptographicSponge,
};
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

use super::config::STIRConfig;

pub struct STIRRoundProof<F: Field, M: MerkleConfig, S: CryptographicSponge> {
    pub coeff: DensePolynomial<F>,
    pub challenge_answers: Vec<Path<M>>,
    pub challenge_values: Vec<Vec<F>>,
    pub config: STIRConfig<M, S>,
    pub is_final_round: bool,
    pub out_of_domain_evaluations: Vec<F>, // Note: empty when is_final_round = true
    pub commitment_digest: M::InnerDigest,
    pub proof_of_work_nonce: Option<usize>,
    pub shake_coeff: DensePolynomial<F>, // Note: empty when is_final_round = true
}

impl<F: Field, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> STIRRoundProof<F, M, S> {
    pub fn new(
        coeff: DensePolynomial<F>,
        challenge_answers: Vec<Path<M>>,
        challenge_values: Vec<Vec<F>>,
        config: STIRConfig<M, S>,
        is_final_round: bool,
        out_of_domain_evaluations: Vec<F>,
        commitment_digest: M::InnerDigest,
        proof_of_work_nonce: Option<usize>,
        shake_coeff: DensePolynomial<F>,
    ) -> Self {
        STIRRoundProof {
            coeff,
            challenge_answers,
            challenge_values,
            config,
            is_final_round,
            out_of_domain_evaluations,
            commitment_digest,
            proof_of_work_nonce,
            shake_coeff,
        }
    }
    pub fn verify_challenge_answers(&self) -> bool {
        for (challenge_value, challenge_answer) in self
            .challenge_values
            .iter()
            .zip(self.challenge_answers.iter())
        {
            if !challenge_answer
                .verify(
                    &self.config.merkle_leaf_hash_param,
                    &self.config.merkle_two_to_one_param,
                    &self.commitment_digest,
                    challenge_value,
                )
                .unwrap()
            {
                return false;
            }
        }
        true
    }
}

impl<F, M, S> Clone for STIRRoundProof<F, M, S>
where
    F: Field,
    M: MerkleConfig,
    S: CryptographicSponge,
    S::Config: Clone,
{
    fn clone(&self) -> Self {
        STIRRoundProof {
            coeff: self.coeff.clone(),
            challenge_answers: self.challenge_answers.clone(),
            challenge_values: self.challenge_values.clone(),
            config: self.config.clone(),
            is_final_round: self.is_final_round,
            out_of_domain_evaluations: self.out_of_domain_evaluations.clone(),
            commitment_digest: self.commitment_digest.clone(),
            proof_of_work_nonce: self.proof_of_work_nonce.clone(),
            shake_coeff: self.shake_coeff.clone(),
        }
    }
}

pub struct STIRProof<F: Field, M: MerkleConfig, S: CryptographicSponge> {
    pub round_proofs: Vec<STIRRoundProof<F, M, S>>,
}

impl<F, M, S> Clone for STIRProof<F, M, S>
where
    F: Field,
    M: MerkleConfig,
    S: CryptographicSponge,
    S::Config: Clone,
{
    fn clone(&self) -> Self {
        STIRProof {
            round_proofs: self.round_proofs.clone(),
        }
    }
}
