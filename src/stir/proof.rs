use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, Path},
    sponge::CryptographicSponge,
};
use ark_ff::{batch_inversion, Field};
use ark_poly::{univariate::DensePolynomial, Polynomial};

use super::config::STIRConfig;

pub struct STIRProofRound<F: Field, M: MerkleConfig, S: CryptographicSponge> {
    pub coeff: DensePolynomial<F>,
    pub challenge_answers: Vec<Path<M>>,
    pub challenge_values: Vec<Vec<F>>,
    pub commitment_digest: M::InnerDigest,
    pub config: STIRConfig<M, S>,
    pub is_final_round: bool,
    pub last_round_commitment_digest: M::InnerDigest,
    pub out_of_domain_evaluations: Vec<F>,
    pub proof_of_work_nonce: Option<usize>,
    pub shake_coeff: DensePolynomial<F>,
}

impl<F: Field, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> STIRProofRound<F, M, S> {
    pub fn new(
        coeff: DensePolynomial<F>,
        challenge_answers: Vec<Path<M>>,
        challenge_values: Vec<Vec<F>>,
        commitment_digest: M::InnerDigest,
        config: STIRConfig<M, S>,
        is_final_round: bool,
        last_round_commitment_digest: M::InnerDigest,
        out_of_domain_evaluations: Vec<F>,
        proof_of_work_nonce: Option<usize>,
        shake_coeff: DensePolynomial<F>,
    ) -> Self {
        STIRProofRound {
            coeff,
            challenge_answers,
            challenge_values,
            commitment_digest,
            config,
            is_final_round,
            last_round_commitment_digest,
            out_of_domain_evaluations,
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
                    &self.last_round_commitment_digest,
                    challenge_value,
                )
                .unwrap()
            {
                return false;
            }
        }
        true
    }
    pub fn verify_quotient_answers(
        &self,
        quotient_answers: &Vec<(F, F)>,
        shake_randomness: &F,
    ) -> bool {
        let ans_eval = self.coeff.evaluate(&shake_randomness);
        let mut denominators: Vec<F> = quotient_answers
            .iter()
            .map(|(x, _)| *shake_randomness - x)
            .collect();
        batch_inversion(&mut denominators);
        let shake_eval = self.shake_coeff.evaluate(&shake_randomness);
        if shake_eval
            != quotient_answers
                .iter()
                .zip(denominators)
                .map(|((_, y), d)| (ans_eval - y) * d)
                .sum()
        {
            return false;
        }
        true
    }
}

impl<F, M, S> Clone for STIRProofRound<F, M, S>
where
    F: Field,
    M: MerkleConfig,
    S: CryptographicSponge,
    S::Config: Clone,
{
    fn clone(&self) -> Self {
        STIRProofRound {
            coeff: self.coeff.clone(),
            challenge_answers: self.challenge_answers.clone(),
            challenge_values: self.challenge_values.clone(),
            commitment_digest: self.commitment_digest.clone(),
            config: self.config.clone(),
            is_final_round: self.is_final_round,
            last_round_commitment_digest: self.last_round_commitment_digest.clone(),
            out_of_domain_evaluations: self.out_of_domain_evaluations.clone(),
            proof_of_work_nonce: self.proof_of_work_nonce.clone(),
            shake_coeff: self.shake_coeff.clone(),
        }
    }
}

pub struct STIRProof<F: Field, M: MerkleConfig, S: CryptographicSponge> {
    pub rounds: Vec<STIRProofRound<F, M, S>>,
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
            rounds: self.rounds.clone(),
        }
    }
}
