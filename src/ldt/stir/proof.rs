use ark_crypto_primitives::merkle_tree::{Config as MerkleConfig, Path};
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

pub struct STIRRoundProof<F: Field, M: MerkleConfig> {
    pub coeff: DensePolynomial<F>,
    pub challenge_answers: Vec<Path<M>>,
    pub committed_values: Vec<Vec<F>>,
    pub is_final_round: bool,
    pub out_of_domain_evaluations: Vec<F>, // Note: empty when is_final_round = true
    pub commitment_digest: M::InnerDigest,
    pub proof_of_work_nonce: Option<usize>,
    pub shake_coeff: DensePolynomial<F>, // Note: empty when is_final_round = true
}

impl<F: Field, M: MerkleConfig> STIRRoundProof<F, M> {
    pub fn new(
        coeff: DensePolynomial<F>,
        challenge_answers: Vec<Path<M>>,
        committed_values: Vec<Vec<F>>,
        is_final_round: bool,
        out_of_domain_evaluations: Vec<F>,
        commitment_digest: M::InnerDigest,
        proof_of_work_nonce: Option<usize>,
        shake_coeff: DensePolynomial<F>,
    ) -> Self {
        STIRRoundProof {
            coeff,
            challenge_answers,
            committed_values,
            is_final_round,
            out_of_domain_evaluations,
            commitment_digest,
            proof_of_work_nonce,
            shake_coeff,
        }
    }
}

impl<F: Field, M: MerkleConfig> Clone for STIRRoundProof<F, M> {
    fn clone(&self) -> Self {
        STIRRoundProof {
            coeff: self.coeff.clone(),
            challenge_answers: self.challenge_answers.clone(),
            committed_values: self.committed_values.clone(),
            is_final_round: self.is_final_round,
            out_of_domain_evaluations: self.out_of_domain_evaluations.clone(),
            commitment_digest: self.commitment_digest.clone(),
            proof_of_work_nonce: self.proof_of_work_nonce.clone(),
            shake_coeff: self.shake_coeff.clone(),
        }
    }
}

pub struct STIRProof<F: Field, M: MerkleConfig> {
    pub round_proofs: Vec<STIRRoundProof<F, M>>,
}

impl<F: Field, M: MerkleConfig> Clone for STIRProof<F, M> {
    fn clone(&self) -> Self {
        STIRProof {
            round_proofs: self.round_proofs.clone(),
        }
    }
}
