use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, MerkleTree, Path, TwoToOneParam},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial};

use crate::{
    domain::Domain,
    poly_utils,
    utils::{dedup, proof_of_work, squeeze_integer, stack_evaluations},
};

use super::proof::STIRRoundProof;

pub struct STIRRound<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    answer_coeff: DensePolynomial<F>,
    domain: Domain<F>,
    challenge_answers: Vec<Path<M>>,
    challenge_values: Vec<Vec<F>>,
    challenges: Vec<usize>,
    commitment: MerkleTree<M>,
    committed_values: Vec<Vec<F>>,
    folding_factor: usize,
    folding_randomness: F,
    last_round_commitment: MerkleTree<M>,
    last_round_committed_values: Vec<Vec<F>>,
    last_round_domain_size: usize,
    merkle_leaf_hash_param: LeafParam<M>,
    merkle_two_to_one_param: TwoToOneParam<M>,
    num_out_of_domain_samples: usize,
    num_proof_of_work_bits: Vec<usize>,
    num_repetitions: Vec<usize>,
    num_rounds: usize,
    out_of_domain_samples: Vec<F>,
    out_of_domain_evaluations: Vec<F>,
    proof_of_work_nonce: Option<usize>,
    proximity_generator_randomness: F,
    quotient_answers: Vec<F>,
    quotient_set: Vec<F>,
    round_num: usize,
    shake_coeff: DensePolynomial<F>,
    sponge: S,
    witness_coeff: DensePolynomial<F>,
}

impl<F, M, S> STIRRound<F, M, S>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
{
    pub fn new(
        domain: Domain<F>,
        commitment: MerkleTree<M>,
        committed_values: Vec<Vec<F>>,
        folding_factor: usize,
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        num_out_of_domain_samples: usize,
        num_proof_of_work_bits: Vec<usize>,
        num_repetitions: Vec<usize>,
        num_rounds: usize,
        sponge_config: S::Config,
        witness_coeff: DensePolynomial<F>,
    ) -> Self {
        let mut sponge = S::new(&sponge_config);
        sponge.absorb(&commitment.root());
        Self {
            answer_coeff: DensePolynomial::from_coefficients_vec(vec![]),
            domain: domain.clone(),
            challenge_answers: vec![],
            challenge_values: vec![],
            challenges: vec![],
            commitment,
            committed_values: committed_values.clone(),
            folding_factor,
            folding_randomness: sponge.squeeze_field_elements(1)[0],
            last_round_domain_size: domain.size(),
            last_round_commitment: MerkleTree::<M>::new(
                &merkle_leaf_hash_param,
                &merkle_two_to_one_param,
                &committed_values,
            )
            .unwrap(),
            last_round_committed_values: committed_values,
            merkle_leaf_hash_param,
            merkle_two_to_one_param,
            num_out_of_domain_samples,
            num_proof_of_work_bits,
            num_repetitions,
            num_rounds,
            out_of_domain_samples: vec![],
            out_of_domain_evaluations: vec![],
            proof_of_work_nonce: None,
            proximity_generator_randomness: F::one(),
            quotient_answers: vec![],
            quotient_set: vec![],
            round_num: 0,
            shake_coeff: DensePolynomial::from_coefficients_vec(vec![]),
            sponge,
            witness_coeff,
        }
    }
    fn fold(&mut self) {
        let folded_coeff = poly_utils::folding::poly_fold(
            &self.witness_coeff.clone(),
            self.folding_factor,
            self.folding_randomness,
        );
        if !self.is_final_round() {
            self.last_round_domain_size = self.domain.size();
            let scaled_domain = self.domain.clone().scale_offset(2);
            let evals = folded_coeff
                .evaluate_over_domain_by_ref(scaled_domain.backing_domain)
                .evals;
            let folded_committed_values = stack_evaluations(evals, self.folding_factor);
            self.witness_coeff = folded_coeff;
            self.domain = scaled_domain;
            self.last_round_committed_values = self.committed_values.clone();
            self.committed_values = folded_committed_values;
        }
    }
    fn is_final_round(&self) -> bool {
        self.round_num == self.num_rounds
    }
    pub fn proof(&self) -> STIRRoundProof<F, M> {
        STIRRoundProof {
            commitment_digest: self.commitment.root(),
            out_of_domain_evaluations: self.out_of_domain_evaluations.clone(),
            challenge_values: self.challenge_values.clone(),
            challenge_answers: self.challenge_answers.clone(),
            coeff: self.answer_coeff.clone(),
            is_final_round: self.is_final_round(),
            shake_coeff: self.shake_coeff.clone(),
            proof_of_work_nonce: self.proof_of_work_nonce,
        }
    }
    fn sponge_absorb(&mut self, element: impl Absorb) {
        self.sponge.absorb(&element);
    }
    fn sponge_squeeze(&mut self) -> F {
        self.sponge.squeeze_field_elements(1)[0]
    }
    fn sponge_squeeze_multiple(&mut self, num_elements: usize) -> Vec<F> {
        self.sponge.squeeze_field_elements(num_elements)
    }
    fn update_challenges(&mut self) {
        let (domain_size, committed_values, commitment) = match self.is_final_round() {
            true => (self.domain.size(), &self.committed_values, &self.commitment),
            false => (
                self.last_round_domain_size,
                &self.last_round_committed_values,
                &self.last_round_commitment,
            ),
        };
        self.challenges = dedup(
            (0..self.num_repetitions[self.round_num])
                .map(|_| squeeze_integer(&mut self.sponge, domain_size / self.folding_factor)),
        );
        self.challenge_values = self
            .challenges
            .iter()
            .map(|index| committed_values[*index].clone())
            .collect();
        self.challenge_answers = Vec::with_capacity(self.challenge_values.len());
        for challenge in &self.challenges {
            self.challenge_answers
                .push(commitment.generate_proof(*challenge).unwrap());
        }
    }
    fn update_coeffs(&mut self) {
        // zip set and answers into Vec<(F, F)>
        let zipped: Vec<(F, F)> = self
            .quotient_set
            .clone()
            .into_iter()
            .zip(self.quotient_answers.clone().into_iter())
            .collect();
        // answer_coeff
        self.answer_coeff = poly_utils::interpolation::naive_interpolation(&zipped);
        // shake_coeff
        let mut shake_coeff = DensePolynomial::from_coefficients_vec(vec![]);
        for (x, y) in &zipped {
            let num_coeff = &self.answer_coeff - &DensePolynomial::from_coefficients_vec(vec![*y]);
            let den_coeff = DensePolynomial::from_coefficients_vec(vec![-*x, F::ONE]);
            shake_coeff = shake_coeff + (&num_coeff / &den_coeff);
        }
        self.shake_coeff = shake_coeff;
        // quotient_coeff
        let quotient_coeff =
            poly_utils::quotient::poly_quotient(&self.witness_coeff, &self.quotient_set);
        // scaling_coeff: 1 + r * x + r^2 * x^2 + ... + r^n * x^n
        let scaling_coeff = DensePolynomial::from_coefficients_vec(
            (0..=self.quotient_set.len())
                .map(|i| self.proximity_generator_randomness.pow([i as u64]))
                .collect(),
        );
        // witness_coeff
        self.witness_coeff = &quotient_coeff * &scaling_coeff;
    }
    fn update_commitment(&mut self) {
        self.last_round_commitment = self.commitment.clone();
        self.commitment = MerkleTree::<M>::new(
            &self.merkle_leaf_hash_param,
            &self.merkle_two_to_one_param,
            &self.committed_values,
        )
        .unwrap();
        // put it in the sponge
        self.sponge_absorb(&self.commitment.root());
    }
    fn update_folding_randomness(&mut self) {
        self.folding_randomness = self.sponge_squeeze();
    }
    fn update_out_of_domain_samples(&mut self) {
        self.out_of_domain_samples = self.sponge_squeeze_multiple(self.num_out_of_domain_samples);
        self.out_of_domain_evaluations = self
            .out_of_domain_samples
            .iter()
            .map(|point| self.witness_coeff.evaluate(point))
            .collect();
        self.sponge_absorb(&self.out_of_domain_evaluations.clone());
    }
    fn update_proof_of_work(&mut self) {
        self.proof_of_work_nonce = proof_of_work(
            &mut self.sponge,
            self.num_proof_of_work_bits[self.round_num],
        );
    }
    fn update_proximity_generator_randomness(&mut self) {
        self.proximity_generator_randomness = self.sponge_squeeze();
    }
    fn update_quotient_answers(&mut self) {
        let stir_randomness: Vec<F> = self
            .challenges
            .iter()
            .map(|index| self.domain.scale(self.folding_factor).element(*index))
            .collect();
        self.quotient_set = self
            .out_of_domain_samples
            .clone()
            .into_iter()
            .chain(stir_randomness.iter().cloned())
            .collect();
        self.quotient_answers = self
            .quotient_set
            .iter()
            .map(|x| self.witness_coeff.evaluate(x))
            .collect();
    }
    fn update_round_num(&mut self) {
        self.round_num = self.round_num + 1;
    }
}

impl<F, M, S> Iterator for STIRRound<F, M, S>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
{
    type Item = Self;

    fn next(&mut self) -> Option<Self::Item> {
        if self.round_num < self.num_rounds {
            // Step 1: Perform fold/scale operation
            self.fold();

            if !self.is_final_round() {
                // Step 2: Generate commitment on the folded stuff
                self.update_commitment();

                // Step 3: Out of domain samples
                self.update_out_of_domain_samples();

                // Step 4: Squeeze some randomness
                self.update_proximity_generator_randomness();
                self.update_folding_randomness();
            }

            // Step 5: Generate challenges and answers
            self.update_challenges();

            // Step 6: Proof of work
            self.update_proof_of_work();

            if !self.is_final_round() {
                // Step 7: Squeeze more randomness (used by only verifier)
                let _shake_randomness: F = self.sponge_squeeze();

                // Step 6: Generate quotient set and answers
                self.update_quotient_answers();

                // Step 7: Compute coeffs
                self.update_coeffs();

                // Step 8: Increment
                self.update_round_num();
            }
            Some(self.clone())
        } else {
            None
        }
    }
}

impl<F, M, S> Clone for STIRRound<F, M, S>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
{
    fn clone(&self) -> Self {
        STIRRound {
            answer_coeff: self.answer_coeff.clone(),
            domain: self.domain.clone(),
            challenge_answers: self.challenge_answers.clone(),
            challenge_values: self.challenge_values.clone(),
            challenges: self.challenges.clone(),
            commitment: self.commitment.clone(),
            committed_values: self.committed_values.clone(),
            folding_factor: self.folding_factor,
            folding_randomness: self.folding_randomness,
            last_round_commitment: self.last_round_commitment.clone(),
            last_round_committed_values: self.last_round_committed_values.clone(),
            last_round_domain_size: self.last_round_domain_size,
            merkle_leaf_hash_param: self.merkle_leaf_hash_param.clone(),
            merkle_two_to_one_param: self.merkle_two_to_one_param.clone(),
            num_out_of_domain_samples: self.num_out_of_domain_samples,
            num_proof_of_work_bits: self.num_proof_of_work_bits.clone(),
            num_repetitions: self.num_repetitions.clone(),
            num_rounds: self.num_rounds,
            out_of_domain_samples: self.out_of_domain_samples.clone(),
            out_of_domain_evaluations: self.out_of_domain_evaluations.clone(),
            proof_of_work_nonce: self.proof_of_work_nonce,
            proximity_generator_randomness: self.proximity_generator_randomness,
            quotient_answers: self.quotient_answers.clone(),
            quotient_set: self.quotient_set.clone(),
            round_num: self.round_num,
            shake_coeff: self.shake_coeff.clone(),
            sponge: self.sponge.clone(),
            witness_coeff: self.witness_coeff.clone(),
        }
    }
}
