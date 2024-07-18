use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{batch_inversion, FftField, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial,
    Radix2EvaluationDomain,
};
use itertools::izip;

use crate::{
    domain::Domain,
    poly_utils,
    utils::{dedup, proof_of_work_verify, squeeze_integer},
};

use super::{config::STIRConfig, proof::STIRProof};

pub struct STIRVerifierState<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    pub comb_randomness: F,
    pub config: STIRConfig<M, S>,
    pub domain_gen: F,
    pub domain_offset: F,
    pub domain_size: usize,
    pub folding_randomness: F,
    pub interpolating_coeff: DensePolynomial<F>,
    pub proof: STIRProof<F, M, S>,
    pub quotient_set: Vec<F>,
    pub root_of_unity: F,
    pub round_num: usize,
    pub sponge: S,
}

impl<F, M, S> STIRVerifierState<F, M, S>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
{
    pub fn new(
        config: STIRConfig<M, S>,
        commitment_digest: M::InnerDigest,
        proof: STIRProof<F, M, S>,
    ) -> Self {
        let mut sponge = S::new(&config.sponge_config);
        sponge.absorb(&commitment_digest);
        let folding_randomness = sponge.squeeze_field_elements(1)[0];

        let domain = Domain::<F>::new(config.starting_degree, config.starting_rate).unwrap();

        let domain_gen = domain.element(1);
        let domain_size = domain.size();
        Self {
            comb_randomness: F::zero(),
            config,
            domain_gen,
            domain_offset: F::one(),
            domain_size,
            folding_randomness,
            proof,
            interpolating_coeff: DensePolynomial::from_coefficients_vec(vec![]),
            quotient_set: vec![],
            root_of_unity: domain_gen,
            round_num: 0,
            sponge,
        }
    }
    fn answer_evaluations(
        &self,
        coset_offsets: Vec<F>,
        coset_offsets_inv: Vec<F>,
        generator: F,
        generator_inv: F,
        interpolating_coeff: DensePolynomial<F>,
        size: F,
        size_inv: F,
    ) -> Vec<Vec<F>> {
        coset_offsets
            .iter()
            .zip(&coset_offsets_inv)
            .map(|(coset_offset, coset_offset_inv)| match self.round_num {
                0 => vec![F::ONE; self.config.folding_factor],
                _ => {
                    let domain = Radix2EvaluationDomain {
                        size: self.config.folding_factor as u64,
                        log_size_of_group: self.config.folding_factor.ilog2(),
                        size_as_field_element: size,
                        size_inv,
                        group_gen: generator,
                        group_gen_inv: generator_inv,
                        offset: *coset_offset,
                        offset_inv: *coset_offset_inv,
                        offset_pow_size: coset_offset.pow([self.config.folding_factor as u64]),
                    };
                    interpolating_coeff
                        .clone()
                        .evaluate_over_domain(domain)
                        .evals
                }
            })
            .collect()
    }
    fn common_factors(&self, query_sets: Vec<Vec<F>>) -> Vec<Vec<F>> {
        let common_factor_scale = self.comb_randomness;
        query_sets
            .into_iter()
            .map(|query_set| {
                query_set
                    .into_iter()
                    .map(|entry| F::ONE - common_factor_scale * entry)
                    .collect()
            })
            .collect()
    }
    fn coset_offsets(&self, randomness_indices: Vec<usize>) -> Vec<F> {
        randomness_indices
            .iter()
            .map(|stir_randomness_index| {
                self.domain_offset * self.domain_gen.pow([*stir_randomness_index as u64])
            })
            .collect()
    }
    fn denominators(&self, query_sets: Vec<Vec<F>>, quotient_set: Vec<F>) -> Vec<Vec<F>> {
        query_sets
            .iter()
            .map(|query_set| match self.round_num {
                0 => vec![F::ONE; query_set.len()],
                _ => query_set
                    .iter()
                    .map(|eval_point| quotient_set.iter().map(|x| *eval_point - x).product::<F>())
                    .collect(),
            })
            .collect()
    }
    fn folded_answers(
        &self,
        answer_evaluations: &Vec<Vec<F>>,
        common_factors_inv: &Vec<Vec<F>>,
        coset_offsets: &Vec<F>,
        coset_offsets_inv: &Vec<F>,
        denominators_inv: &Vec<Vec<F>>,
        domain_gen: F,
        domain_offset: F,
        generator: F,
        generator_inv: F,
        oracle_answers: Vec<Vec<F>>,
        query_sets: &Vec<Vec<F>>,
        randomness_indices: &Vec<usize>,
        size_inv: F,
    ) -> Vec<(F, F)> {
        let scaled_offset = domain_offset.pow([self.config.folding_factor as u64]);
        let lil_map = izip!(
            0..,
            randomness_indices,
            coset_offsets,
            coset_offsets_inv,
            query_sets,
            common_factors_inv,
            denominators_inv,
            answer_evaluations
        );
        lil_map
            .map(
                |(
                    index,
                    randomness_index,
                    coset_offset,
                    coset_offset_inv,
                    query_set,
                    common_factors_inv,
                    denominators_inv,
                    evaluation_of_ans,
                )| {
                    // This is the point that we are querying at
                    let stir_randomness = scaled_offset
                        * domain_gen.pow([(self.config.folding_factor * randomness_index) as u64]);
                    let f_answers: Vec<_> = query_set
                        .into_iter()
                        .enumerate()
                        .map(|(j, x)| {
                            self.query(
                                *x,
                                oracle_answers[index][j],
                                common_factors_inv[j],
                                denominators_inv[j],
                                evaluation_of_ans[j],
                            )
                        })
                        .collect();
                    // This is the folding
                    let folded_answer = poly_utils::interpolation::fft_interpolate(
                        generator,
                        *coset_offset,
                        generator_inv,
                        *coset_offset_inv,
                        size_inv,
                        &f_answers,
                    )
                    .evaluate(&self.folding_randomness);

                    // Return the folded answer
                    (stir_randomness, folded_answer)
                },
            )
            .collect()
    }
    pub fn folded_evaluations(
        &self,
        randomness_indices: Vec<usize>,
        oracle_answers: Vec<Vec<F>>,
    ) -> Vec<(F, F)> {
        // Step 1: Generator
        let generator: F = self.generator();

        // Step 2: Coset offsets
        let coset_offsets: Vec<F> = self.coset_offsets(randomness_indices.clone());

        // Step 3: Query sets
        let query_sets: Vec<Vec<F>> = self.query_sets(coset_offsets.clone(), generator);

        // Step 4: Common Factors
        let common_factors: Vec<Vec<F>> = self.common_factors(query_sets.clone());

        // Step 5: Denominators
        let denominators = self.denominators(query_sets.clone(), self.quotient_set.clone());

        // Step 6:Invert
        let (
            common_factors_inv,
            coset_offsets_inv,
            denominators_inv,
            generator_inv,
            size,
            size_inv,
        ) = self.invert(
            common_factors.clone(),
            coset_offsets.clone(),
            denominators,
            generator,
        );

        // Step 7: Answer evaluations
        let answer_evaluations = self.answer_evaluations(
            coset_offsets.clone(),
            coset_offsets_inv.clone(),
            generator,
            generator_inv,
            self.interpolating_coeff.clone(),
            size,
            size_inv,
        );

        // Step 8: Folded answer
        self.folded_answers(
            &answer_evaluations,
            &common_factors_inv,
            &coset_offsets,
            &coset_offsets_inv,
            &denominators_inv,
            self.domain_gen,
            self.domain_offset,
            generator,
            generator_inv,
            oracle_answers,
            &query_sets,
            &randomness_indices,
            size_inv,
        )
    }
    fn generator(&self) -> F {
        let scaling_factor = self.domain_size / self.config.folding_factor;
        self.domain_gen.pow([scaling_factor as u64])
    }
    fn invert(
        &self,
        common_factors: Vec<Vec<F>>,
        coset_offsets: Vec<F>,
        denominators: Vec<Vec<F>>,
        generator: F,
    ) -> (Vec<Vec<F>>, Vec<F>, Vec<Vec<F>>, F, F, F) {
        let mut to_invert: Vec<F> = common_factors
            .iter()
            .flatten()
            .chain(denominators.iter().flatten())
            .chain(coset_offsets.iter())
            .cloned()
            .collect();
        to_invert.push(generator);
        let size = F::from(self.config.folding_factor as u64);
        to_invert.push(size);
        batch_inversion(&mut to_invert);
        let size_inv = to_invert.pop().unwrap();
        let generator_inv = to_invert.pop().unwrap();
        let coset_offsets_inv = to_invert.split_off(to_invert.len() - coset_offsets.len());
        let common_factors_len = common_factors.len();
        let chunked: Vec<Vec<F>> = to_invert
            .chunks(self.config.folding_factor)
            .map(|x| x.to_vec())
            .collect();
        let common_factors_inv = chunked[..common_factors_len].to_vec();
        let denominators_inv = chunked[common_factors_len..].to_vec();
        (
            common_factors_inv,
            coset_offsets_inv,
            denominators_inv,
            generator_inv,
            size,
            size_inv,
        )
    }
    pub fn randomness(
        &mut self,
        commitment_digest: M::InnerDigest,
        out_of_domain_evaluations: Vec<F>,
    ) -> (Vec<F>, F, F, Vec<usize>) {
        self.sponge_absorb(&commitment_digest);
        let out_of_domain = self.sponge_squeeze_multiple(self.config.num_out_of_domain_samples);
        self.sponge_absorb(&out_of_domain_evaluations);
        let comb = self.sponge_squeeze();
        let folding = self.sponge_squeeze();
        let scaling_factor = self.domain_size / self.config.folding_factor;
        let num_repetitions = self.config.num_repetitions[self.round_num];
        let indices =
            dedup((0..num_repetitions).map(|_| squeeze_integer(&mut self.sponge, scaling_factor)));
        (out_of_domain, comb, folding, indices)
    }
    pub fn sponge_absorb(&mut self, element: impl Absorb) {
        self.sponge.absorb(&element);
    }
    pub fn sponge_squeeze(&mut self) -> F {
        self.sponge.squeeze_field_elements(1)[0]
    }
    pub fn sponge_squeeze_multiple(&mut self, num_elements: usize) -> Vec<F> {
        self.sponge.squeeze_field_elements(num_elements)
    }
    // TODO: Nuke this
    fn query(
        &self,
        evaluation_point: F,
        value_of_prev_oracle: F,
        common_factors_inverse: F,
        denom_hint: F,
        ans_eval: F,
    ) -> F {
        match &self.round_num {
            0 => value_of_prev_oracle, // In case this is the initial function, we just return the value of the previous oracle
            _ => {
                let num_terms = self.quotient_set.len();
                let quotient_evaluation = poly_utils::quotient::quotient_with_hint(
                    value_of_prev_oracle,
                    evaluation_point,
                    &self.quotient_set,
                    denom_hint,
                    ans_eval,
                );

                let common_factor = evaluation_point * self.comb_randomness;

                let scale_factor = if common_factor != F::ONE {
                    (F::ONE - common_factor.pow([(num_terms + 1) as u64])) * common_factors_inverse
                } else {
                    F::from((num_terms + 1) as u64)
                };

                quotient_evaluation * scale_factor
            }
        }
    }
    fn query_sets(&self, coset_offsets: Vec<F>, generator: F) -> Vec<Vec<F>> {
        let scales: Vec<F> = self.scales(generator);
        coset_offsets
            .iter()
            .map(|coset_offset| {
                (0..self.config.folding_factor)
                    .map(|j| *coset_offset * scales[j])
                    .collect::<Vec<_>>()
            })
            .collect()
    }
    pub fn quotient_answers(
        &self,
        challenge_values: &Vec<Vec<F>>,
        out_of_domain_randomness: &Vec<F>,
        out_of_domain_evaluations: &Vec<F>,
        randomness_indices: &Vec<usize>,
    ) -> Vec<(F, F)> {
        // Step 1: for random indices compute folding of previous oracle TODO: check indices?
        let folded_answers: Vec<(F, F)> =
            self.folded_evaluations(randomness_indices.clone(), challenge_values.clone());

        // Step 2:
        out_of_domain_randomness
            .into_iter()
            .zip(out_of_domain_evaluations)
            .map(|(alpha, beta)| (*alpha, *beta))
            .chain(folded_answers)
            .collect()
    }
    // pub fn randomness_indices(&mut self) -> Vec<usize> {
    //     let final_repetitions = self.config.num_repetitions[self.config.num_rounds];
    //     let scaling_factor = self.domain_size / self.config.folding_factor;
    //     dedup((0..final_repetitions).map(|_| squeeze_integer(&mut self.sponge, scaling_factor)))
    // }
    fn scales(&self, generator: F) -> Vec<F> {
        let scale = generator;
        let mut temp = F::ONE;
        let mut scales = vec![];
        for _ in 0..self.config.folding_factor {
            scales.push(temp);
            temp *= scale;
        }
        scales
    }
    pub fn verify_proof_of_work(&mut self, proof: &STIRProof<F, M, S>) -> bool {
        proof_of_work_verify(
            &mut self.sponge,
            self.config.num_proof_of_work_bits[self.config.num_rounds],
            proof
                .rounds
                .get(self.round_num)
                .unwrap()
                .proof_of_work_nonce,
        )
    }
}
