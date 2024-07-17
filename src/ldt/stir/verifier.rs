use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{batch_inversion, FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, Polynomial, Radix2EvaluationDomain};
use ark_std::marker::PhantomData;
use itertools::izip;

use crate::{
    ldt::{
        stir::{
            config::STIRConfig,
            proof::{STIRProof, STIRProofRound},
            verifier_state::STIRVerifierState,
        },
        Verifier,
    },
    poly_utils,
    statement::single::SingleStatement,
    utils::{dedup, proof_of_work_verify, squeeze_integer},
    witness::Witness,
};

pub struct STIRVerifier<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
    W: Witness<F, M, MerkleConfig = M>,
{
    config: STIRConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
    _witness: PhantomData<W>,
}
impl<F, M, S, W> Verifier<F> for STIRVerifier<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M> + Clone,
    W::ChallengeAnswers: Clone,
{
    type Statement = SingleStatement<M>;
    type VerifierConfig = STIRConfig<M, S>;
    type Proof = STIRProof<F, M, S>;

    fn new(config: STIRConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
            _witness: PhantomData::<W>,
        }
    }
    fn verify(&self, claim: &Self::Statement, proof: &Self::Proof) -> bool {
        // if proof.final_round_proof.coeff.degree() + 1 > self.config.stopping_degree { // TODO: fix this
        //     return false;
        // }

        // Step 1: Verify merkle paths for all rounds
        if !proof.verify_challenge_answers() {
            return false;
        }

        // Step 2: Recompute
        let mut state = STIRVerifierState::new(self.config.clone(), claim.commitment_digest());

        for round_proof in &proof.rounds {
            if !round_proof.is_final_round {
                let round_result = self.round(round_proof, state);
                if round_result.is_none() {
                    return false;
                }
                state = round_result.unwrap();
            }
        }

        // Now, we sample the last points that we want to check consisntency at
        let final_repetitions = self.config.num_repetitions[self.config.num_rounds];
        let scaling_factor = state.domain_size / self.config.folding_factor;
        let final_randomness_indexes = dedup(
            (0..final_repetitions).map(|_| squeeze_integer(&mut state.sponge, scaling_factor)),
        );

        if !proof_of_work_verify(
            &mut state.sponge,
            self.config.num_proof_of_work_bits[self.config.num_rounds],
            proof.rounds.last().unwrap().proof_of_work_nonce,
        ) {
            return false;
        }

        // First, we want to query back the last oracle at this point, which is, again, just a
        // lookup
        let oracle_answers = proof.rounds.last().unwrap().challenge_values.clone();

        let folded_answers =
            self.folded_evaluations(&state, final_randomness_indexes, oracle_answers);

        folded_answers
            .into_iter()
            .all(|(point, value)| proof.rounds.last().unwrap().coeff.evaluate(&point) == value)
    }
}

impl<F, M, S, W> STIRVerifier<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M> + Clone,
    W::ChallengeAnswers: Clone,
{
    fn generator(domain_gen: F, domain_size: usize, folding_factor: usize) -> F {
        let scaling_factor = domain_size / folding_factor;
        domain_gen.pow([scaling_factor as u64])
    }
    fn coset_offsets(domain_gen: F, domain_offset: F, randomness_indices: Vec<usize>) -> Vec<F> {
        randomness_indices
            .iter()
            .map(|stir_randomness_index| {
                domain_offset * domain_gen.pow([*stir_randomness_index as u64])
            })
            .collect()
    }
    fn scales(folding_factor: usize, generator: F) -> Vec<F> {
        let scale = generator;
        let mut temp = F::ONE;
        let mut scales = vec![];
        for _ in 0..folding_factor {
            scales.push(temp);
            temp *= scale;
        }
        scales
    }
    fn query_sets(coset_offsets: Vec<F>, folding_factor: usize, generator: F) -> Vec<Vec<F>> {
        let scales: Vec<F> = Self::scales(folding_factor, generator);
        coset_offsets
            .iter()
            .map(|coset_offset| {
                (0..folding_factor)
                    .map(|j| *coset_offset * scales[j])
                    .collect::<Vec<_>>()
            })
            .collect()
    }
    fn answer_evaluations(
        coset_offsets: Vec<F>,
        coset_offsets_inv: Vec<F>,
        folding_factor: usize,
        generator: F,
        generator_inv: F,
        interpolating_coeff: DensePolynomial<F>,
        round_num: usize,
        size: F,
        size_inv: F,
    ) -> Vec<Vec<F>> {
        coset_offsets
            .iter()
            .zip(&coset_offsets_inv)
            .map(|(coset_offset, coset_offset_inv)| match round_num {
                0 => vec![F::ONE; folding_factor],
                _ => {
                    let domain = Radix2EvaluationDomain {
                        size: folding_factor as u64,
                        log_size_of_group: folding_factor.ilog2(),
                        size_as_field_element: size,
                        size_inv,
                        group_gen: generator,
                        group_gen_inv: generator_inv,
                        offset: *coset_offset,
                        offset_inv: *coset_offset_inv,
                        offset_pow_size: coset_offset.pow([folding_factor as u64]),
                    };
                    interpolating_coeff
                        .clone()
                        .evaluate_over_domain(domain)
                        .evals
                }
            })
            .collect()
    }
    fn common_factors(common_factor_scale: F, query_sets: Vec<Vec<F>>) -> Vec<Vec<F>> {
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
    fn denominators(
        query_sets: Vec<Vec<F>>,
        quotient_set: Vec<F>,
        round_num: usize,
    ) -> Vec<Vec<F>> {
        query_sets
            .iter()
            .map(|query_set| match round_num {
                0 => vec![F::ONE; query_set.len()],
                _ => query_set
                    .iter()
                    .map(|eval_point| quotient_set.iter().map(|x| *eval_point - x).product::<F>())
                    .collect(),
            })
            .collect()
    }
    fn folded_answers(
        answer_evaluations: &Vec<Vec<F>>,
        common_factors_inv: &Vec<Vec<F>>,
        coset_offsets: &Vec<F>,
        coset_offsets_inv: &Vec<F>,
        denominators_inv: &Vec<Vec<F>>,
        domain_gen: F,
        domain_offset: F,
        folding_factor: usize,
        folding_randomness: F,
        generator: F,
        generator_inv: F,
        oracle_answers: Vec<Vec<F>>,
        query_sets: &Vec<Vec<F>>,
        randomness_indices: &Vec<usize>,
        size_inv: F,
        state: &STIRVerifierState<F, M, S>,
    ) -> Vec<(F, F)> {
        let scaled_offset = domain_offset.pow([folding_factor as u64]);
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
                        * domain_gen.pow([(folding_factor * randomness_index) as u64]);
                    let f_answers: Vec<_> = query_set
                        .into_iter()
                        .enumerate()
                        .map(|(j, x)| {
                            state.query(
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
                    .evaluate(&folding_randomness);

                    // Return the folded answer
                    (stir_randomness, folded_answer)
                },
            )
            .collect()
    }
    fn invert(
        common_factors: Vec<Vec<F>>,
        coset_offsets: Vec<F>,
        denominators: Vec<Vec<F>>,
        folding_factor: usize,
        generator: F,
        size: F,
    ) -> (Vec<Vec<F>>, Vec<F>, Vec<Vec<F>>, F, F) {
        let mut to_invert: Vec<F> = common_factors
            .iter()
            .flatten()
            .chain(denominators.iter().flatten())
            .chain(coset_offsets.iter())
            .cloned()
            .collect();
        to_invert.push(generator);
        to_invert.push(size);
        batch_inversion(&mut to_invert);
        let size_inv = to_invert.pop().unwrap();
        let generator_inv = to_invert.pop().unwrap();
        let coset_offsets_inv = to_invert.split_off(to_invert.len() - coset_offsets.len());
        let common_factors_len = common_factors.len();
        let chunked: Vec<Vec<F>> = to_invert
            .chunks(folding_factor)
            .map(|x| x.to_vec())
            .collect();
        let common_factors_inv = chunked[..common_factors_len].to_vec();
        let denominators_inv = chunked[common_factors_len..].to_vec();
        (
            common_factors_inv,
            coset_offsets_inv,
            denominators_inv,
            generator_inv,
            size_inv,
        )
    }
    fn folded_evaluations(
        &self,
        state: &STIRVerifierState<F, M, S>,
        randomness_indices: Vec<usize>,
        oracle_answers: Vec<Vec<F>>,
    ) -> Vec<(F, F)> {
        // Step 1: Generator
        let generator = Self::generator(
            state.domain_gen,
            state.domain_size,
            self.config.folding_factor,
        );
        // Step 2: Coset offsets
        let coset_offsets: Vec<F> = Self::coset_offsets(
            state.domain_gen,
            state.domain_offset,
            randomness_indices.clone(),
        );
        // Step 3: Query sets
        let query_sets: Vec<Vec<F>> =
            Self::query_sets(coset_offsets.clone(), self.config.folding_factor, generator);

        // Step 4: Common Factors
        let common_factors = Self::common_factors(state.comb_randomness, query_sets.clone());

        // Step 5: Denominators
        let denominators = Self::denominators(
            query_sets.clone(),
            state.quotient_set.clone(),
            state.round_num,
        );

        // Step 6:Invert
        let size = F::from(self.config.folding_factor as u64);
        let (common_factors_inv, coset_offsets_inv, denominators_inv, generator_inv, size_inv) =
            Self::invert(
                common_factors.clone(),
                coset_offsets.clone(),
                denominators,
                self.config.folding_factor,
                generator,
                size,
            );

        // Step 7: Answer evaluations
        let answer_evaluations = Self::answer_evaluations(
            coset_offsets.clone(),
            coset_offsets_inv.clone(),
            self.config.folding_factor,
            generator,
            generator_inv,
            state.interpolating_polynomial.clone(),
            state.round_num,
            size,
            size_inv,
        );

        // Step 8: Folded answer
        Self::folded_answers(
            &answer_evaluations,
            &common_factors_inv,
            &coset_offsets,
            &coset_offsets_inv,
            &denominators_inv,
            state.domain_gen,
            state.domain_offset,
            self.config.folding_factor,
            state.folding_randomness,
            generator,
            generator_inv,
            oracle_answers,
            &query_sets,
            &randomness_indices,
            size_inv,
            state,
        )
    }
    fn round(
        &self,
        round_proof: &STIRProofRound<F, M, S>,
        mut state: STIRVerifierState<F, M, S>,
    ) -> Option<STIRVerifierState<F, M, S>> {
        // Redo FS
        state.sponge_absorb(&round_proof.commitment_digest);
        let ood_randomness = state.sponge_squeeze_multiple(self.config.num_out_of_domain_samples);
        state.sponge_absorb(&round_proof.out_of_domain_evaluations);
        let comb_randomness = state.sponge_squeeze();
        let new_folding_randomness = state.sponge_squeeze();
        let scaling_factor = state.domain_size / self.config.folding_factor;

        let num_repetitions = self.config.num_repetitions[state.round_num];
        let stir_randomness_indexes =
            dedup((0..num_repetitions).map(|_| squeeze_integer(&mut state.sponge, scaling_factor)));

        // PoW verification
        if !proof_of_work_verify(
            &mut state.sponge,
            self.config.num_proof_of_work_bits[state.round_num],
            round_proof.proof_of_work_nonce,
        ) {
            return None;
        }

        let shake_randomness = state.sponge_squeeze();

        // Now, we are starting to define the next function.
        // First, we need to query the previous oracle (which is either f_0 or g_i)
        // At the indexes B_i for i in stir_randomness_indexes
        // Since we previously verified the Merkle paths, this is easy
        // TODO: We should probably check the indexes
        let oracle_answers = round_proof.challenge_values.clone();

        // Now, for each of the selected random points, we need to compute the folding of the
        // previous oracle
        let folded_answers =
            self.folded_evaluations(&state, stir_randomness_indexes, oracle_answers);

        // The quotient definining the function
        let quotient_answers: Vec<_> = ood_randomness
            .into_iter()
            .zip(&round_proof.out_of_domain_evaluations)
            .map(|(alpha, beta)| (alpha, *beta))
            .chain(folded_answers)
            .collect();
        let interpolating_polynomial = round_proof.coeff.clone();

        let ans_eval = interpolating_polynomial.evaluate(&shake_randomness);
        let shake_eval = round_proof.shake_coeff.evaluate(&shake_randomness);

        let mut denoms: Vec<_> = quotient_answers
            .iter()
            .map(|(x, _)| shake_randomness - x)
            .collect();

        batch_inversion(&mut denoms);
        // TODO: This maybe should be better
        if shake_eval
            != quotient_answers
                .iter()
                .zip(denoms)
                .map(|((_, y), d)| (ans_eval - y) * d)
                .sum()
        {
            return None;
        }

        let quotient_set = quotient_answers
            .into_iter()
            .map(|(x, _)| x)
            .collect::<Vec<_>>();

        Some(STIRVerifierState {
            comb_randomness: comb_randomness.clone(),
            config: self.config.clone(),
            domain_gen: state.domain_gen * state.domain_gen,
            domain_offset: state.domain_offset * state.domain_offset * state.root_of_unity,
            domain_size: state.domain_size / 2,
            folding_randomness: new_folding_randomness,
            interpolating_polynomial: interpolating_polynomial.clone(),
            quotient_set,
            root_of_unity: state.root_of_unity,
            round_num: state.round_num + 1,
            sponge: state.sponge,
        })
    }
}
