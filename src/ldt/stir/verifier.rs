use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{batch_inversion, FftField, PrimeField};
use ark_poly::{Polynomial, Radix2EvaluationDomain};
use ark_std::marker::PhantomData;

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
            self.compute_folded_evaluations(&state, final_randomness_indexes, oracle_answers);

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
    fn compute_folded_evaluations(
        &self,
        verification_state: &STIRVerifierState<F, M, S>,
        stir_randomness_indexes: Vec<usize>,
        oracle_answers: Vec<Vec<F>>,
    ) -> Vec<(F, F)> {
        let scaling_factor = verification_state.domain_size / self.config.folding_factor;
        let generator = verification_state.domain_gen.pow([scaling_factor as u64]);

        // We do a single batch inversion
        let coset_offsets: Vec<_> = stir_randomness_indexes
            .iter()
            .map(|stir_randomness_index| {
                verification_state.domain_offset
                    * verification_state
                        .domain_gen
                        .pow([*stir_randomness_index as u64])
            })
            .collect();

        // We use this to more efficiently compute query_sets
        let scale = generator;
        let mut temp = F::ONE;
        let mut scales = vec![];
        for _ in 0..self.config.folding_factor {
            scales.push(temp);
            temp *= scale;
        }

        let query_sets: Vec<_> = coset_offsets
            .iter()
            .map(|coset_offset| {
                (0..self.config.folding_factor)
                    .map(|j| *coset_offset * scales[j])
                    .collect::<Vec<_>>()
            })
            .collect();

        let common_factor_scale = verification_state.comb_randomness;

        let global_common_factors = query_sets
            .iter()
            .map(|query_set| query_set.iter().map(|x| F::ONE - common_factor_scale * x));

        let global_denominators =
            query_sets
                .iter()
                .map(|query_set| match &verification_state.round_num {
                    0 => vec![F::ONE; query_set.len()],
                    _ => query_set
                        .iter()
                        .map(|eval_point| {
                            verification_state
                                .quotient_set
                                .iter()
                                .map(|x| *eval_point - x)
                                .product::<F>()
                        })
                        .collect::<Vec<_>>(),
                });

        // To invert contains a bunch of stuff offsets, generator, size, and common factors
        let size = F::from(self.config.folding_factor as u64);
        let mut to_invert = vec![];
        let global_common_factors_len = global_common_factors.len();
        for common_factors in global_common_factors {
            to_invert.extend(common_factors);
        }
        for denominators in global_denominators {
            to_invert.extend(denominators);
        }
        to_invert.extend(coset_offsets.iter());
        to_invert.push(generator);
        to_invert.push(size);
        batch_inversion(&mut to_invert);
        let size_inv = to_invert.pop().unwrap();
        let generator_inv = to_invert.pop().unwrap();
        let coset_offsets_inv = to_invert.split_off(to_invert.len() - coset_offsets.len());
        let chunked: Vec<Vec<_>> = to_invert
            .chunks(self.config.folding_factor)
            .map(|x| x.to_vec())
            .collect();

        // TODO: Could be split_off
        let common_factors_inv = chunked[0..global_common_factors_len].to_vec();
        let denominators_inv = chunked[global_common_factors_len..].to_vec();

        let evaluations_of_ans: Vec<_> = coset_offsets
            .iter()
            .zip(&coset_offsets_inv)
            .map(
                |(coset_offset, coset_offset_inv)| match &verification_state.round_num {
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

                        verification_state
                            .interpolating_polynomial
                            .clone()
                            .evaluate_over_domain(domain)
                            .evals
                    }
                },
            )
            .collect();

        let scaled_offset = verification_state
            .domain_offset
            .pow([self.config.folding_factor as u64]);

        stir_randomness_indexes
            .iter()
            .zip(coset_offsets)
            .zip(coset_offsets_inv)
            .zip(query_sets)
            .zip(common_factors_inv)
            .zip(denominators_inv)
            .zip(evaluations_of_ans)
            .enumerate()
            // Just restructure
            .map(
                |(
                    i,
                    (
                        (
                            (
                                (
                                    ((stir_randomness_index, coset_offset), coset_offset_inv),
                                    query_set,
                                ),
                                common_factors_inv,
                            ),
                            denominators_inv,
                        ),
                        evaluation_of_ans,
                    ),
                )| {
                    (
                        i,
                        stir_randomness_index,
                        coset_offset,
                        coset_offset_inv,
                        query_set,
                        common_factors_inv,
                        denominators_inv,
                        evaluation_of_ans,
                    )
                },
            )
            .map(
                |(
                    i,
                    stir_randomness_index,
                    coset_offset,
                    coset_offset_inv,
                    query_set,
                    common_factors_inv,
                    denominators_inv,
                    evaluation_of_ans,
                )| {
                    // This is the point that we are querying at
                    let stir_randomness = scaled_offset
                        * verification_state
                            .domain_gen
                            .pow([(self.config.folding_factor * stir_randomness_index) as u64]);

                    let f_answers: Vec<_> = query_set
                        .into_iter()
                        .enumerate()
                        .map(|(j, x)| {
                            verification_state.query(
                                x,
                                oracle_answers[i][j],
                                common_factors_inv[j],
                                denominators_inv[j],
                                evaluation_of_ans[j],
                            )
                        })
                        .collect();

                    // This is the folding
                    let folded_answer = poly_utils::interpolation::fft_interpolate(
                        generator,
                        coset_offset,
                        generator_inv,
                        coset_offset_inv,
                        size_inv,
                        &f_answers,
                    )
                    .evaluate(&verification_state.folding_randomness);

                    // Return the folded answer
                    (stir_randomness, folded_answer)
                },
            )
            .collect()
    }
    fn round(
        &self,
        round_proof: &STIRProofRound<F, M, S>,
        mut verification_state: STIRVerifierState<F, M, S>,
    ) -> Option<STIRVerifierState<F, M, S>> {
        // Redo FS
        verification_state.sponge_absorb(&round_proof.commitment_digest);
        let ood_randomness =
            verification_state.sponge_squeeze_multiple(self.config.num_out_of_domain_samples);
        verification_state.sponge_absorb(&round_proof.out_of_domain_evaluations);
        let comb_randomness = verification_state.sponge_squeeze();
        let new_folding_randomness = verification_state.sponge_squeeze();
        let scaling_factor = verification_state.domain_size / self.config.folding_factor;

        let num_repetitions = self.config.num_repetitions[verification_state.round_num];
        let stir_randomness_indexes = dedup(
            (0..num_repetitions)
                .map(|_| squeeze_integer(&mut verification_state.sponge, scaling_factor)),
        );

        // PoW verification
        if !proof_of_work_verify(
            &mut verification_state.sponge,
            self.config.num_proof_of_work_bits[verification_state.round_num],
            round_proof.proof_of_work_nonce,
        ) {
            return None;
        }

        let shake_randomness = verification_state.sponge_squeeze();

        // Now, we are starting to define the next function.
        // First, we need to query the previous oracle (which is either f_0 or g_i)
        // At the indexes B_i for i in stir_randomness_indexes
        // Since we previously verified the Merkle paths, this is easy
        // TODO: We should probably check the indexes
        let oracle_answers = round_proof.challenge_values.clone();

        // Now, for each of the selected random points, we need to compute the folding of the
        // previous oracle
        let folded_answers = self.compute_folded_evaluations(
            &verification_state,
            stir_randomness_indexes,
            oracle_answers,
        );

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
            domain_gen: verification_state.domain_gen * verification_state.domain_gen,
            domain_offset: verification_state.domain_offset
                * verification_state.domain_offset
                * verification_state.root_of_unity,
            domain_size: verification_state.domain_size / 2,
            folding_randomness: new_folding_randomness,
            interpolating_polynomial: interpolating_polynomial.clone(),
            quotient_set,
            root_of_unity: verification_state.root_of_unity,
            round_num: verification_state.round_num + 1,
            sponge: verification_state.sponge,
        })
    }
}
