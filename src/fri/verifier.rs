use crate::fri::FRIParameters;
use ark_ff::PrimeField;
use ark_sponge::{FieldBasedCryptographicSponge, FieldElementSize};
use std::marker::PhantomData;
use crate::domain::Radix2CosetDomain;
use ark_poly::Polynomial;
use crate::direct::DirectLDT;

/// Implements FRI verifier.
pub struct FRIVerifier<F: PrimeField> {
    _verifier: PhantomData<F>,
}

impl<F: PrimeField> FRIVerifier<F> {
    /// ## Step 1: Interative Phase
    /// Sample alpha in interactive phase.
    pub fn interactive_phase_single_round<S: FieldBasedCryptographicSponge<F>>(
        sponge: &mut S,
    ) -> F {
        sponge.squeeze_native_field_elements(1)[0]
    }

    /// ## Step 2: Sample Queried Coset
    /// Sample the coset to be queried.
    pub fn sample_coset_index<S: FieldBasedCryptographicSponge<F>>(
        sponge: &mut S,
        fri_parameters: &FRIParameters<F>
    ) -> usize{
        let log_num_cosets = fri_parameters.domain.dim();
        // we use the fact that number of cosets is always power of two
        let rand_coset_index = le_bits_array_to_usize(&sponge
            .squeeze_bits(log_num_cosets));
        rand_coset_index
    }

    /// ## Step 2: Query Phase (Prepare Query)
    /// Prepare all the queries in query phase. The returned value `queries[i]` is the coset query
    /// of the `ith` round polynomial (including codeword polynomial but does not include final polynomial).
    /// Final polynomial is not queried. Instead, verifier will get
    /// the whole final polynomial in evaluation form, and do direct LDT.
    pub fn prepare_queries(
        rand_coset_index: usize,
        fri_parameters: &FRIParameters<F>,
    ) -> (Vec<Radix2CosetDomain<F>>) {
        let num_fri_rounds = fri_parameters.localization_parameters.len();

        // this does not include domain.offset! multiply domain.offset in use
        let mut curr_coset_index = rand_coset_index;
        let mut queries = Vec::with_capacity(num_fri_rounds);
        let mut curr_domain_coset_size = fri_parameters.domain.size();
        // sample a coset index
        for i in 0..num_fri_rounds {
            // current coset index = last coset index % (distance between coset at current round)
            // edge case: at first round, this still applies

            let dist_between_coset_elems = curr_domain_coset_size / (1 << fri_parameters.localization_parameters[i]);
            curr_coset_index = curr_coset_index & ((1 << dist_between_coset_elems) - 1);

            let (_, query_coset) = fri_parameters.domain.query_position_to_coset(curr_coset_index,
                                                                     fri_parameters.localization_parameters[i] as usize);

            queries.push(query_coset);

            // get next round coset size
            curr_domain_coset_size = dist_between_coset_elems;
        }

        queries
    }

    /// ## Step 3: Query Phase (Check query)
    /// After preparing the query, verifier get the evaluations of corresponding query. Those evaluations needs
    /// to be checked by merkle tree. Then verifier calls this method to check if polynomial sent in each round
    /// is consistent with each other, and the final polynomial is low-degree.
    ///
    /// `queries[i]` is the coset query of the `ith` round polynomial, including the codeword polynomial.
    /// `queried_evaluations` stores the result of corresponding query.
    pub fn consistency_check(
        fri_parameters: &FRIParameters<F>,
        queried_coset_index: usize,
        queries: &[Radix2CosetDomain<F>],
        queried_evaluations: &[Vec<F>],
        alphas: &[F],
        final_polynomial_domain: &Radix2CosetDomain<F>,
        final_polynomial: &[F]
    ) -> bool {
        let mut curr_coset_index = queried_coset_index;
        // this does not include domain.offset! multiply domain.offset in use
        let mut curr_coset_offset = fri_parameters.domain.gen().pow(&[curr_coset_index as u64]);
        let mut expected_round_eval = F::zero();

        debug_assert_eq!(fri_parameters.localization_parameters.len(), queries.len());
        for i in 0..queries.len() {
            expected_round_eval = FRIVerifier::expected_evaluation(&queries[i], queried_evaluations[i].clone(), alphas[i]);

            // check if current round result is consistent with next round polynomial (if next round is not final)
            if i < queries.len() - 1 {
                let next_localization_param = fri_parameters.localization_parameters[i + 1] as usize;
                let next_intra_coset_index = &curr_coset_index >> next_localization_param;

                let actual = queried_evaluations[i+1][next_intra_coset_index];
                if expected_round_eval != actual {
                    return false
                }
                curr_coset_index = curr_coset_index & ((1 << next_localization_param) - 1)
            }

            curr_coset_offset = curr_coset_offset.pow(&[1 << fri_parameters.localization_parameters[i] as u64]);
        }

        // check final polynomial (low degree & consistency check)
        // We assume degree_bound is power of 2.
        assert!(fri_parameters.tested_degree.is_power_of_two());
        let total_shrink_factor: u64 = fri_parameters.localization_parameters.iter().sum();
        let final_poly_degree_bound = fri_parameters.tested_degree >> total_shrink_factor;

        let final_ldt = DirectLDT::generate_low_degree_coefficients(final_polynomial_domain.clone(),
                                                                    final_polynomial.to_vec(),
                                                                    final_poly_degree_bound as usize);
        DirectLDT::verify_low_degree_single_round(curr_coset_offset, expected_round_eval, &final_ldt)
    }

    /// Map coset in current round to a single point in next round.
    ///
    /// Essentially, this function interpolates the polynomial and evaluate on `alpha`.
    #[inline]
    fn expected_evaluation(
        coset: &Radix2CosetDomain<F>,
        queried_evaluations: Vec<F>,
        alpha: F
    ) -> F {
        let poly = coset.interpolate(queried_evaluations);
        poly.evaluate(&alpha)
    }
}

fn le_bits_array_to_usize(bits: &[bool]) -> usize {
    let mut result = 0;
    for &bit in bits{
        result += bit as usize;
        result *= 2;
    }
    result
}