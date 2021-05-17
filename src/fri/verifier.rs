use crate::fri::FRIParameters;
use ark_ff::PrimeField;
use ark_sponge::{FieldBasedCryptographicSponge, FieldElementSize};
use std::marker::PhantomData;
use crate::domain::Radix2CosetDomain;
use ark_poly::Polynomial;

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
    /// of the `ith` round polynomial. Final polynomial is not queried. Instead, verifier will get
    /// the whole final polynomial in evaluation form, and do direct LDT.
    pub fn prepare_queries(
        rand_coset_index: usize,
        fri_parameters: &FRIParameters<F>,
    ) -> Vec<Radix2CosetDomain<F>> {
        let num_fri_rounds = fri_parameters.localization_parameters.len();

        // this does not include domain.offset! multiply domain.offset in use
        let mut curr_gen_offset = fri_parameters.domain.gen().pow(&[rand_coset_index as u64]);
        let mut queries = Vec::with_capacity(num_fri_rounds);
        // sample a coset index
        for i in 0..num_fri_rounds {
            let coset_size = 1 << fri_parameters.localization_parameters[i];
            let coset_offset = fri_parameters.domain.offset * curr_gen_offset;

            // // TODO: debug only. seems that it's not necessary
            // // coset_generator = root of unity of order coset_size
            // // for L's generator g, it equals g^{1 << (L.dim - localization_parameter)}
            // let coset_generator = fri_parameters.domain.gen.pow(&[
            //     1 << (fri_parameters.domain.dim - fri_parameters.localization_parameters[i])
            // ]);

            let mut c = Radix2CosetDomain::new_radix2_coset(coset_size, coset_offset);
            c.base_domain.group_gen = fri_parameters.domain.gen().pow(&[
                    1 << (fri_parameters.domain.dim() - fri_parameters.localization_parameters[i] as usize)]);
            debug_assert_eq!(c.base_domain.group_gen.pow(&[coset_size as u64]), F::one());
            queries.push(c);

            curr_gen_offset = curr_gen_offset.pow(&[coset_size as u64]); // set next coset where current current folded coset is in
        }

        queries
    }

    /// ## Step 3: Query Phase (Check query)
    /// After preparing the query, verifier get the evaluations of corresponding query. Those evaluations needs
    /// to be checked by merkle tree. Then verifier calls this method to check if polynomial sent in each round
    /// is consistent with each other, and the final polynomial is low-degree.
    ///
    /// `queries[i]` is the coset query of the `ith` round polynomial. `queried_evaluations` stores the result of
    /// corresponding query.
    pub fn consistency_check(
        fri_parameters: &FRIParameters<F>,
        queried_coset_index: usize,
        queries: &[Radix2CosetDomain<F>],
        queried_evaluations: &[Vec<F>],
        alphas: &[F],
        final_polynomial: &[F]
    ) {
        let mut curr_coset_index = queried_coset_index;
        let mut curr_coset_offset = fri_parameters.domain.gen().pow(&[curr_coset_index as u64]);
        let mut expected_round_eval = F::zero();

        for i in 0..queries.len(){
            expected_round_eval = FRIVerifier::expected_evaluation(&queries[i], queried_evaluations[i].clone(), alphas[i]);

            // check if current round result is consistent with next round polynomial
            if i != queries.len() - 1 {

            }
        }
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