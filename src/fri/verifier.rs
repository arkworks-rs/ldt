use ark_std::marker::PhantomData;

use crate::direct::DirectLDT;
use crate::domain::Radix2CosetDomain;
use crate::fri::FRIParameters;
use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::Polynomial;
use ark_sponge::FieldBasedCryptographicSponge;
use ark_std::vec::Vec;

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
        fri_parameters: &FRIParameters<F>,
    ) -> usize {
        let log_num_cosets =
            fri_parameters.domain.dim() - fri_parameters.localization_parameters[0] as usize;
        // we use the fact that number of cosets is always power of two
        let rand_coset_index = le_bits_array_to_usize(&sponge.squeeze_bits(log_num_cosets));
        rand_coset_index
    }

    /// ## Step 2: Query Phase (Prepare Query)
    /// Prepare one query given the random coset index. The returned value `queries[i]` is the coset query
    /// of the `ith` round polynomial (including codeword but does not include final polynomial).
    /// Final polynomial is not queried. Instead, verifier will get
    /// the whole final polynomial in evaluation form, and do direct LDT.
    ///
    /// Returns the all query domains, and query coset index, final polynomial domain
    pub fn prepare_query(
        rand_coset_index: usize,
        fri_parameters: &FRIParameters<F>,
    ) -> (Vec<Radix2CosetDomain<F>>, Vec<usize>, Radix2CosetDomain<F>) {
        let num_fri_rounds = fri_parameters.localization_parameters.len();
        let mut coset_indices = Vec::new();
        let mut curr_coset_index = rand_coset_index;
        let mut queries = Vec::with_capacity(num_fri_rounds);
        let mut curr_round_domain = fri_parameters.domain;
        // sample a coset index
        for i in 0..num_fri_rounds {
            // current coset index = last coset index % (distance between coset at current round)
            // edge case: at first round, this still applies

            let dist_between_coset_elems =
                curr_round_domain.size() / (1 << fri_parameters.localization_parameters[i]);
            curr_coset_index = curr_coset_index % dist_between_coset_elems;

            coset_indices.push(curr_coset_index);

            let (_, query_coset) = curr_round_domain.query_position_to_coset(
                curr_coset_index,
                fri_parameters.localization_parameters[i] as usize,
            );

            queries.push(query_coset);

            // get next round coset size, and next round domain
            curr_round_domain = curr_round_domain.fold(fri_parameters.localization_parameters[i]);
        }

        (queries, coset_indices, curr_round_domain)
    }

    /// ## Step 2: Query Phase (Prepare Query)
    /// Prepare all queries given the sampled random coset indices.
    ///
    /// The first returned value `queries[i][j]` is the coset query
    /// of the `j`th round polynomial (including codeword but does not include final polynomial) for `i`th query.
    ///
    /// The second returned value `indices[i][j]` is the coset index
    /// of the `j`th round polynomial (including codeword but does not include final polynomial) for `i`th query.
    ///
    /// The last returned value `final[i]` is the final polynomial domain at round `i`.
    pub fn batch_prepare_queries(
        rand_coset_indices: &[usize],
        fri_parameters: &FRIParameters<F>,
    ) -> (
        Vec<Vec<Radix2CosetDomain<F>>>,
        Vec<Vec<usize>>,
        Vec<Radix2CosetDomain<F>>,
    ) {
        let mut queries = Vec::with_capacity(rand_coset_indices.len());
        let mut indices = Vec::with_capacity(rand_coset_indices.len());
        let mut finals = Vec::with_capacity(rand_coset_indices.len());

        rand_coset_indices
            .iter()
            .map(|&i| Self::prepare_query(i, fri_parameters))
            .for_each(|(query, index, fp)| {
                queries.push(query);
                indices.push(index);
                finals.push(fp);
            });

        (queries, indices, finals)
    }

    /// ## Step 3: Decision Phase (Check query)
    /// After preparing the query, verifier get the evaluations of corresponding query. Those evaluations needs
    /// to be checked by merkle tree. Then verifier calls this method to check if polynomial sent in each round
    /// is consistent with each other, and the final polynomial is low-degree.
    ///
    /// `queries[i]` is the coset query of the `ith` round polynomial, including the codeword polynomial.
    /// `queried_evaluations` stores the result of corresponding query.
    pub fn consistency_check(
        fri_parameters: &FRIParameters<F>,
        queried_coset_indices: &[usize],
        queries: &[Radix2CosetDomain<F>],
        queried_evaluations: &[Vec<F>],
        alphas: &[F],
        final_polynomial_domain: &Radix2CosetDomain<F>,
        final_polynomial: &DensePolynomial<F>,
    ) -> bool {
        let mut expected_next_round_eval = F::zero();

        debug_assert_eq!(fri_parameters.localization_parameters.len(), queries.len());
        for i in 0..queries.len() {
            expected_next_round_eval = FRIVerifier::expected_evaluation(
                &queries[i],
                queried_evaluations[i].clone(),
                alphas[i],
            );

            // check if current round result is consistent with next round polynomial (if next round is not final)
            if i < queries.len() - 1 {
                let next_localization_param =
                    fri_parameters.localization_parameters[i + 1] as usize;
                let log_next_dist_between_coset_elems =
                    fri_parameters.log_round_coset_sizes[i + 1] - next_localization_param;
                let next_intra_coset_index =
                    queried_coset_indices[i] >> log_next_dist_between_coset_elems;

                let actual = queried_evaluations[i + 1][next_intra_coset_index];
                if expected_next_round_eval != actual {
                    return false;
                }
            }
        }

        // check final polynomial (low degree & consistency check)
        // We assume degree_bound is power of 2.
        assert!(fri_parameters.tested_degree.is_power_of_two());
        let total_shrink_factor: u64 = fri_parameters.localization_parameters.iter().sum();
        let final_poly_degree_bound = fri_parameters.tested_degree >> total_shrink_factor;

        let final_element_index = *queried_coset_indices.last().unwrap();

        assert!(
            final_polynomial.degree() <= final_poly_degree_bound as usize,
            "final polynomial degree is too large!"
        );
        DirectLDT::verify_low_degree_single_round(
            final_polynomial_domain.element(final_element_index),
            expected_next_round_eval,
            &final_polynomial,
        )
    }

    /// ## Step 3: Decision Phase (Check query)
    /// After preparing all queries, verifier gets the evaluations of corresponding query. Those evaluations needs
    /// to be checked by merkle tree. Then verifier calls this method to check if polynomial sent in each round
    /// is consistent with each other, and the final polynomial is low-degree.
    ///
    /// * `all_queried_coset_indices[i][j]` is the `j`th round query coset index of `i`th query
    /// * `all_queries_domains[i][j]` is the `j`th round query coset of `i`th query
    /// * `all_queried_evaluations[i][j]` is a vector storing corresponding evaluations at `all_queries_domains[i][j]`
    /// * `alphas[i]` is the randomness used by the polynomial
    /// * `all_final_polynomial_domain[i]` is the final polynomial domain for `i`th query
    /// * `all_final_polynomials` is the final polynomial for `i`th query
    pub fn batch_consistency_check(
        fri_parameters: &FRIParameters<F>,
        all_queried_coset_indices: &[Vec<usize>],
        all_queries_domains: &[Vec<Radix2CosetDomain<F>>],
        all_queried_evaluations: &[Vec<Vec<F>>],
        alphas: &[F],
        all_final_polynomial_domain: &[Radix2CosetDomain<F>],
        all_final_polynomials: &[DensePolynomial<F>],
    ) -> bool {
        for i in 0..all_queried_coset_indices.len() {
            let result = Self::consistency_check(
                fri_parameters,
                &all_queried_coset_indices[i],
                &all_queries_domains[i],
                &all_queried_evaluations[i],
                alphas,
                &all_final_polynomial_domain[i],
                &all_final_polynomials[i],
            );
            if !result {
                return false;
            }
        }
        true
    }

    /// Map coset in current round to a single point in next round.
    ///
    /// Essentially, this function interpolates the polynomial and evaluate on `alpha`.
    #[inline]
    fn expected_evaluation(
        coset: &Radix2CosetDomain<F>,
        queried_evaluations: Vec<F>,
        alpha: F,
    ) -> F {
        let poly = coset.interpolate(queried_evaluations);
        poly.evaluate(&alpha)
    }
}

fn le_bits_array_to_usize(bits: &[bool]) -> usize {
    let mut result = 0;
    for &bit in bits {
        result += bit as usize;
        result *= 2;
    }
    result
}

#[cfg(test)]
mod tests {
    use ark_ff::UniformRand;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{Polynomial, UVPolynomial};
    use ark_std::test_rng;
    use ark_test_curves::bls12_381::Fr;

    use crate::direct::DirectLDT;
    use crate::domain::Radix2CosetDomain;
    use crate::fri::prover::FRIProver;
    use crate::fri::verifier::FRIVerifier;
    use crate::fri::FRIParameters;

    #[test]
    fn two_rounds_fri_test() {
        // First, generate a low degree polynomial, and its evaluations.
        let mut rng = test_rng();
        let poly = DensePolynomial::rand(32, &mut rng);
        let offset = Fr::rand(&mut rng);
        let domain_input = Radix2CosetDomain::new_radix2_coset(128, offset);
        let evaluations_input = domain_input.evaluate(&poly);

        // Set up verifier parameter
        let fri_parameters = FRIParameters::new(32, vec![1, 2, 1], domain_input);
        let alphas: Vec<_> = (0..3).map(|_| Fr::rand(&mut rng)).collect();

        // prover commits all round polynomial
        let (domain_round_0, evaluations_round_0) = FRIProver::interactive_phase_single_round(
            domain_input,
            evaluations_input.clone(),
            fri_parameters.localization_parameters[0],
            alphas[0],
        );

        let (domain_round_1, evaluations_round_1) = FRIProver::interactive_phase_single_round(
            domain_round_0,
            evaluations_round_0.clone(),
            fri_parameters.localization_parameters[1],
            alphas[1],
        );

        let (expected_domain_final, evaluations_final) = FRIProver::interactive_phase_single_round(
            domain_round_1,
            evaluations_round_1.clone(),
            fri_parameters.localization_parameters[2],
            alphas[2],
        );

        // verifier prepare queries
        let rand_coset_index = 31;
        let (query_cosets, query_indices, domain_final) =
            FRIVerifier::prepare_query(rand_coset_index, &fri_parameters);
        assert_eq!(query_indices.len(), 3);
        assert_eq!(domain_final, expected_domain_final);

        // prover generate answers to queries
        let (indices, qi) = domain_input.query_position_to_coset(
            query_indices[0],
            fri_parameters.localization_parameters[0] as usize,
        );
        let answer_input: Vec<_> = indices.iter().map(|&i| evaluations_input[i]).collect();
        assert_eq!(qi, query_cosets[0]);

        let (indices, q0) = domain_round_0.query_position_to_coset(
            query_indices[1],
            fri_parameters.localization_parameters[1] as usize,
        );
        let answer_round_0: Vec<_> = indices.iter().map(|&i| evaluations_round_0[i]).collect();
        assert_eq!(q0, query_cosets[1]);

        let (indices, q1) = domain_round_1.query_position_to_coset(
            query_indices[2],
            fri_parameters.localization_parameters[2] as usize,
        );
        let answer_round_1: Vec<_> = indices.iter().map(|&i| evaluations_round_1[i]).collect();
        assert_eq!(q1, query_cosets[2]);

        // sanity check: answer_round_i interpolate version contained in answer_round_i+1
        assert!(answer_round_0.contains(&qi.interpolate(answer_input.clone()).evaluate(&alphas[0])));
        assert!(
            answer_round_1.contains(&q0.interpolate(answer_round_0.clone()).evaluate(&alphas[1]))
        );
        assert!(evaluations_final
            .contains(&q1.interpolate(answer_round_1.clone()).evaluate(&alphas[2])));

        let total_shrink_factor: u64 = fri_parameters.localization_parameters.iter().sum();
        let final_poly_degree_bound = fri_parameters.tested_degree >> total_shrink_factor;
        let final_polynomial = DirectLDT::generate_low_degree_coefficients(
            domain_final,
            evaluations_final,
            final_poly_degree_bound as usize,
        );

        // verifier verifies consistency
        let result = FRIVerifier::consistency_check(
            &fri_parameters,
            &query_indices,
            &query_cosets,
            &vec![answer_input, answer_round_0, answer_round_1],
            &alphas,
            &domain_final,
            &final_polynomial,
        );

        assert!(result)
    }
}
