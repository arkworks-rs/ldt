#![allow(unused)] // temporary
use crate::direct::constraints::DirectLDTGadget;
use crate::domain::Radix2CosetDomain;
use crate::fri::FRIParameters;
use ark_crypto_primitives::sponge::constraints::CryptographicSpongeVar;
use ark_crypto_primitives::sponge::FieldBasedCryptographicSponge;
use ark_ff::PrimeField;
use ark_r1cs_std::boolean::*;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::poly::domain::Radix2DomainVar;
use ark_r1cs_std::poly::evaluations::univariate::EvaluationsVar;
use ark_r1cs_std::poly::polynomial::univariate::dense::DensePolynomialVar;
use ark_r1cs_std::prelude::CondSelectGadget;
use ark_relations::r1cs::SynthesisError;
use ark_std::marker::PhantomData;
use ark_std::vec::*;

/// Constraints for FRI verifier.
pub struct FRIVerifierGadget<F: PrimeField> {
    _marker: PhantomData<F>,
}

impl<F: PrimeField> FRIVerifierGadget<F> {
    /// ## Step 1: Interative Phase
    /// Sample alpha in interactive phase.
    pub fn interactive_phase_single_round<
        S: FieldBasedCryptographicSponge<F>,
        SV: CryptographicSpongeVar<F, S>,
    >(
        sponge_var: &mut SV,
    ) -> Result<FpVar<F>, SynthesisError> {
        Ok(sponge_var
            .squeeze_field_elements(1)?
            .first()
            .unwrap()
            .clone())
    }

    /// ## Step 2: Sample Queried Coset
    /// Sample the coset to be queried.
    pub fn sample_coset_index<
        S: FieldBasedCryptographicSponge<F>,
        SV: CryptographicSpongeVar<F, S>,
    >(
        sponge_var: &mut SV,
        fri_parameters: &FRIParameters<F>,
    ) -> Result<Vec<Boolean<F>>, SynthesisError> {
        let log_num_cosets =
            fri_parameters.domain.dim() - fri_parameters.localization_parameters[0] as usize;
        sponge_var.squeeze_bits(log_num_cosets)
    }

    /// ## Step 2: Query Phase (Prepare Query)
    /// Prepare one query given the random coset index. The returned value `queries[i]` is the coset query
    /// of the `ith` round polynomial (including codeword but does not include final polynomial).
    /// Final polynomial is not queried. Instead, verifier will get
    /// the whole final polynomial in evaluation form, and do direct LDT.
    ///
    /// Returns the all query domains, and query coset index, final polynomial domain
    pub fn prepare_query(
        rand_coset_index: Vec<Boolean<F>>,
        fri_parameters: &FRIParameters<F>,
    ) -> Result<
        (
            Vec<Radix2DomainVar<F>>,
            Vec<Vec<Boolean<F>>>,
            Radix2CosetDomain<F>,
        ),
        SynthesisError,
    > {
        let num_fri_rounds = fri_parameters.localization_parameters.len();
        let mut coset_indices = Vec::new();
        let mut curr_coset_index = rand_coset_index;
        let mut queries = Vec::with_capacity(num_fri_rounds);
        let mut curr_round_domain = fri_parameters.domain;

        // sample coset index
        for i in 0..num_fri_rounds {
            let log_dist_between_coset_elems =
                curr_round_domain.dim() - fri_parameters.localization_parameters[i] as usize;
            curr_coset_index = curr_coset_index[..log_dist_between_coset_elems].to_vec();

            coset_indices.push(curr_coset_index.clone());

            // get the query coset from coset index
            let query_gen = fri_parameters.domain.gen().pow(&[1
                << (fri_parameters.domain.dim()
                    - fri_parameters.localization_parameters[i] as usize)]);
            debug_assert_eq!(
                query_gen.pow(&[1 << fri_parameters.localization_parameters[i]]),
                F::one()
            );

            let query_offset = &FpVar::constant(curr_round_domain.offset)
                * &(FpVar::constant(curr_round_domain.gen()).pow_le(&curr_coset_index)?);

            let query_coset = Radix2DomainVar::new(
                query_gen,
                fri_parameters.localization_parameters[i],
                query_offset,
            )?;

            queries.push(query_coset);

            curr_round_domain = curr_round_domain.fold(fri_parameters.localization_parameters[i])
        }

        Ok((queries, coset_indices, curr_round_domain))
    }

    /// Map coset in current round to a single point in next round.
    ///
    /// Essentially, this function interpolates the polynomial and evaluate on `alpha`.
    fn expected_evaluation(
        coset: &Radix2DomainVar<F>,
        queried_evaluations: Vec<FpVar<F>>,
        alpha: FpVar<F>,
    ) -> Result<FpVar<F>, SynthesisError> {
        let evaluations =
            EvaluationsVar::from_vec_and_domain(queried_evaluations, coset.clone(), true);
        evaluations.interpolate_and_evaluate(&alpha)
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
        queried_coset_indices: &[Vec<Boolean<F>>],
        queries: &[Radix2DomainVar<F>],
        queried_evaluations: &[Vec<FpVar<F>>],
        alphas: &[FpVar<F>],
        final_polynomial_domain: &Radix2CosetDomain<F>,
        final_polynomial: &DensePolynomialVar<F>,
    ) -> Result<Boolean<F>, SynthesisError> {
        let mut expected_next_round_eval = FpVar::zero();

        debug_assert_eq!(fri_parameters.localization_parameters.len(), queries.len());
        let mut check_result = Boolean::constant(true);
        for i in 0..queries.len() {
            expected_next_round_eval = FRIVerifierGadget::expected_evaluation(
                &queries[i],
                queried_evaluations[i].clone(),
                alphas[i].clone(),
            )?;

            // check if current round result is consistent with next round polynomial (if next round is not final)
            if i < queries.len() - 1 {
                let next_localization_param =
                    fri_parameters.localization_parameters[i + 1] as usize;
                let log_next_dist_between_coset_elems =
                    fri_parameters.log_round_coset_sizes[i + 1] - next_localization_param;
                // native code: queried_coset_indices[i] >> log_next_dist_between_coset_elems
                let next_intra_coset_index =
                    &queried_coset_indices[i][log_next_dist_between_coset_elems..];

                let actual = FpVar::<F>::conditionally_select_power_of_two_vector(
                    next_intra_coset_index,
                    &queried_evaluations[i + 1],
                )?;

                check_result = Boolean::<F>::kary_and(&[
                    check_result,
                    expected_next_round_eval.is_eq(&actual)?,
                ])
                .unwrap();
            }
        }

        // check final polynomial (low degree & consistency check)
        // We assume degree_bound is power of 2.
        assert!(fri_parameters.tested_degree.is_power_of_two());
        let total_shrink_factor: u64 = fri_parameters.localization_parameters.iter().sum();
        let final_poly_degree_bound = fri_parameters.tested_degree >> total_shrink_factor;

        let final_element_index = queried_coset_indices.last().unwrap();

        DirectLDTGadget::verify(
            final_polynomial_domain.element_var(final_element_index)?,
            expected_next_round_eval,
            final_polynomial,
            final_poly_degree_bound as usize,
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::direct::DirectLDT;
    use crate::domain::Radix2CosetDomain;
    use crate::fri::constraints::FRIVerifierGadget;
    use crate::fri::prover::FRIProver;
    use crate::fri::verifier::FRIVerifier;
    use crate::fri::FRIParameters;
    use ark_poly::polynomial::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::convert::ToBitsGadget;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::poly::polynomial::univariate::dense::DensePolynomialVar;
    use ark_r1cs_std::uint64::UInt64;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_relations::*;
    use ark_std::{test_rng, UniformRand};
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_prepare_query() {
        let mut rng = test_rng();
        let offset = Fr::rand(&mut rng);
        let domain_input = Radix2CosetDomain::new_radix2_coset(1 << 7, offset);

        let fri_parameters = FRIParameters::new(32, vec![1, 2, 1], domain_input);

        let rand_coset_index = 31usize;
        let cs = ConstraintSystem::new_ref();
        let rand_coset_index_var =
            UInt64::new_witness(ns!(cs, "rand_coset_index"), || Ok(rand_coset_index as u64))
                .unwrap();
        let rand_coset_index_var_arr =
            rand_coset_index_var.to_bits_le().unwrap()[..(1 << 6)].to_vec();

        let rand_coset_index = 31;
        let (query_cosets, query_indices, domain_final) =
            FRIVerifier::prepare_query(rand_coset_index, &fri_parameters);
        let (query_cosets_actual, query_indices_actual, domain_final_actual) =
            FRIVerifierGadget::prepare_query(rand_coset_index_var_arr, &fri_parameters).unwrap();

        for i in 0..query_cosets.len() {
            assert_eq!(
                query_cosets_actual[i].offset().value().unwrap(),
                query_cosets[i].offset
            );
            assert_eq!(query_cosets_actual[i].gen, query_cosets[i].gen());
            assert_eq!(query_cosets_actual[i].dim as usize, query_cosets[i].dim());
        }

        assert_eq!(domain_final, domain_final_actual)
    }

    #[test]
    fn two_rounds_fri_test() {
        let cs = ConstraintSystem::new_ref();

        let mut rng = test_rng();
        let poly = DensePolynomial::rand(64, &mut rng);
        let offset = Fr::rand(&mut rng);
        let domain_input = Radix2CosetDomain::new_radix2_coset(128, offset);
        let evaluations_input = domain_input.evaluate(&poly);

        // set up verifier parameters
        let fri_parameters = FRIParameters::new(64, vec![1, 2, 2], domain_input);
        let alphas: Vec<_> = (0..3).map(|_| Fr::rand(&mut rng)).collect();
        let alphas_var: Vec<_> = alphas
            .iter()
            .map(|x| FpVar::new_witness(ns!(cs, "alphas"), || Ok(x)).unwrap())
            .collect();

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

        let rand_coset_index = 31;
        let rand_coset_index_var =
            UInt64::new_witness(ns!(cs, "rand_coset_index"), || Ok(rand_coset_index))
                .unwrap()
                .to_bits_le();

        let (query_cosets, query_indices, domain_final) =
            FRIVerifierGadget::prepare_query(rand_coset_index_var.unwrap(), &fri_parameters)
                .unwrap();
        let (_, query_indices_native, _) =
            FRIVerifier::prepare_query(rand_coset_index as usize, &fri_parameters);

        assert_eq!(query_indices.len(), 3);
        assert_eq!(domain_final, expected_domain_final);

        let (indices, qi) = domain_input.query_position_to_coset(
            query_indices_native[0],
            fri_parameters.localization_parameters[0] as usize,
        );
        assert_eq!(qi.offset, query_cosets[0].offset().value().unwrap());
        let answer_input: Vec<_> = indices
            .iter()
            .map(|&i| {
                FpVar::new_witness(ns!(cs, "answer_input"), || Ok(evaluations_input[i])).unwrap()
            })
            .collect();

        let (indices, q0) = domain_round_0.query_position_to_coset(
            query_indices_native[1],
            fri_parameters.localization_parameters[1] as usize,
        );
        assert_eq!(q0.offset, query_cosets[1].offset().value().unwrap());
        let answer_round_0: Vec<_> = indices
            .iter()
            .map(|&i| {
                FpVar::new_witness(
                    ns!(cs, "evaluations_round_0"),
                    || Ok(evaluations_round_0[i]),
                )
                .unwrap()
            })
            .collect();

        let (indices, q1) = domain_round_1.query_position_to_coset(
            query_indices_native[2],
            fri_parameters.localization_parameters[2] as usize,
        );
        let answer_round_1: Vec<_> = indices
            .iter()
            .map(|&i| {
                FpVar::new_witness(
                    ns!(cs, "evaluations_round_1"),
                    || Ok(evaluations_round_1[i]),
                )
                .unwrap()
            })
            .collect();
        assert_eq!(q1.offset, query_cosets[2].offset().value().unwrap());

        let total_shrink_factor: u64 = fri_parameters.localization_parameters.iter().sum();
        let final_poly_degree_bound = fri_parameters.tested_degree >> total_shrink_factor;
        let final_polynomial = DirectLDT::generate_low_degree_coefficients(
            domain_final,
            evaluations_final,
            final_poly_degree_bound as usize,
        );
        let final_polynomial_coeffs: Vec<_> = final_polynomial
            .coeffs()
            .iter()
            .map(|x| FpVar::new_witness(ns!(cs, "final_poly_coeff"), || Ok(*x)).unwrap())
            .collect();
        let final_polynomial_var =
            DensePolynomialVar::from_coefficients_slice(&final_polynomial_coeffs);

        let result = FRIVerifierGadget::consistency_check(
            &fri_parameters,
            &query_indices,
            &query_cosets,
            &vec![answer_input, answer_round_0, answer_round_1],
            &alphas_var,
            &domain_final,
            &final_polynomial_var,
        )
        .unwrap();

        assert!(result.value().unwrap());
        assert!(cs.is_satisfied().unwrap());
    }
}
