use crate::domain::Radix2CosetDomain;
use ark_ff::{batch_inversion_and_mul, PrimeField};
use ark_r1cs_std::poly::evaluations::univariate::lagrange_interpolator::LagrangeInterpolator;
use ark_std::marker::PhantomData;
use ark_std::vec::Vec;
/// FRI Prover
pub struct FRIProver<F: PrimeField> {
    _prover: PhantomData<F>,
}

impl<F: PrimeField> FRIProver<F> {
    /// Single round prover in commit phase. Returns the polynomial for next round
    /// represented by evaluations over domain in next round.
    ///
    /// The prover is inefficient. TODO: Adapt code from libiop.
    ///
    /// Returns domain for next round polynomial and evaluations over the domain.
    pub fn interactive_phase_single_round_naive(
        domain: Radix2CosetDomain<F>,
        poly_over_domain: Vec<F>,
        localization_param: u64,
        alpha: F,
    ) -> (Radix2CosetDomain<F>, Vec<F>) {
        let coset_size = 1 << localization_param;
        let domain_size = domain.base_domain.size;
        let dist_between_coset_elems = domain_size / coset_size;
        let mut new_evals = Vec::with_capacity(dist_between_coset_elems as usize);
        let coset_generator = domain
            .gen()
            .pow(&[1 << (domain.dim() as u64 - localization_param)]);
        let mut cur_coset_offset = domain.offset;

        for coset_index in 0..dist_between_coset_elems {
            let mut poly_evals = Vec::new();
            for intra_coset_index in 0..coset_size {
                poly_evals.push(
                    poly_over_domain
                        [(coset_index + intra_coset_index * dist_between_coset_elems) as usize],
                );
            }

            let interpolator = LagrangeInterpolator::new(
                cur_coset_offset,
                coset_generator,
                localization_param,
                poly_evals,
            );
            new_evals.push(interpolator.interpolate(alpha));
            cur_coset_offset *= domain.gen();
        }

        let c = Radix2CosetDomain::new_radix2_coset(new_evals.len(), domain.offset);
        // c.base_domain.group_gen = coset_generator;
        // c.base_domain.group_gen_inv = coset_generator.inverse().unwrap();
        debug_assert_eq!(coset_generator.pow(&[new_evals.len() as u64]), F::one());
        debug_assert_eq!(c.size(), new_evals.len());
        (c, new_evals)
    }

    /// Single round prover in commit phase. Returns the polynomial for next round
    /// represented by evaluations over domain in next round.
    ///
    /// Returns domain for next round polynomial and evaluations over the domain.
    pub fn interactive_phase_single_round(
        domain: Radix2CosetDomain<F>,
        evals_over_domain: Vec<F>,
        localization_param: u64,
        alpha: F,
    ) -> (Radix2CosetDomain<F>, Vec<F>) {
        let coset_size = 1 << localization_param;
        let num_cosets = domain.size() / coset_size;
        let mut next_f_i = Vec::with_capacity(num_cosets); // new_evals

        let h_inc = domain.gen();
        let h_inc_to_coset_inv_plus_one =
            h_inc.pow(&[coset_size as u64]).inverse().unwrap() * h_inc;

        let shiftless_coset = Radix2CosetDomain::new_radix2_coset(coset_size, F::one());
        let g = shiftless_coset.gen();
        let g_inv = g.inverse().unwrap();
        let x_to_order_coset = alpha.pow(&[coset_size as u64]);

        // x * g^{-k}
        let mut shifted_x_elements = Vec::with_capacity(coset_size);
        shifted_x_elements.push(alpha);
        for i in 1..coset_size {
            shifted_x_elements.push(shifted_x_elements[i - 1] * g_inv);
        }

        let mut cur_h = domain.offset;
        let first_h_to_coset_inv_plus_one =
            cur_h.pow(&[coset_size as u64]).inverse().unwrap() * cur_h;
        let mut cur_coset_constant_plus_h = x_to_order_coset * first_h_to_coset_inv_plus_one;

        /* x * g^{-k} - h, for all combinations of k, h.  */
        let mut elements_to_invert = Vec::with_capacity(evals_over_domain.len());

        /* constant for each coset, equal to
         *  vp_coset(x) / h^{|coset| - 1} = x^{|coset|} h^{-|coset| + 1} - h */
        let mut constant_for_each_coset = Vec::with_capacity(num_cosets);

        let constant_for_all_cosets = F::from(coset_size as u128).inverse().unwrap();
        let mut x_ever_in_domain = false;
        let mut x_coset_index = 0;
        let mut x_index_in_domain = 0;

        /* First we create all the constants for each coset,
          and the entire vector of elements to invert, xg^{-k} - h.
        */

        for j in 0..num_cosets {
            /* coset constant = x^|coset| * h^{1 - |coset|} - h */
            let coset_constant: F = cur_coset_constant_plus_h - cur_h;
            constant_for_each_coset.push(coset_constant);
            /* coset_constant = vp_coset(x) * h^{-|coset| + 1},
            since h is non-zero, coset_constant is zero iff vp_coset(x) is zero.
            If vp_coset(x) is zero, then x is in the coset. */
            let x_in_coset = coset_constant.is_zero();
            /* if x is in the coset, we mark which position x is within f_i_domain,
            and we pad elements to invert to simplify inversion later. */
            if x_in_coset {
                x_ever_in_domain = true;
                x_coset_index = j;
                // find which element in the coset x belongs to.
                // also pad elements_to_invert to simplify indexing
                let mut cur_elem = cur_h;
                for k in 0..coset_size {
                    if cur_elem == alpha {
                        x_index_in_domain = k * num_cosets + j;
                    }
                    cur_elem *= g;
                    elements_to_invert.push(F::one());
                }
                continue;
            }

            /* Append all elements to invert, (xg^{-k} - h) */
            for k in 0..coset_size {
                elements_to_invert.push(shifted_x_elements[k] - cur_h);
            }

            cur_h *= h_inc;
            /* coset constant = x^|coset| * h^{1 - |coset|} - h
            So we can efficiently increment x^|coset| * h^{1 - |coset|} */
            cur_coset_constant_plus_h *= h_inc_to_coset_inv_plus_one;
        }
        /* Technically not lagrange coefficients, its missing the constant for each coset */
        batch_inversion_and_mul(&mut elements_to_invert, &constant_for_all_cosets);
        let lagrange_coefficients = elements_to_invert;
        for j in 0..num_cosets {
            let mut interpolation = F::zero();
            for k in 0..coset_size {
                interpolation += evals_over_domain[k * num_cosets + j]
                    * lagrange_coefficients[j * coset_size + k];
            }
            /* Multiply the constant for each coset, to get the correct interpolation */
            interpolation *= constant_for_each_coset[j];
            next_f_i.push(interpolation);
        }

        /* if x ever in domain, correct that evaluation. */
        if x_ever_in_domain {
            next_f_i[x_coset_index] = evals_over_domain[x_index_in_domain];
        }

        // domain definition
        let c = Radix2CosetDomain::new_radix2_coset(next_f_i.len(), domain.offset);

        (c, next_f_i)
    }
}

#[cfg(test)]
pub mod tests {
    use crate::direct::DirectLDT;
    use crate::domain::Radix2CosetDomain;
    use crate::fri::prover::FRIProver;
    use crate::fri::verifier::FRIVerifier;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::UVPolynomial;
    use ark_std::{test_rng, UniformRand};
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn efficient_prover_consistency_test() {
        let degree = 32;

        let mut rng = test_rng();
        let poly = DensePolynomial::<Fr>::rand(degree, &mut rng);
        let domain_coset = Radix2CosetDomain::new_radix2_coset(64, Fr::rand(&mut rng));
        let evaluations = domain_coset.evaluate(&poly);

        // fri prover should reduce its degree
        let alpha = Fr::rand(&mut rng);
        let localization = 2;
        let (expected_domain_next_round, expected_eval_next_round) =
            FRIProver::interactive_phase_single_round_naive(
                domain_coset,
                evaluations.to_vec(),
                localization,
                alpha,
            );

        let (actual_domain_next_round, actual_eval_next_round) =
            FRIProver::interactive_phase_single_round(
                domain_coset,
                evaluations.to_vec(),
                localization,
                alpha,
            );

        assert_eq!(actual_domain_next_round, expected_domain_next_round);
        assert_eq!(actual_eval_next_round, expected_eval_next_round);
    }

    #[test]
    fn degree_reduction_test() {
        let degree = 32;

        let mut rng = test_rng();
        let poly = DensePolynomial::<Fr>::rand(degree, &mut rng);
        let domain_coset = Radix2CosetDomain::new_radix2_coset(64, Fr::rand(&mut rng));
        let evaluations = domain_coset.evaluate(&poly);

        // fri prover should reduce its degree
        let alpha = Fr::rand(&mut rng);
        let localization = 2;
        let (domain_next_round, eval_next_round) = FRIProver::interactive_phase_single_round(
            domain_coset.clone(),
            evaluations.to_vec(),
            localization,
            alpha,
        );

        let low_degree_poly = DirectLDT::generate_low_degree_coefficients(
            domain_next_round.clone(),
            eval_next_round.to_vec(),
            degree / (1 << localization),
        );

        let sampled_element = domain_next_round.element(15);
        let sampled_evaluation = eval_next_round[15];

        assert!(DirectLDT::verify_low_degree_single_round(
            sampled_element,
            sampled_evaluation,
            &low_degree_poly
        ));

        // test `fold_domain`
        let fold_domain = FRIVerifier::fold_domain(domain_coset.clone(), localization);
        assert_eq!(fold_domain, domain_next_round);
    }
}
