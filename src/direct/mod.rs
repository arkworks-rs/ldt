/// R1CS constraints for DirectLDT
#[cfg(feature = "r1cs")]
pub mod constraints;

use crate::domain::Radix2CosetDomain;
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::Polynomial;
use ark_std::marker::PhantomData;
use ark_std::vec::Vec;
/// Direct LDT by interpolating evaluations and truncating coefficients to low degree.
pub struct DirectLDT<F: PrimeField> {
    marker: PhantomData<F>,
}

/// A linear-communication protocol for testing if a function is a polynomial of certain degree.
/// Method is described in Aurora appendix C.1.
///
/// For now, the domain of the function needs to support IFFT.
impl<F: PrimeField> DirectLDT<F> {
    /// ### Prover Side
    ///
    /// Generate the coefficient of the low-degree polynomial obtained by interpolating the domain evaluations.
    /// The polynomial is trimmed to `degree_bound` as necessary.
    pub fn generate_low_degree_coefficients(
        domain: Radix2CosetDomain<F>,
        evaluations: Vec<F>,
        degree_bound: usize,
    ) -> DensePolynomial<F> {
        let mut poly = domain.interpolate(evaluations);
        // trim higher degree: if poly is higher degree, then the soundness should fail
        poly.coeffs.truncate(degree_bound + 1);
        poly
    }

    /// ### Verifier Side
    ///
    /// Verifier sample one element from domain and get its evaluation. Check if that evaluation
    /// agrees with low-degree polynomial.
    pub fn verify_low_degree_single_round(
        sampled_domain_element: F,
        sampled_evaluation_element: F,
        coefficients: &DensePolynomial<F>,
    ) -> bool {
        return coefficients.evaluate(&sampled_domain_element) == sampled_evaluation_element;
    }
}

#[cfg(test)]
mod tests {
    use crate::direct::{DirectLDT, Radix2CosetDomain};
    use ark_ff::UniformRand;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, Polynomial, Radix2EvaluationDomain, UVPolynomial};
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_r1cs_std::fields::FieldVar;
    use ark_r1cs_std::poly::evaluations::univariate::EvaluationsVar;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::test_rng;
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_native_coset() {
        let mut rng = test_rng();
        let degree = 51;
        let poly = DensePolynomial::<Fr>::rand(degree, &mut rng);
        let base_domain = Radix2EvaluationDomain::new(degree + 1).unwrap();
        let offset = Fr::rand(&mut rng);
        let coset = Radix2CosetDomain::new(base_domain, offset);

        // test evaluation
        let expected_eval: Vec<_> = coset
            .base_domain
            .elements()
            .map(|x| poly.evaluate(&(offset * x)))
            .collect();
        let actual_eval = coset.evaluate(&poly);
        assert_eq!(actual_eval, expected_eval);

        // test interpolation
        let interpolated_poly = coset.interpolate(expected_eval.to_vec());
        assert_eq!(interpolated_poly, poly);

        // test consistency with r1cs-std
        let cs = ConstraintSystem::new_ref();
        let eval_var: Vec<_> = expected_eval
            .iter()
            .map(|x| FpVar::new_witness(ark_relations::ns!(cs, "eval_var"), || Ok(*x)).unwrap())
            .collect();
        let r1cs_coset = ark_r1cs_std::poly::domain::Radix2DomainVar {
            gen: base_domain.group_gen,
            offset: FpVar::constant(offset),
            dim: ark_std::log2(degree.next_power_of_two()) as u64,
        };
        let eval_var = EvaluationsVar::from_vec_and_domain(eval_var, r1cs_coset, true);

        let pt = Fr::rand(&mut rng);
        let pt_var = FpVar::new_witness(ark_relations::ns!(cs, "random point"), || Ok(pt)).unwrap();

        let expected = poly.evaluate(&pt);
        let actual = eval_var.interpolate_and_evaluate(&pt_var).unwrap();

        assert_eq!(actual.value().unwrap(), expected);
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn test_direct_ldt() {
        let degree = 51;

        let mut rng = test_rng();
        let poly = DensePolynomial::<Fr>::rand(degree, &mut rng);
        let domain_coset = Radix2CosetDomain::new_radix2_coset(52, Fr::rand(&mut rng));
        let evaluations = domain_coset.evaluate(&poly);

        let low_degree_poly = DirectLDT::generate_low_degree_coefficients(
            domain_coset.clone(),
            evaluations.to_vec(),
            degree,
        );

        let sampled_element = domain_coset.element(15);
        let sampled_evaluation = evaluations[15];

        assert!(DirectLDT::verify_low_degree_single_round(
            sampled_element,
            sampled_evaluation,
            &low_degree_poly
        ))
    }
}
