use ark_ff::PrimeField;
use std::marker::PhantomData;
use ark_poly::{EvaluationDomain, Evaluations, Polynomial};
use ark_poly::univariate::DensePolynomial;

pub struct DirectLDT<F: PrimeField>{
    marker: PhantomData<F>
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
    pub fn generate_low_degree_coefficients<D: EvaluationDomain<F>>(domain: D,
                                                                    evaluations: Vec<F>,
                                                                    degree_bound: usize) -> DensePolynomial<F> {
       let mut evals = Evaluations::from_vec_and_domain(evaluations, domain);
        let mut poly = evals.interpolate();
        // trim higher degree
        poly.coeffs.truncate(degree_bound+1);
        poly
    }

    /// ### Verifier Side
    ///
    /// Verifier sample one element from domain and get its evaluation. Check if that evaluation
    /// agrees with low-degree polynomial.
    pub fn verify_low_degree_single_round(sampled_domain_element: F,
                                          sampled_evaluation_element: F,
                                          coefficients: &DensePolynomial<F>) -> bool {
        return coefficients.evaluate(&sampled_domain_element) == sampled_evaluation_element;
    }

}

#[cfg(test)]
mod tests{
    use ark_std::test_rng;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{UVPolynomial, Radix2EvaluationDomain, EvaluationDomain};
    use ark_test_curves::bls12_381::Fr;
    use crate::direct::DirectLDT;

    #[test]
    fn test_direct_ldt(){
        let degree = 51;

        let mut rng = test_rng();
        let poly = DensePolynomial::<Fr>::rand(degree, &mut rng);
        let domain = Radix2EvaluationDomain::new(64).unwrap();
        let evaluations = poly.evaluate_over_domain(domain);

        let low_degree_poly = DirectLDT::generate_low_degree_coefficients(domain.clone(),
                                                                evaluations.evals.to_vec(), degree);

        let sampled_element = domain.element(15);
        let sampled_evaluation = evaluations.evals[15];

        assert!(DirectLDT::verify_low_degree_single_round(sampled_element,
                                                          sampled_evaluation,
                                                          &low_degree_poly))

    }
}