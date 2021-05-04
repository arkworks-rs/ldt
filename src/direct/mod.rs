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
        // trim to degree bound
        poly.coeffs.truncate(degree_bound);
        poly
    }

    /// ### Verifier Side
    ///
    /// Verifier sample one element from domain and get its evaluation. Check if that evaluation
    /// agrees with low-degree polynomial.
    pub fn verify_low_degree_single_round<D: EvaluationDomain<F>>(sampled_domain_element: F,
                                                     sampled_evaluation_element: F,
                                                     low_degree_interpolated_polynomial: &DensePolynomial<F>) -> bool {
        return low_degree_interpolated_polynomial.evaluate(&sampled_domain_element) == sampled_evaluation_element;
    }

}

