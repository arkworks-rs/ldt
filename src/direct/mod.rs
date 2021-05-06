use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, Evaluations, Polynomial, Radix2EvaluationDomain, UVPolynomial};
use std::marker::PhantomData;

pub struct DirectLDT<F: PrimeField> {
    marker: PhantomData<F>,
}

/// Given domain as `<g>`, `CosetOfDomain` represents `h<g>`
///
/// Constraint equivalent is in `r1cs_std::poly::domain`.
#[derive(Clone)]
pub struct CosetDomain<F: PrimeField, D: EvaluationDomain<F>> {
    pub base_domain: D,
    pub offset: F,
}

impl<F: PrimeField> CosetDomain<F, Radix2EvaluationDomain<F>> {
    /// Returns a coset of size of power of two that is large enough for evaluations of a polynomial
    /// having `num_coeffs` coefficients.
    pub fn new_radix2_coset(num_coeffs: usize, offset: F) -> Self {
        Self::new(Radix2EvaluationDomain::new(num_coeffs).unwrap(), offset)
    }
}

impl<F: PrimeField, D: EvaluationDomain<F>> CosetDomain<F, D> {
    /// Returns a new coset domain.
    pub fn new(base_domain: D, offset: F) -> Self {
        CosetDomain {
            base_domain,
            offset,
        }
    }

    /// Given f(x) = \sum a_i x^i. Returns g(x) = \sum a_i h^i x^i
    ///
    /// Note that g(x) = f(hx)
    fn add_offset_to_coeffs(&self, poly: &DensePolynomial<F>) -> DensePolynomial<F> {
        let mut r = F::one();
        let mut transformed_coeff = Vec::with_capacity(poly.coeffs.len());
        for &coeff in poly.coeffs.iter() {
            transformed_coeff.push(coeff * r);
            r *= self.offset
        }
        DensePolynomial::from_coefficients_vec(transformed_coeff)
    }

    /// Given g(x) = \sum a_i h^i x^i. Returns f(x) = \sum a_i x^i
    ///
    /// Note that g(x) = f(hx)
    fn remove_offset_from_coeffs(&self, poly: &DensePolynomial<F>) -> DensePolynomial<F> {
        let mut r = F::one();
        let h_inv = self.offset.inverse().unwrap();
        let mut transformed_coeff = Vec::with_capacity(poly.coeffs.len());
        for &coeff in poly.coeffs.iter() {
            transformed_coeff.push(coeff * r);
            r *= h_inv
        }
        DensePolynomial::from_coefficients_vec(transformed_coeff)
    }

    /// Evaluate polynomial on this coset
    pub fn evaluate(&self, poly: &DensePolynomial<F>) -> Vec<F> {
        // g(x) = f(hx). So, f(coset) = g(base_domain)
        let gx = self.add_offset_to_coeffs(poly);
        gx.evaluate_over_domain(self.base_domain.clone()).evals
    }

    /// given evaluation over this coset. Interpolate and returns coefficients.
    fn interpolate(&self, evaluations: Vec<F>) -> DensePolynomial<F> {
        assert_eq!(evaluations.len(), self.base_domain.size());
        // first get g(x)
        let gx = Evaluations::from_vec_and_domain(evaluations, self.base_domain).interpolate();
        // g(x) = f(hx). Let g(x) = \sum a_i h^i x^i. Then f(x) = \sum a_i x^i
        let fx = self.remove_offset_from_coeffs(&gx);
        fx
    }

    /// Returns an element of the coset
    fn element(&self, i: usize) -> F {
        self.base_domain.element(i) * self.offset
    }
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
    pub fn generate_low_degree_coefficients<D: EvaluationDomain<F>>(
        domain: CosetDomain<F, D>,
        evaluations: Vec<F>,
        degree_bound: usize,
    ) -> DensePolynomial<F> {
        let mut poly = domain.interpolate(evaluations);
        // trim higher degree
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
    use crate::direct::{CosetDomain, DirectLDT};
    use ark_ff::UniformRand;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, Polynomial, Radix2EvaluationDomain, UVPolynomial};
    use ark_std::test_rng;
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_native_coset() {
        let mut rng = test_rng();
        let degree = 51;
        let poly = DensePolynomial::<Fr>::rand(degree, &mut rng);
        let base_domain = Radix2EvaluationDomain::new(64).unwrap();
        let offset = Fr::rand(&mut rng);
        let coset = CosetDomain::new(base_domain, offset);

        // test evaluation
        let expected_eval: Vec<_> = coset
            .base_domain
            .elements()
            .map(|x| poly.evaluate(&(offset * x)))
            .collect();
        let actual_eval = coset.evaluate(&poly);
        assert_eq!(actual_eval, expected_eval);

        // test interpolation
        let interpolated_poly = coset.interpolate(expected_eval);
        assert_eq!(interpolated_poly, poly);
    }

    #[test]
    fn test_direct_ldt() {
        let degree = 51;

        let mut rng = test_rng();
        let poly = DensePolynomial::<Fr>::rand(degree, &mut rng);
        let domain_coset = CosetDomain::new_radix2_coset(52, Fr::rand(&mut rng));
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
