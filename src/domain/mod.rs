use ark_ff::PrimeField;
use ark_poly::{Radix2EvaluationDomain, EvaluationDomain, UVPolynomial, Evaluations};
use ark_poly::polynomial::univariate::DensePolynomial;

/// Given domain as `<g>`, `CosetOfDomain` represents `h<g>`
///
/// Constraint equivalent is in `r1cs_std::poly::domain`.
#[derive(Clone)]
pub struct Radix2CosetDomain<F: PrimeField> {
    pub base_domain: Radix2EvaluationDomain<F>,
    pub offset: F,
}

impl<F: PrimeField> Radix2CosetDomain<F> {
    /// Returns a new coset domain.
    pub fn new(base_domain: D, offset: F) -> Self {
        Radix2CosetDomain {
            base_domain,
            offset,
        }
    }

    /// Returns a coset of size of power of two that is large enough for evaluations of a polynomial
    /// having `num_coeffs` coefficients.
    pub fn new_radix2_coset(num_coeffs: usize, offset: F) -> Self {
        Self::new(Radix2EvaluationDomain::new(num_coeffs).unwrap(), offset)
    }



    /// returns the size of the domain
    pub fn size(&self) -> usize {
        self.base_domain.size()
    }

    /// return the log 2 size of domain
    pub fn dim(&self) -> usize {
        self.base_domain.log_size_of_group as usize
    }

    /// returns generator of the coset
    pub fn gen(&self) -> F {
        self.base_domain.group_gen
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
    pub fn interpolate(&self, evaluations: Vec<F>) -> DensePolynomial<F> {
        assert_eq!(evaluations.len(), self.base_domain.size());
        // first get g(x)
        let gx = Evaluations::from_vec_and_domain(evaluations, self.base_domain).interpolate();
        // g(x) = f(hx). Let g(x) = \sum a_i h^i x^i. Then f(x) = \sum a_i x^i
        let fx = self.remove_offset_from_coeffs(&gx);
        fx
    }

    /// Returns an element of the coset
    pub fn element(&self, i: usize) -> F {
        self.base_domain.element(i) * self.offset
    }
}