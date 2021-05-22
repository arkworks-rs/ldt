use ark_ff::PrimeField;
use ark_poly::{Radix2EvaluationDomain, EvaluationDomain, UVPolynomial, Evaluations, Polynomial};
use ark_poly::polynomial::univariate::DensePolynomial;

/// Given domain as `<g>`, `CosetOfDomain` represents `h<g>`
///
/// Constraint equivalent is in `r1cs_std::poly::domain`.
#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub struct Radix2CosetDomain<F: PrimeField> {
    pub base_domain: Radix2EvaluationDomain<F>,
    pub offset: F,
}

impl<F: PrimeField> Radix2CosetDomain<F> {
    /// Returns a new coset domain.
    pub fn new(base_domain: Radix2EvaluationDomain<F>, offset: F) -> Self {
        Radix2CosetDomain {
            base_domain,
            offset,
        }
    }

    /// Returns a coset of size of power of two.
    pub fn new_radix2_coset(coset_size: usize, offset: F) -> Self {
        Self::new(Radix2EvaluationDomain::new(coset_size).unwrap(), offset)
    }

    /// Query a coset according to its position. Returns the positions of coset elements in `self`,
    /// and the result coset represented as domain.
    pub fn query_position_to_coset(&self, coset_position: usize, log_coset_size: usize) -> (Vec<usize>, Self){
        // make sure coset position is not out of range
        assert!(log_coset_size < self.base_domain.log_size_of_group as usize, "query coset size too large");
        assert!(coset_position < (1 << (self.base_domain.log_size_of_group - log_coset_size as u32)), "coset position out of range");

        let dist_between_coset_elems = 1 << (self.base_domain.log_size_of_group as usize - log_coset_size);

        // generate coset
        let mut c = Self::new_radix2_coset(1 << log_coset_size, self.offset * self.gen().pow(&[coset_position as u64]));
        c.base_domain.group_gen = self.gen().pow(&[
            1 << (self.dim() - log_coset_size as usize)]);
        c.base_domain.group_gen_inv = c.base_domain.group_gen.inverse().unwrap();

        // generate positions
        let mut indices = Vec::with_capacity(1 << log_coset_size);
        for i in 0..(1 << (log_coset_size as usize)) {
            indices.push(coset_position + i * dist_between_coset_elems)
        }

        (indices, c)
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
        if self.size() < poly.degree() + 1 {
            // we use naive method. TODO: use a more efficient method
            return self.base_domain.elements().map(|g|poly.evaluate(&(self.offset * g))).collect()
        }
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

#[cfg(test)]
mod tests{
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{UVPolynomial, Polynomial, Radix2EvaluationDomain, EvaluationDomain};
    use ark_std::{test_rng, UniformRand};
    use crate::domain::Radix2CosetDomain;
    use ark_test_curves::bls12_381::Fr;
    use ark_ff::{Field, One};

    #[test]
    fn query_coset_test(){
        let mut rng = test_rng();
        let poly = DensePolynomial::rand(4, &mut rng);

        let offset = Fr::rand(&mut rng);
        let domain_coset = Radix2CosetDomain::new_radix2_coset(15, offset);

        let evals_on_domain_coset = domain_coset.evaluate(&poly);
        assert_eq!(poly.evaluate(&domain_coset.element(2)), evals_on_domain_coset[2]);

        let (query_coset_pos, query_coset)
            = domain_coset.query_position_to_coset(2, 2);

        assert_eq!(query_coset_pos, vec![2,6,10,14]);

        assert_eq!(query_coset.element(0), domain_coset.element(2));
        assert_eq!(query_coset.element(1), domain_coset.element(6));
        assert_eq!(query_coset.element(2), domain_coset.element(10));
        assert_eq!(query_coset.element(3), domain_coset.element(14));

        assert_eq!(query_coset.evaluate(&poly), vec![evals_on_domain_coset[2], evals_on_domain_coset[6], evals_on_domain_coset[10], evals_on_domain_coset[14]])
    }
}