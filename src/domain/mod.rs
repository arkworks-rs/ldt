use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, Evaluations, Polynomial, Radix2EvaluationDomain, UVPolynomial};
use ark_std::vec::Vec;

#[cfg(feature = "r1cs")]
use ark_r1cs_std::bits::boolean::Boolean;
#[cfg(feature = "r1cs")]
use ark_r1cs_std::fields::fp::FpVar;
#[cfg(feature = "r1cs")]
use ark_r1cs_std::fields::FieldVar;
#[cfg(feature = "r1cs")]
use ark_relations::r1cs::SynthesisError;

/// Given domain as `<g>`, `CosetOfDomain` represents `h<g>`
///
/// Constraint equivalent is in `r1cs_std::poly::domain`.
#[derive(Clone, Copy, Eq, PartialEq, Debug)]
pub struct Radix2CosetDomain<F: PrimeField> {
    /// A non-coset radix 2 domain: `<g>`
    pub base_domain: Radix2EvaluationDomain<F>,
    /// offset `h`
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
    pub fn query_position_to_coset(
        &self,
        query_position: usize,
        log_coset_size: usize,
    ) -> (Vec<usize>, Self) {
        // make sure coset position is not out of range
        assert!(
            log_coset_size < self.base_domain.log_size_of_group as usize,
            "query coset size too large"
        );
        assert!(
            query_position < (1 << (self.base_domain.log_size_of_group - log_coset_size as u32)),
            "coset position out of range"
        );

        let dist_between_coset_elems =
            1 << (self.base_domain.log_size_of_group as usize - log_coset_size);

        // generate coset
        let c = Self::new_radix2_coset(
            1 << log_coset_size,
            self.offset * self.gen().pow(&[query_position as u64]),
        );
        // c.base_domain.group_gen = self.gen().pow(&[1 << (self.dim() - log_coset_size)]);
        // c.base_domain.group_gen_inv = c.base_domain.group_gen.inverse().unwrap(); // not necessary

        // generate positions
        let mut indices = Vec::with_capacity(1 << log_coset_size);
        for i in 0..(1 << log_coset_size) {
            indices.push(query_position + i * dist_between_coset_elems)
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
            return self
                .base_domain
                .elements()
                .map(|g| poly.evaluate(&(self.offset * g)))
                .collect();
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

    #[cfg(feature = "r1cs")]
    pub fn element_var(&self, index: &[Boolean<F>]) -> Result<FpVar<F>, SynthesisError> {
        Ok(FpVar::constant(self.offset) * FpVar::constant(self.gen()).pow_le(index)?)
    }

    /// Shrink the domain size such that new domain size = `self.size() / (1 << log_shrink_factor)`
    /// and has same offset.
    pub fn fold(&self, log_shrink_factor: u64) -> Radix2CosetDomain<F> {
        let coset_size = 1 << log_shrink_factor;
        let domain_size = self.base_domain.size;
        let dist_between_coset_elems = domain_size / coset_size;
        Radix2CosetDomain::new_radix2_coset(dist_between_coset_elems as usize, self.offset)
    }
}

#[cfg(test)]
mod tests {
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{Polynomial, UVPolynomial};
    use ark_std::{test_rng, UniformRand};
    use ark_test_curves::bls12_381::Fr;

    use crate::domain::Radix2CosetDomain;

    #[test]
    fn query_coset_test() {
        let mut rng = test_rng();
        let poly = DensePolynomial::rand(4, &mut rng);

        let offset = Fr::rand(&mut rng);
        let domain_coset = Radix2CosetDomain::new_radix2_coset(15, offset);

        let evals_on_domain_coset = domain_coset.evaluate(&poly);
        assert_eq!(
            poly.evaluate(&domain_coset.element(2)),
            evals_on_domain_coset[2]
        );

        let (query_coset_pos, query_coset) = domain_coset.query_position_to_coset(2, 2);

        assert_eq!(query_coset_pos, vec![2, 6, 10, 14]);

        assert_eq!(query_coset.element(0), domain_coset.element(2));
        assert_eq!(query_coset.element(1), domain_coset.element(6));
        assert_eq!(query_coset.element(2), domain_coset.element(10));
        assert_eq!(query_coset.element(3), domain_coset.element(14));

        assert_eq!(
            query_coset.evaluate(&poly),
            vec![
                evals_on_domain_coset[2],
                evals_on_domain_coset[6],
                evals_on_domain_coset[10],
                evals_on_domain_coset[14]
            ]
        )
    }

    #[test]
    #[cfg(feature = "r1cs")]
    fn element_var_test() {
        use ark_r1cs_std::alloc::AllocVar;
        use ark_r1cs_std::uint64::UInt64;
        use ark_r1cs_std::R1CSVar;
        use ark_relations::r1cs::ConstraintSystem;
        use ark_relations::*;

        let mut rng = test_rng();
        let offset = Fr::rand(&mut rng);
        let domain_coset = Radix2CosetDomain::new_radix2_coset(15, offset);

        let cs = ConstraintSystem::new_ref();
        let index = 11;
        let index_var = UInt64::new_witness(ns!(cs, "index"), || Ok(index))
            .unwrap()
            .to_bits_le();

        let expected = domain_coset.element(index as usize);
        let actual = domain_coset
            .element_var(&index_var)
            .unwrap()
            .value()
            .unwrap();

        assert_eq!(expected, actual);
        assert!(cs.is_satisfied().unwrap())
    }
}
