use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{
    UVPolynomial, EvaluationDomain, Evaluations, Polynomial, Radix2EvaluationDomain,
};
#[cfg(feature = "r1cs")]
use ark_r1cs_std::bits::boolean::Boolean;
#[cfg(feature = "r1cs")]
use ark_r1cs_std::fields::fp::FpVar;
#[cfg(feature = "r1cs")]
use ark_r1cs_std::fields::FieldVar;
#[cfg(feature = "r1cs")]
use ark_relations::r1cs::SynthesisError;
use ark_std::vec::Vec;

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

// TODO: Move this to algebra, per https://github.com/arkworks-rs/algebra/issues/88#issuecomment-734963835
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

    /// Converts a query position to the elements of the unique coset of size `log_coset_size`
    /// within this domain that the query lies in.
    /// `query_position` is an index within this domain.
    /// Returns the positions of coset elements in `self`,
    /// and the coset represented as a Radix2CosetDomain.
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
            // we use naive method for evaluating a polynomial larger than the domain size.
            // TODO: use a more efficient method using the fact that:
            // (hg)^{|base_domain|} = h^{|base_domain|},
            // so we can efficiently fold the polynomial's coefficients on itself,
            // into a single polynomial of degree `self.size() - 1`
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
    /// Returns an element fo the coset, given the index as a variable.
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
    use ark_poly::{UVPolynomial, Polynomial};
    use ark_std::{test_rng, UniformRand};
    use ark_test_curves::bls12_381::Fr;

    use crate::domain::Radix2CosetDomain;

    #[cfg(feature = "r1cs")]
    mod consistency_with_constraints {
        use ark_poly::univariate::DensePolynomial;
        use ark_poly::Radix2EvaluationDomain;
        use ark_poly::{DenseUVPolynomial, EvaluationDomain, Polynomial};
        use ark_r1cs_std::alloc::AllocVar;
        use ark_r1cs_std::fields::fp::FpVar;
        use ark_r1cs_std::fields::FieldVar;
        use ark_r1cs_std::poly::domain::Radix2DomainVar;
        use ark_r1cs_std::poly::evaluations::univariate::EvaluationsVar;
        use ark_r1cs_std::R1CSVar;
        use ark_relations::r1cs::ConstraintSystem;
        use ark_std::{test_rng, UniformRand};
        use ark_test_curves::bls12_381::Fr;

        use crate::domain::Radix2CosetDomain;

        #[test]
        fn test_consistency_with_coset_constraints() {
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

            let r1cs_coset = Radix2DomainVar::new(
                base_domain.group_gen,
                ark_std::log2(degree.next_power_of_two()) as u64,
                FpVar::constant(offset),
            )
            .unwrap();
            let eval_var = EvaluationsVar::from_vec_and_domain(eval_var, r1cs_coset, true);

            let pt = Fr::rand(&mut rng);
            let pt_var =
                FpVar::new_witness(ark_relations::ns!(cs, "random point"), || Ok(pt)).unwrap();

            let expected = poly.evaluate(&pt);
            let actual = eval_var.interpolate_and_evaluate(&pt_var).unwrap();

            assert_eq!(actual.value().unwrap(), expected);
            assert!(cs.is_satisfied().unwrap());
        }
    }

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
