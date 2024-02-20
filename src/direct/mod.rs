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
/// /// This requires communication linear in the degree bound; use FRI for better communication complexity.
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
    /// The polynomial is trimmed to `degree_bound` when necessary.
    pub fn generate_low_degree_coefficients(
        domain: Radix2CosetDomain<F>,
        codewords: Vec<F>,
        degree_bound: usize,
    ) -> DensePolynomial<F> {
        let mut poly = domain.interpolate(codewords);
        // trim higher degree: if poly is higher degree, then the soundness should fail
        poly.coeffs.truncate(degree_bound + 1);
        poly
    }

    /// ### Verifier Side
    ///
    /// The Direct LDT Verify function tests that given a list of coefficients `a_0, a_1, ..., a_{d-1}`
    /// an evaluation point `x`, and claimed evaluation `y`, that `y = \sum_{i =0}^{d} a_i x^i`.
    /// This proves that the provided coefficients of a degree `d` polynomial agree with the claimed
    /// `(evaluation_point, claimed_evaluation)` pair.
    /// This is used to construct a low degree test for an oracle to a claimed polynomials evaluations over a domain.
    /// By sampling enough (domain_element, claimed_evaluation) pairs from the oracle, and testing them
    /// via this method, you become convinced w.h.p. that the oracle is sufficiently close to the claimed coefficients list.
    pub fn verify(
        evaluation_point: F,
        claimed_evaluation: F,
        bounded_coefficients: &DensePolynomial<F>,
    ) -> bool {
        return bounded_coefficients.evaluate(&evaluation_point) == claimed_evaluation;
    }
}

#[cfg(test)]
mod tests {
    use crate::direct::{DirectLDT, Radix2CosetDomain};
    use ark_ff::UniformRand;
    use ark_poly::{univariate::DensePolynomial, UVPolynomial};
    use ark_std::test_rng;
    use ark_test_curves::bls12_381::Fr;

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

        assert!(DirectLDT::verify(
            sampled_element,
            sampled_evaluation,
            &low_degree_poly
        ))
    }
}
