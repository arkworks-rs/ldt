use ark_ff::PrimeField;
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::poly::polynomial::univariate::dense::DensePolynomialVar;
use ark_relations::r1cs::SynthesisError;
use ark_std::marker::PhantomData;

/// Constraints for direct ldt.
pub struct DirectLDTGadget<CF: PrimeField> {
    _marker: PhantomData<CF>,
}

impl<CF: PrimeField> DirectLDTGadget<CF> {
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
        evaluation_point: FpVar<CF>,
        claimed_evaluation: FpVar<CF>,
        coefficients: &DensePolynomialVar<CF>,
        degree_bound: usize,
    ) -> Result<Boolean<CF>, SynthesisError> {
        // make sure the degree is within degree_bound. No need to include degree_bound check
        // in constraints because the verifier can just verify the size of circuit.
        assert!(
            coefficients.coeffs.len() <= degree_bound + 1,
            "polynomial degree out of bound"
        );
        coefficients
            .evaluate(&evaluation_point)?
            .is_eq(&claimed_evaluation)
    }
}
