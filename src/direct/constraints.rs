use ark_ff::PrimeField;
use ark_r1cs_std::boolean::Boolean;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::poly::polynomial::univariate::dense::DensePolynomialVar;
use ark_relations::r1cs::SynthesisError;
use ark_std::marker::PhantomData;

pub struct DirectLDTGadget<CF: PrimeField> {
    _marker: PhantomData<F>,
}

impl<CF: PrimeField> DirectLDTGadget<CF> {
    /// ### Verifier Side
    ///
    /// Verifier sample one element from domain and get its evaluation. Check if that evaluation
    /// agrees with low-degree polynomial. Assume that `DensePolynomial` is low-degree, and verifier
    /// need to check that.
    #[tracing::instrument(target = "r1cs", skip(self, other))]
    pub fn verify_low_degree_single_round(
        sampled_domain_element: FpVar<CF>,
        sampled_evaluation_element: FpVar<CF>,
        coefficients: &DensePolynomialVar<CF>,
    ) -> Result<Boolean<CF>, SynthesisError> {
        coefficients
            .evaluate(&sampled_domain_element)?
            .is_eq(&sampled_evaluation_element)
    }
}
