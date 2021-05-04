use ark_ff::PrimeField;
use std::marker::PhantomData;
use ark_sponge::FieldBasedCryptographicSponge;
use ark_r1cs_std::poly::domain::EvaluationDomain;

/// Implements FRI verifier.
pub struct FRIVerifier<F: PrimeField>  {
    _verifier: PhantomData<F>
}

impl<F: PrimeField> FRIVerifier<F> {
    /// Sample alpha in interactive phase.
    pub fn interactive_phase_single_round<S: FieldBasedCryptographicSponge<F>>(sponge: &mut S) -> F {
        sponge.squeeze_native_field_elements(1)[0]
    }

    /// Prepare all the queries in query phase.
    pub fn prepare_queries<S: FieldBasedCryptographicSponge<F>>(
        sponge: &mut S,
        localization_parameters: Vec<u64>,
        alphas: Vec<F>,
        domain: EvaluationDomain<F>
    ) -> Vec<EvaluationDomain<F>> {
        todo!()
    }
}