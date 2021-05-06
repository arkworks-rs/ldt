use crate::fri::FRIParameters;
use ark_ff::PrimeField;
use ark_r1cs_std::poly::domain::EvaluationDomain;
use ark_sponge::FieldBasedCryptographicSponge;
use std::marker::PhantomData;

/// Implements FRI verifier.
pub struct FRIVerifier<F: PrimeField> {
    _verifier: PhantomData<F>,
}

impl<F: PrimeField> FRIVerifier<F> {
    /// ## Step 1: Interative Phase
    /// Sample alpha in interactive phase.
    pub fn interactive_phase_single_round<S: FieldBasedCryptographicSponge<F>>(
        sponge: &mut S,
    ) -> F {
        sponge.squeeze_native_field_elements(1)[0]
    }

    /// ## Step 2: Query Phase (Prepare)
    /// Prepare all the queries in query phase.
    pub fn prepare_queries<S: FieldBasedCryptographicSponge<F>>(
        sponge: &mut S,
        fri_parameters: &FRIParameters<F>,
    ) -> Vec<EvaluationDomain<F>> {
        let num_rounds = fri_parameters.localization_parameters.len();
        for i in 0..num_rounds {
            let coset_size = 1 << fri_parameters.localization_parameters[i];
            let coset_offset = todo!();

            // coset_generator = root of unity of order coset_size
            // for L's generator g, it equals g^{1 << (L.dim - localization_parameter)}
            let coset_generator = fri_parameters.domain.gen.pow(&[
                1 << (fri_parameters.domain.dim - fri_parameters.localization_parameters[i])
            ]);
        }

        todo!()
    }
}
