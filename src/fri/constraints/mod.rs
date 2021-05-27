#![warn(unused)]
use crate::fri::FRIParameters;
use ark_ff::PrimeField;
use ark_r1cs_std::bits::boolean::Boolean;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::poly::domain::Radix2Domain;
use ark_r1cs_std::poly::polynomial::univariate::dense::DensePolynomialVar;
use ark_relations::r1cs::SynthesisError;
use ark_sponge::constraints::CryptographicSpongeVar;
use ark_sponge::FieldBasedCryptographicSponge;
use ark_std::marker::PhantomData;

pub struct FRIVerifierGadget<F: PrimeField> {
    _marker: PhantomData<F>,
}

impl<F: PrimeField> FRIVerifierGadget<F> {
    pub fn interactive_phase_single_round<
        S: FieldBasedCryptographicSponge<F>,
        SV: CryptographicSpongeVar<F, S>,
    >(
        sponge_var: &mut SV,
    ) -> Result<FpVar<F>, SynthesisError> {
        // currently r1cs-std is not in cargo, so sponge has some compatability issue
        todo!()
    }

    pub fn sample_coset_index<
        S: FieldBasedCryptographicSponge<F>,
        SV: CryptographicSpongeVar<F, S>,
    >(
        sponge_var: &mut SV,
        fri_parameters: &FRIParameters<F>,
    ) -> Result<Vec<Boolean<F>>, SynthesisError> {
        todo!()
    }

    pub fn prepare_queries(
        rand_coset_index: Vec<Boolean<F>>,
        fri_parameters: &FRIParameters<F>,
    ) -> Result<(Vec<Radix2Domain<F>>, Vec<Vec<Boolean<F>>>, Radix2Domain<F>), SynthesisError> {
        todo!()
    }

    pub fn consistency_check(
        fri_parameters: &FRIParameters<F>,
        queried_coset_indices: &[Vec<Boolean<F>>],
        queries: &[Radix2Domain<F>],
        queried_evaluations: &[Vec<FpVar<F>>],
        alphas: &[FpVar<F>],
        final_polynomial_domain: &Radix2Domain<F>,
        final_polynomial: &DensePolynomialVar<F>,
    ) -> Result<Boolean<F>, SynthesisError> {
        todo!()
    }
}
