#![allow(unused)] // temporary
use crate::domain::Radix2CosetDomain;
use crate::fri::FRIParameters;
use ark_ff::PrimeField;
use ark_r1cs_std::bits::boolean::Boolean;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::poly::domain::Radix2DomainVar;
use ark_r1cs_std::poly::polynomial::univariate::dense::DensePolynomialVar;
use ark_relations::r1cs::SynthesisError;
use ark_sponge::constraints::CryptographicSpongeVar;
use ark_sponge::FieldBasedCryptographicSponge;
use ark_std::marker::PhantomData;

pub struct FRIVerifierGadget<F: PrimeField> {
    _marker: PhantomData<F>,
}

impl<F: PrimeField> FRIVerifierGadget<F> {
    /// Doc: TODO
    pub fn interactive_phase_single_round<
        S: FieldBasedCryptographicSponge<F>,
        SV: CryptographicSpongeVar<F, S>,
    >(
        sponge_var: &mut SV,
    ) -> Result<FpVar<F>, SynthesisError> {
        Ok(sponge_var
            .squeeze_field_elements(1)?
            .first()
            .unwrap()
            .clone())
    }

    /// Doc: TODO
    pub fn sample_coset_index<
        S: FieldBasedCryptographicSponge<F>,
        SV: CryptographicSpongeVar<F, S>,
    >(
        sponge_var: &mut SV,
        fri_parameters: &FRIParameters<F>,
    ) -> Result<Vec<Boolean<F>>, SynthesisError> {
        let log_num_cosets =
            fri_parameters.domain.dim() - fri_parameters.localization_parameters[0] as usize;
        sponge_var.squeeze_bits(log_num_cosets)
    }

    /// Doc: TODO
    pub fn prepare_query(
        rand_coset_index: Vec<Boolean<F>>,
        fri_parameters: &FRIParameters<F>,
    ) -> Result<
        (
            Vec<Radix2DomainVar<F>>,
            Vec<Vec<Boolean<F>>>,
            Radix2CosetDomain<F>,
        ),
        SynthesisError,
    > {
        let num_fri_rounds = fri_parameters.localization_parameters.len();
        let mut coset_indices = Vec::new();
        let mut curr_coset_index = rand_coset_index;
        let mut queries = Vec::with_capacity(num_fri_rounds);
        let mut curr_round_domain = fri_parameters.domain;

        // sample coset index
        for i in 0..num_fri_rounds {
            let log_dist_between_coset_elems =
                curr_round_domain.dim() - fri_parameters.localization_parameters[i] as usize;
            curr_coset_index = curr_coset_index[..log_dist_between_coset_elems].to_vec();

            coset_indices.push(curr_coset_index.clone());

            // get the query coset from coset index
            let query_gen = fri_parameters.domain.gen().pow(&[1
                << (fri_parameters.domain.dim()
                    - fri_parameters.localization_parameters[i] as usize)]);
            debug_assert_eq!(
                query_gen.pow(&[1 << fri_parameters.localization_parameters[i]]),
                F::one()
            );

            let query_offset = &FpVar::constant(curr_round_domain.offset)
                * &(FpVar::constant(curr_round_domain.gen()).pow_le(&curr_coset_index)?);

            let query_coset = Radix2DomainVar {
                gen: query_gen,
                offset: query_offset,
                dim: fri_parameters.localization_parameters[i],
            };

            queries.push(query_coset);

            curr_round_domain = curr_round_domain.fold(fri_parameters.localization_parameters[i])
        }

        Ok((queries, coset_indices, curr_round_domain))
    }

    /// Doc: TODO
    pub fn consistency_check(
        fri_parameters: &FRIParameters<F>,
        queried_coset_indices: &[Vec<Boolean<F>>],
        queries: &[Radix2DomainVar<F>],
        queried_evaluations: &[Vec<FpVar<F>>],
        alphas: &[FpVar<F>],
        final_polynomial_domain: &Radix2CosetDomain<F>,
        final_polynomial: &DensePolynomialVar<F>,
    ) -> Result<Boolean<F>, SynthesisError> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::domain::Radix2CosetDomain;
    use crate::fri::constraints::FRIVerifierGadget;
    use crate::fri::verifier::FRIVerifier;
    use crate::fri::FRIParameters;
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::bits::uint64::UInt64;
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_relations::*;
    use ark_std::{test_rng, UniformRand};
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_prepare_query() {
        let mut rng = test_rng();
        let offset = Fr::rand(&mut rng);
        let domain_input = Radix2CosetDomain::new_radix2_coset(1 << 7, offset);

        let fri_parameters = FRIParameters::new(32, vec![1, 2, 1], domain_input);

        let rand_coset_index = 31usize;
        let cs = ConstraintSystem::new_ref();
        let rand_coset_index_var =
            UInt64::new_witness(ns!(cs, "rand_coset_index"), || Ok(rand_coset_index as u64))
                .unwrap();
        let rand_coset_index_var_arr = rand_coset_index_var.to_bits_le()[..(1 << 6)].to_vec();

        let rand_coset_index = 31;
        let (query_cosets, query_indices, domain_final) =
            FRIVerifier::prepare_query(rand_coset_index, &fri_parameters);
        let (query_cosets_actual, query_indices_actual, domain_final_actual) =
            FRIVerifierGadget::prepare_query(rand_coset_index_var_arr, &fri_parameters).unwrap();

        for i in 0..query_cosets.len() {
            assert_eq!(
                query_cosets_actual[i].offset.value().unwrap(),
                query_cosets[i].offset
            );
            assert_eq!(query_cosets_actual[i].gen, query_cosets[i].gen());
            assert_eq!(query_cosets_actual[i].dim as usize, query_cosets[i].dim());
        }

        assert_eq!(domain_final, domain_final_actual)
    }
}
