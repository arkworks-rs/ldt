use ark_ff::PrimeField;
use ark_r1cs_std::poly::domain::EvaluationDomain;
use ark_r1cs_std::poly::evaluations::univariate::lagrange_interpolator::LagrangeInterpolator;
use std::marker::PhantomData;
use crate::direct::Radix2CosetDomain;
use crate::domain::Radix2CosetDomain;

pub struct FRIProver<F: PrimeField> {
    _prover: PhantomData<F>,
}

impl<F: PrimeField> FRIProver<F> {
    /// Single round prover in commit phase. Returns the polynomial for next round
    /// represented by evaluations over domain in next round.
    ///
    /// The prover is inefficient. TODO: Adapt code from libiop.
    ///
    /// Returns domain for next round polynomial and evaluations over the domain.
    pub fn interactive_phase_single_round(
        domain: Radix2CosetDomain<F>,
        poly_over_domain: Vec<F>,
        localization_param: u64,
        alpha: F,
    ) -> (Radix2CosetDomain<F>, Vec<F>) {
        let coset_size = 1 << localization_param;
        let domain_size = domain.base_domain.size;
        let dist_between_coset_elems = domain_size / coset_size;
        let mut new_evals = Vec::with_capacity(dist_between_coset_elems as usize);
        let coset_generator = domain.gen.pow(&[1 << (domain.dim - localization_param)]);
        let mut cur_coset_offset = domain.offset;

        for coset_index in 0..dist_between_coset_elems {
            let mut poly_evals = Vec::new();
            for intra_coset_index in 0..coset_size {
                poly_evals.push(
                    poly_over_domain[coset_index + intra_coset_index * dist_between_coset_elems],
                );
            }

            let interpolator = LagrangeInterpolator::new(
                cur_coset_offset,
                coset_generator,
                localization_param,
                poly_evals,
            );
            new_evals.push(interpolator.interpolate(alpha));
            cur_coset_offset *= domain.gen;
        }

        let mut c = Radix2CosetDomain::new(coset_size, domain.offset);
        c.base_domain.group_gen = coset_generator;

        (c, new_evals)
    }
}
