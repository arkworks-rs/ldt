use ark_ff::PrimeField;
use ark_r1cs_std::poly::evaluations::univariate::lagrange_interpolator::LagrangeInterpolator;
use std::marker::PhantomData;
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
        let coset_generator = domain.gen().pow(&[1 << (domain.dim() as u64 - localization_param)]);
        let mut cur_coset_offset = domain.offset;

        for coset_index in 0..dist_between_coset_elems {
            let mut poly_evals = Vec::new();
            for intra_coset_index in 0..coset_size {
                poly_evals.push(
                    poly_over_domain[(coset_index + intra_coset_index * dist_between_coset_elems) as usize],
                );
            }

            let interpolator = LagrangeInterpolator::new(
                cur_coset_offset,
                coset_generator,
                localization_param,
                poly_evals,
            );
            new_evals.push(interpolator.interpolate(alpha));
            cur_coset_offset *= domain.gen();
        }

        let c = Radix2CosetDomain::new_radix2_coset(new_evals.len(), domain.offset);
        // c.base_domain.group_gen = coset_generator;
        // c.base_domain.group_gen_inv = coset_generator.inverse().unwrap();
        debug_assert_eq!(coset_generator.pow(&[new_evals.len() as u64]), F::one());
        debug_assert_eq!(c.size(), new_evals.len());
        (c, new_evals)
    }

    /// Returns the domain returned by `interactive_phase_single_round` without supplying evaluations
    /// over domain. This method is useful for verifiers.
    pub fn fold_domain(domain: Radix2CosetDomain<F>, localization_param: u64) -> Radix2CosetDomain<F> {
        let coset_size = 1 << localization_param;
        let domain_size = domain.base_domain.size;
        let dist_between_coset_elems = domain_size / coset_size;
        Radix2CosetDomain::new_radix2_coset(dist_between_coset_elems as usize, domain.offset)
    }
}


#[cfg(test)]
pub mod tests{
    use ark_std::{test_rng, UniformRand};
    use ark_test_curves::bls12_381::Fr;
    use crate::domain::Radix2CosetDomain;
    use crate::direct::DirectLDT;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::UVPolynomial;
    use crate::fri::prover::FRIProver;

    #[test]
    fn degree_reduction_test(){
        let degree = 32;

        let mut rng = test_rng();
        let poly = DensePolynomial::<Fr>::rand(degree, &mut rng);
        let domain_coset = Radix2CosetDomain::new_radix2_coset(64, Fr::rand(&mut rng));
        let evaluations = domain_coset.evaluate(&poly);

        // fri prover should reduce its degree
        let alpha = Fr::rand(&mut rng);
        let localization = 2;
        let (domain_next_round, eval_next_round) =
            FRIProver::interactive_phase_single_round(domain_coset.clone(),
                                                               evaluations.to_vec(),
                                                               localization,
                                                               alpha);

        let low_degree_poly = DirectLDT::generate_low_degree_coefficients(
            domain_next_round.clone(),
            eval_next_round.to_vec(),
            degree / (1 << localization),
        );


        let sampled_element = domain_next_round.element(15);
        let sampled_evaluation = eval_next_round[15];

        assert!(DirectLDT::verify_low_degree_single_round(
            sampled_element,
            sampled_evaluation,
            &low_degree_poly
        ));

        // test `fold_domain`
        let fold_domain = FRIProver::fold_domain(domain_coset.clone(), localization);
        assert_eq!(fold_domain, domain_next_round);
    }
}