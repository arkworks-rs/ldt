use ark_ff::FftField;
use ark_poly::univariate::DensePolynomial;

use crate::{domain::Domain, utils::stack_evaluations};

pub trait Witness<F: FftField> {
    fn folded_evaluations_over_domain(
        &self,
        domain: Domain<F>,
        folding_factor: usize,
    ) -> Vec<Vec<F>>;
}

pub struct SingleWitness<F: FftField> {
    pub polynomial: DensePolynomial<F>,
}

impl<F: FftField> Witness<F> for SingleWitness<F> {
    fn folded_evaluations_over_domain(
        &self,
        domain: Domain<F>,
        folding_factor: usize,
    ) -> Vec<Vec<F>> {
        let evals: Vec<F> = self
            .polynomial
            .evaluate_over_domain_by_ref(domain.backing_domain)
            .evals;
        stack_evaluations(evals, folding_factor)
    }
}
