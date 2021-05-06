#![forbid(unsafe_code)]

use ark_ff::PrimeField;
use ark_r1cs_std::poly::domain::EvaluationDomain;
use std::marker::PhantomData;

pub mod prover;
#[cfg(test)]
mod test;
pub mod verifier;

pub struct FRIParameters<F: PrimeField> {
    /// The degree
    pub tested_degree: u64,
    /// At each round `i`, domain size will shrink to `last_round_domain_size` / `localization_parameters[i]`^2
    pub localization_parameters: Vec<u64>,
    /// Evaluation domain, which is represented as a coset. Note that this domain is different
    /// from evaluation domain defined in `ark-poly`, which is strictly a subgroup with no offset.
    pub domain: EvaluationDomain<F>,
}

pub struct FRI<F: PrimeField> {
    _protocol: PhantomData<F>,
}
