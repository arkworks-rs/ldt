#![forbid(unsafe_code)]

use ark_ff::PrimeField;
use ark_r1cs_std::poly::domain::EvaluationDomain;
use std::marker::PhantomData;
use crate::direct::Radix2CosetDomain;
use crate::domain::Radix2CosetDomain;

pub mod prover;
#[cfg(test)]
mod test;
pub mod verifier;

pub struct FRIParameters<F: PrimeField> {
    /// The degree
    pub tested_degree: u64,
    /// At each round `i`, domain size will shrink to `last_round_domain_size` / `localization_parameters[i]`^2
    pub localization_parameters: Vec<u64>,
    /// Evaluation domain, which is represented as a coset.
    pub domain: Radix2CosetDomain<F>,
}

pub struct FRI<F: PrimeField> {
    _protocol: PhantomData<F>,
}
