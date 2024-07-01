use ark_ff::FftField;

use crate::witness::Witness;

pub trait Prover<F: FftField> {
    type Config;
    type Proof;
    fn new(config: Self::Config) -> Self;
    fn prove(&self, argument: impl Witness<F>) -> Self::Proof;
}
pub trait Verifier<F: FftField> {
    type Config;
    type Proof;
    fn new(config: Self::Config) -> Self;
    fn verify(&self, proof: &Self::Proof) -> bool;
}
pub trait LowDegreeTest<F: FftField> {
    type Config;
    type Proof;
    type Prover;
    type Verifier;
    fn new(config: Self::Config) -> (Self::Prover, Self::Verifier);
}
