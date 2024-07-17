pub mod direct;
pub mod fri;
pub mod stir;

use ark_ff::FftField;

pub trait Prover<F: FftField> {
    type Proof;
    type ProverConfig;
    type Witness;
    fn new(prover_config: Self::ProverConfig) -> Self;
    fn prove(&self, witness: &Self::Witness) -> Self::Proof;
}
pub trait Verifier<F: FftField> {
    type Statement;
    type Proof;
    type VerifierConfig;
    fn new(verifier_config: Self::VerifierConfig) -> Self;
    fn verify(&self, commitment: &Self::Statement, proof: &Self::Proof) -> bool;
}
pub trait LowDegreeTest<F: FftField> {
    type LDTConfig;
    type Proof;
    type Prover;
    type Verifier;
    fn new(ldt_config: Self::LDTConfig) -> (Self::Prover, Self::Verifier);
}
