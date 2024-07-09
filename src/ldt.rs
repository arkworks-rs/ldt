use ark_ff::FftField;

pub trait Prover<F: FftField> {
    type Proof;
    type ProverConfig;
    type Witness;
    fn new(config: Self::ProverConfig) -> Self;
    fn prove(&self, witness: &Self::Witness) -> Self::Proof;
}
pub trait Verifier<F: FftField> {
    type Proof;
    type VerifierConfig;
    fn new(config: Self::VerifierConfig) -> Self;
    fn verify(&self, proof: &Self::Proof) -> bool;
}
pub trait LowDegreeTest<F: FftField> {
    type LDTConfig;
    type Proof;
    type Prover;
    type Verifier;
    fn new(config: Self::LDTConfig) -> (Self::Prover, Self::Verifier);
}
