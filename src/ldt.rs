use ark_ff::FftField;
use ark_poly::univariate::DensePolynomial;

pub trait Prover<F: FftField> {
    type Config;
    type Commitment;
    type Proof;
    fn new(config: Self::Config) -> Self;
    fn commit(&self, polynomial: DensePolynomial<F>) -> Self::Commitment;
    fn prove(&self, commitment: Self::Commitment) -> Self::Proof;
}
pub trait Verifier<F: FftField> {
    
}
pub trait LowDegreeTest<F: FftField> {
    type Proof;
    type Config;
    type Prover;
    type Verifier;
    fn prover(config: Self::Config) -> Self::Prover;
    fn verifier(config: Self::Config) -> Self::Verifier;
}
