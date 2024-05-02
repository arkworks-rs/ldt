use ark_ff::FftField;
use ark_std::marker::PhantomData;

// Config
pub struct LDTConfig<F: FftField> {
    pub degree: usize,
    pub rate: usize,
    _field: PhantomData<F>,
}
impl<F: FftField> LDTConfig<F> {
    pub fn new(degree: usize, rate: usize) -> Self {
        Self {
            degree,
            rate,
            _field: PhantomData::<F>,
        }
    }
}

// LowDegreeTest
pub trait Config<F: FftField> {}
pub trait Prover<F: FftField> {
    fn commit(&self, evals: Vec<Vec<F>>);
    fn prove(&self, evals: Vec<Vec<F>>);
}
pub trait Verifier<F: FftField> {}
pub trait LowDegreeTest<F: FftField> {
    type Proof;
    type Config;
    type Prover;
    type Verifier;
    fn config(ldt_config: LDTConfig<F>) -> Self::Config;
    fn prover(config: Self::Config) -> Self::Prover;
    fn verifier(config: Self::Config) -> Self::Verifier;
}
