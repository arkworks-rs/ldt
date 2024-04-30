// Config
pub struct LDTConfig {
    pub degree: usize,
    pub rate: usize,
}
impl LDTConfig {
    pub fn new(degree: usize, rate: usize) -> Self {
        Self {
            degree,
            rate
        }
    }
}

// LowDegreeTest
pub trait Config<Field> {}
pub trait Prover<Field> {}
pub trait Verifier<Field> {}
pub trait LowDegreeTest<Field> {
    type Config;
    type Prover;
    type Verifier;
    fn config(ldt_config: LDTConfig) -> Self::Config;
    fn prover(config: Self::Config) -> Self::Prover;
    fn verifier(config: Self::Config) -> Self::Verifier;
}
