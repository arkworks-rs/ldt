use crate::domain::Radix2CosetDomain;
use ark_ff::PrimeField;
use ark_std::marker::PhantomData;
use ark_std::vec::Vec;
/// R1CS constraints for FRI Verifier.
#[cfg(feature = "r1cs")]
pub mod constraints;
/// Prover used by FRI protocol.
pub mod prover;
/// Verifier used by FRI protocol.
pub mod verifier;

/// Some parameters used by FRI verifiers.
#[derive(Clone)]
pub struct FRIParameters<F: PrimeField> {
    /// The degree
    pub tested_degree: u64,
    /// At each round `i`, domain size will shrink to `last_round_domain_size` / `localization_parameters[i]`^2
    pub localization_parameters: Vec<u64>,
    /// Evaluation domain, which is represented as a coset.
    pub domain: Radix2CosetDomain<F>,
    /// coset sizes in each round (first round is input coset)
    log_round_coset_sizes: Vec<usize>,
}

impl<F: PrimeField> FRIParameters<F> {
    /// Check parameter validity and returns new `FRIParameters`.
    pub fn new(
        tested_degree: u64,
        localization_parameters: Vec<u64>,
        domain: Radix2CosetDomain<F>,
    ) -> Self {
        assert!(
            domain.size() >= tested_degree as usize + 1,
            "Evaluations is not low degree!\
                Domain size needs to be >= tested_degree + 1"
        );
        let mut log_round_coset_sizes = Vec::new();
        log_round_coset_sizes.push(domain.dim());
        for i in 0..localization_parameters.len() {
            log_round_coset_sizes
                .push(log_round_coset_sizes[i] - localization_parameters[i] as usize)
        }
        FRIParameters {
            tested_degree,
            localization_parameters,
            domain,
            log_round_coset_sizes,
        }
    }
}

/// Fast Reed-Solomon Interactive Oracle Proof of Proximity
pub struct FRI<F: PrimeField> {
    _protocol: PhantomData<F>,
}
