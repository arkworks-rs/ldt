use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::EvaluationDomain;

use crate::domain::Domain;

use super::config::STIRConfig;

pub struct STIRVerifierState<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    pub config: STIRConfig<M, S>,
    pub domain_gen: F,
    pub domain_offset: F,
    pub domain_size: usize,
    pub folding_randomness: F,
    pub root_of_unity: F,
    pub round_num: usize,
    pub sponge: S,
}

impl<F, M, S> STIRVerifierState<F, M, S>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
{
    pub fn new(config: STIRConfig<M, S>, commitment_digest: M::InnerDigest) -> Self {
        let mut sponge = S::new(&config.sponge_config);
        sponge.absorb(&commitment_digest);
        let folding_randomness = sponge.squeeze_field_elements(1)[0];

        let domain = Domain::<F>::new(config.starting_degree, config.starting_rate).unwrap();

        let domain_gen = domain.element(1);
        let domain_size = domain.size();
        Self {
            config,
            domain_gen,
            domain_offset: F::one(),
            domain_size,
            folding_randomness,
            root_of_unity: domain_gen,
            round_num: 0,
            sponge,
        }
    }
}
