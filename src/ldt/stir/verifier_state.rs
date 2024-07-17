use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain};

use crate::{domain::Domain, poly_utils};

use super::config::STIRConfig;

pub struct STIRVerifierState<F, M, S>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
{
    pub comb_randomness: F,
    pub config: STIRConfig<M, S>,
    pub domain_gen: F,
    pub domain_offset: F,
    pub domain_size: usize,
    pub folding_randomness: F,
    pub interpolating_polynomial: DensePolynomial<F>,
    pub quotient_set: Vec<F>,
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
            comb_randomness: F::zero(),
            config,
            domain_gen,
            domain_offset: F::one(),
            domain_size,
            folding_randomness,
            interpolating_polynomial: DensePolynomial::from_coefficients_vec(vec![]),
            quotient_set: vec![],
            root_of_unity: domain_gen,
            round_num: 0,
            sponge,
        }
    }
    pub fn sponge_absorb(&mut self, element: impl Absorb) {
        self.sponge.absorb(&element);
    }
    pub fn sponge_squeeze(&mut self) -> F {
        self.sponge.squeeze_field_elements(1)[0]
    }
    pub fn sponge_squeeze_multiple(&mut self, num_elements: usize) -> Vec<F> {
        self.sponge.squeeze_field_elements(num_elements)
    }
    // TODO: Nuke this
    pub fn query(
        &self,
        evaluation_point: F,
        value_of_prev_oracle: F,
        common_factors_inverse: F,
        denom_hint: F,
        ans_eval: F,
    ) -> F {
        match &self.round_num {
            0 => value_of_prev_oracle, // In case this is the initial function, we just return the value of the previous oracle
            _ => {
                let num_terms = self.quotient_set.len();
                let quotient_evaluation = poly_utils::quotient::quotient_with_hint(
                    value_of_prev_oracle,
                    evaluation_point,
                    &self.quotient_set,
                    denom_hint,
                    ans_eval,
                );

                let common_factor = evaluation_point * self.comb_randomness;

                let scale_factor = if common_factor != F::ONE {
                    (F::ONE - common_factor.pow([(num_terms + 1) as u64])) * common_factors_inverse
                } else {
                    F::from((num_terms + 1) as u64)
                };

                quotient_evaluation * scale_factor
            }
        }
    }
}
