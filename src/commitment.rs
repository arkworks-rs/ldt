use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, MerkleTree, TwoToOneParam},
    sponge::Absorb,
};
use ark_ff::FftField;
use ark_poly::univariate::DensePolynomial;

use crate::{domain::Domain, utils};

// BTW I don't see how this works without adding MerkleConfig to LDT
// pub trait Commitment<F: FftField, M: MerkleConfig> {
//     fn new(
//         degree: usize,
//         merkle_leaf_hash_param: LeafParam<M>,
//         merkle_two_to_one_param: TwoToOneParam<M>,
//         polynomials: &[DensePolynomial<F>],
//     ) -> Self;
//     fn get_p_commitment(&self) -> MerkleTree<M>;
//     fn get_p_evaluations_over_domain(&self) -> Vec<Vec<F>>;
// }

pub struct Commitment<F: FftField, M: MerkleConfig> {
    pub domain: Domain<F>,
    pub polynomials: Vec<DensePolynomial<F>>,
    pub p_evaluations: Vec<Vec<F>>,
    pub p_commitment: MerkleTree<M>,
}

impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>> Commitment<F, M>
where
    M::InnerDigest: Absorb,
{
    pub fn new(
        starting_degree: usize,
        starting_rate: usize,
        folding_factor: usize,
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        polynomials: Vec<DensePolynomial<F>>, // TODO: What if we just keep this as Vec?
    ) -> Self {
        // get evaluations over a domain
        let domain: Domain<F> = Domain::<F>::new(starting_degree, starting_rate).unwrap();
        let evals: Vec<F> = polynomials[0]
            .evaluate_over_domain_by_ref(domain.backing_domain)
            .evals;
        let p_evaluations: Vec<Vec<F>> = utils::stack_evaluations(evals, folding_factor);

        // generate the committment
        let p_commitment = MerkleTree::<M>::new(
            &merkle_leaf_hash_param,
            &merkle_two_to_one_param,
            &p_evaluations,
        )
        .unwrap();

        Self {
            domain,
            polynomials,
            p_commitment,
            p_evaluations,
        }
    }
}
