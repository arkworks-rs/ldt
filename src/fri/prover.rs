use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree, Path},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, Polynomial};
use ark_std::marker::PhantomData;

use crate::{
    domain::Domain,
    fri::{
        config::FRIConfig,
        proof::{FRIProof, FRIRoundProof},
    },
    ldt::Prover,
    poly_utils,
    utils::{dedup, proof_of_work, squeeze_integer, stack_evaluations},
    commitment::Witness,
};

pub struct FRIProver<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    config: FRIConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}

impl<F: FftField + PrimeField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> Prover<F>
    for FRIProver<F, M, S>
where
    M::InnerDigest: Absorb,
{
    type Config = FRIConfig<M, S>;
    type Proof = FRIProof<F, M>;

    fn new(config: FRIConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn prove(&self, witness: impl Witness<F>) -> Self::Proof {
        // TODO fix this
        // assert!(commitment.polynomials[0].degree() < self.config.starting_degree);

        // get evaluations over a domain
        let domain: Domain<F> =
            Domain::<F>::new(self.config.starting_degree, self.config.starting_rate).unwrap();
        // let evals: Vec<F> = polynomials[0]
        //     .evaluate_over_domain_by_ref(domain.backing_domain)
        //     .evals;
        // let p_evaluations: Vec<Vec<F>> = utils::stack_evaluations(evals, folding_factor);
        let committed_values =
            witness.folded_evaluations_over_domain(domain.clone(), self.config.folding_factor);

        // generate the committment
        let p_commitment = MerkleTree::<M>::new(
            &self.config.merkle_leaf_hash_param,
            &self.config.merkle_two_to_one_param,
            &committed_values,
        )
        .unwrap();

        // assert_eq!(commitment.p_commitment.root(), p_commitment.root());

        let mut sponge: S = S::new(&self.config.sponge_config);
        sponge.absorb(&p_commitment.root());

        let mut g_domain: crate::domain::Domain<F> = domain.clone();
        let mut g_poly: DensePolynomial<F> = witness.polynomial();

        // Commit phase
        let mut commitments: Vec<<M>::InnerDigest> = vec![];
        let mut merkle_trees: Vec<MerkleTree<M>> = vec![p_commitment.clone()];
        let mut folded_evals: Vec<Vec<Vec<F>>> = vec![committed_values.clone()];

        let mut folding_randomness = sponge.squeeze_field_elements(1)[0];
        for _ in 0..self.config.num_rounds {
            // Fold the initial polynomial
            g_poly = poly_utils::folding::poly_fold(
                &g_poly,
                self.config.folding_factor,
                folding_randomness,
            );

            let prev_evals = folded_evals.last().unwrap();

            // The following lines are just precomputations, to avoid having to do inversion
            // and exponentiations in the inner loop
            let domain_size = g_domain.size();
            let generator = g_domain
                .backing_domain
                .element(domain_size / self.config.folding_factor);
            let generator_inv = generator.inverse().unwrap();
            let size_inv = F::from(self.config.folding_factor as u64)
                .inverse()
                .unwrap();
            let coset_offsets: Vec<_> = g_domain
                .backing_domain
                .elements()
                .take(prev_evals.len())
                .collect();
            let mut counter = F::ONE;
            let scale = g_domain.backing_domain.element(1).inverse().unwrap();
            let mut coset_offsets_inv: Vec<_> = vec![];
            for _ in 0..prev_evals.len() {
                coset_offsets_inv.push(counter);
                counter *= scale;
            }

            // Compute the evalations of the folded polynomial
            let g_evaluations: Vec<_> = prev_evals
                .iter()
                .zip(coset_offsets.into_iter())
                .zip(coset_offsets_inv.into_iter())
                .map(|((e, c), ci)| (e, c, ci))
                .map(|(evals, coset_offset, coset_offset_inv)| {
                    poly_utils::interpolation::fft_interpolate(
                        generator,
                        coset_offset,
                        generator_inv,
                        coset_offset_inv,
                        size_inv,
                        evals,
                    )
                    .evaluate(&folding_randomness)
                })
                .collect();

            g_domain = g_domain.scale(self.config.folding_factor);
            //let g_evaluations = g_poly.evaluate_over_domain_by_ref(g_domain.backing_domain).evals;

            let g_folded_evaluations = stack_evaluations(g_evaluations, self.config.folding_factor);
            let g_merkle = MerkleTree::<M>::new(
                &self.config.merkle_leaf_hash_param,
                &self.config.merkle_two_to_one_param,
                &g_folded_evaluations,
            )
            .unwrap();
            let g_root = g_merkle.root();
            sponge.absorb(&g_root);

            folding_randomness = sponge.squeeze_field_elements(1)[0];

            commitments.push(g_root);
            merkle_trees.push(g_merkle);
            folded_evals.push(g_folded_evaluations);
        }

        g_poly =
            poly_utils::folding::poly_fold(&g_poly, self.config.folding_factor, folding_randomness);

        // Query phase
        let mut folded_evals_len = domain.size() / self.config.folding_factor;
        let mut query_indexes = dedup(
            (0..self.config.repetitions).map(|_| squeeze_integer(&mut sponge, folded_evals_len)),
        );

        // Note that we include final round as well
        let mut round_proofs = vec![];
        for round in 0..=self.config.num_rounds {
            let queries_to_prev_ans = query_indexes
                .iter()
                .map(|&index| folded_evals[round][index].clone())
                .collect();
            // let queries_to_prev_proof = merkle_trees[round]
            //     .generate_multi_proof(query_indexes.clone())
            //     .unwrap();
            // get the openings
            let mut queries_to_prev_proof: Vec<Path<M>> = Vec::with_capacity(query_indexes.len());
            for query in query_indexes.clone() {
                queries_to_prev_proof.push(merkle_trees[round].generate_proof(query).unwrap());
            }
            let queries_to_prev: (Vec<Vec<F>>, Vec<Path<M>>) =
                (queries_to_prev_ans, queries_to_prev_proof);

            folded_evals_len = folded_evals_len / self.config.folding_factor;
            query_indexes = dedup(query_indexes.into_iter().map(|i| i % folded_evals_len));

            round_proofs.push(FRIRoundProof { queries_to_prev });
        }

        Self::Proof {
            polynomial: g_poly,
            commitments,
            round_proofs,
            proof_of_work_nonce: proof_of_work(&mut sponge, self.config.proof_of_work_bits),
            initial_p_commitment_root: p_commitment.root(),
        }
    }
}
