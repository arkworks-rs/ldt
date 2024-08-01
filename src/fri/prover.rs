use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree, Path},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, Polynomial};
use ark_std::marker::PhantomData;

#[cfg(not(feature = "std"))]
use ark_std::{vec, vec::Vec};

use crate::{
    fri::{
        config::FRIConfig,
        proof::{FRIProof, FRIProofRound},
    },
    ldt::Prover,
    poly_utils,
    utils::{dedup, proof_of_work, squeeze_integer, stack_evaluations},
    witness::Witness,
};

pub struct FRIProver<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
    W: Witness<F, M>,
{
    prover_config: FRIConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge: PhantomData<S>,
    _witness: PhantomData<W>,
}

impl<F, M, S, W> Prover<F> for FRIProver<F, M, S, W>
where
    F: FftField + PrimeField,
    M: MerkleConfig<Leaf = Vec<F>> + Clone,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M, Commitment = MerkleTree<M>, CommittedValues = Vec<Vec<F>>>
        + Clone,
    W::ChallengeAnswers: Clone,
{
    type Witness = W;
    type ProverConfig = FRIConfig<W::MerkleConfig, S>;
    type Proof = FRIProof<F, W::MerkleConfig, S>;

    fn new(prover_config: FRIConfig<W::MerkleConfig, S>) -> Self {
        Self {
            prover_config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<W::MerkleConfig>,
            _sponge: PhantomData::<S>,
            _witness: PhantomData::<W>,
        }
    }
    fn prove(&self, witness: &W) -> Self::Proof {
        // assert!(witness.coeff_degree() < self.config.starting_degree);

        // Initialize a sponge with the commitment digest
        let mut sponge: S = S::new(&self.prover_config.sponge_config);
        sponge.absorb(&witness.commitment_digest());

        let mut g_domain: crate::domain::Domain<F> = witness.domain();
        let mut g_poly: DensePolynomial<F> = witness.coeff();

        // Commit phase
        let mut commitments: Vec<M::InnerDigest> = vec![];
        let mut merkle_trees: Vec<MerkleTree<M>> = vec![witness.commitment()];
        let mut folded_evals: Vec<Vec<Vec<F>>> = vec![witness.committed_values()];

        let mut folding_randomness = sponge.squeeze_field_elements(1)[0];
        for _ in 0..self.prover_config.num_rounds {
            // Fold the initial polynomial
            g_poly = poly_utils::folding::poly_fold(
                &g_poly,
                self.prover_config.folding_factor,
                folding_randomness,
            );

            let prev_evals = folded_evals.last().unwrap();

            // The following lines are just precomputations, to avoid having to do inversion
            // and exponentiations in the inner loop
            let domain_size = g_domain.size();
            let generator = g_domain
                .backing_domain
                .element(domain_size / self.prover_config.folding_factor);
            let generator_inv = generator.inverse().unwrap();
            let size_inv = F::from(self.prover_config.folding_factor as u64)
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

            g_domain = g_domain.scale(self.prover_config.folding_factor);
            //let g_evaluations = g_poly.evaluate_over_domain_by_ref(g_domain.backing_domain).evals;

            let g_folded_evaluations =
                stack_evaluations(g_evaluations, self.prover_config.folding_factor);
            let g_merkle = MerkleTree::<W::MerkleConfig>::new(
                &self.prover_config.merkle_leaf_hash_param,
                &self.prover_config.merkle_two_to_one_param,
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

        g_poly = poly_utils::folding::poly_fold(
            &g_poly,
            self.prover_config.folding_factor,
            folding_randomness,
        );

        // Query phase
        let mut folded_evals_len = witness.domain().size() / self.prover_config.folding_factor;
        let mut query_indexes = dedup(
            (0..self.prover_config.repetitions)
                .map(|_| squeeze_integer(&mut sponge, folded_evals_len)),
        );

        // Note that we include final round as well
        let mut round_proofs = vec![];
        for round in 0..=self.prover_config.num_rounds {
            let queries_to_prev_ans: Vec<Vec<F>> = query_indexes
                .iter()
                .map(|&index| folded_evals[round][index].clone())
                .collect();
            let challenge_answers = merkle_trees[round]
                .generate_multi_proof(query_indexes.clone())
                .unwrap();
            // get the openings
            let mut queries_to_prev_proof: Vec<Path<W::MerkleConfig>> =
                Vec::with_capacity(query_indexes.len());
            for query in query_indexes.clone() {
                queries_to_prev_proof.push(merkle_trees[round].generate_proof(query).unwrap());
            }
            let clone = queries_to_prev_ans.clone();
            let queries_to_prev: (Vec<Vec<F>>, Vec<Path<W::MerkleConfig>>) =
                (clone, queries_to_prev_proof);

            folded_evals_len = folded_evals_len / self.prover_config.folding_factor;
            query_indexes = dedup(query_indexes.into_iter().map(|i| i % folded_evals_len));

            let last_round_commitment_digest = if round == 0 {
                witness.commitment_digest()
            } else {
                commitments[round - 1].clone()
            };
            round_proofs.push(FRIProofRound {
                queries_to_prev,
                challenge_answers,
                challenge_values: queries_to_prev_ans,
                config: self.prover_config.clone(),
                last_round_commitment_digest,
            });
        }

        Self::Proof {
            config: self.prover_config.clone(),
            coeff: g_poly,
            commitment_digests: commitments,
            round_proofs,
            proof_of_work_nonce: proof_of_work(&mut sponge, self.prover_config.proof_of_work_bits),
            initial_commitment_digest: witness.commitment_digest(),
        }
    }
}
