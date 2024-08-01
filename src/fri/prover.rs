use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{EvaluationDomain, Polynomial};
use ark_std::marker::PhantomData;

#[cfg(not(feature = "std"))]
use ark_std::{vec, vec::Vec};

use crate::{
    fri::{
        config::FRIConfig,
        proof::{FRIProof, FRIProofRound},
        prover_state::FRIProverState,
    },
    ldt::Prover,
    poly_utils,
    utils::{dedup, stack_evaluations},
    witness::Witness,
};

pub struct FRIProver<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
    W: Witness<F, M>,
{
    config: FRIConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge: PhantomData<S>,
    _witness: PhantomData<W>,
}

impl<F, M, S, W> Prover<F> for FRIProver<F, M, S, W>
where
    F: Absorb + FftField + PrimeField,
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

    fn new(config: FRIConfig<W::MerkleConfig, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<W::MerkleConfig>,
            _sponge: PhantomData::<S>,
            _witness: PhantomData::<W>,
        }
    }
    fn prove(&self, witness: &W) -> Self::Proof {
        assert!(witness.coeff().degree() < self.config.starting_degree);

        let mut state = FRIProverState::<F, M, S>::new(
            witness.coeff(),
            witness.commitment_digest(),
            witness.committed_values(),
            self.config.clone(),
            witness.domain(),
        );

        // Commit phase
        let mut commitments: Vec<M::InnerDigest> = vec![];
        let mut merkle_trees: Vec<MerkleTree<M>> = vec![witness.commitment()];
        let mut folded_evals: Vec<Vec<Vec<F>>> = vec![witness.committed_values()];

        let mut folding_randomness = state.sponge_squeeze();

        for _ in 0..self.config.num_rounds {
            let g_folded_evaluations =
                state.folded_evaluations(folding_randomness, folded_evals.last().unwrap().to_vec());

            state.domain = state.domain.scale(self.config.folding_factor);

            let g_merkle = MerkleTree::<W::MerkleConfig>::new(
                &self.config.merkle_leaf_hash_param,
                &self.config.merkle_two_to_one_param,
                &g_folded_evaluations,
            )
            .unwrap();
            let g_root = g_merkle.root();
            state.sponge_absorb(&g_root);

            folding_randomness = state.sponge_squeeze();

            commitments.push(g_root);
            merkle_trees.push(g_merkle);
            folded_evals.push(g_folded_evaluations);
        }

        // Query phase

        let mut folded_evals_len = witness.domain().size() / self.config.folding_factor;
        let mut challenges = state.challenges();

        // Note that we include final round as well
        let mut round_proofs = vec![];
        for round in 0..=self.config.num_rounds {
            let challenge_values: Vec<Vec<F>> = challenges
                .iter()
                .map(|&index| folded_evals[round][index].clone())
                .collect();
            let challenge_answers = merkle_trees[round]
                .generate_multi_proof(challenges.clone())
                .unwrap();

            folded_evals_len = folded_evals_len / self.config.folding_factor;
            challenges = dedup(challenges.into_iter().map(|i| i % folded_evals_len));

            let last_round_commitment_digest = if round == 0 {
                witness.commitment_digest()
            } else {
                commitments[round - 1].clone()
            };
            round_proofs.push(FRIProofRound {
                challenge_answers,
                challenge_values,
                config: self.config.clone(),
                last_round_commitment_digest,
            });
        }

        state.coeff = poly_utils::folding::poly_fold(
            &state.coeff,
            self.config.folding_factor,
            folding_randomness,
        );

        Self::Proof {
            coeff: state.coeff.clone(),
            rounds: round_proofs,
            proof_of_work_nonce: state.proof_of_work(),
        }
    }
}
