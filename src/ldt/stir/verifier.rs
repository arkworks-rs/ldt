use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::Polynomial;
use ark_std::marker::PhantomData;

use crate::{
    ldt::{
        stir::{
            config::STIRConfig,
            proof::{STIRProof, STIRProofRound},
            verifier_state::STIRVerifierState,
        },
        Verifier,
    },
    statement::single::SingleStatement,
    witness::Witness,
};

pub struct STIRVerifier<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
    W: Witness<F, M, MerkleConfig = M>,
{
    config: STIRConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
    _witness: PhantomData<W>,
}
impl<F, M, S, W> Verifier<F> for STIRVerifier<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M> + Clone,
    W::ChallengeAnswers: Clone,
{
    type Statement = SingleStatement<M>;
    type VerifierConfig = STIRConfig<M, S>;
    type Proof = STIRProof<F, M, S>;

    fn new(config: STIRConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
            _witness: PhantomData::<W>,
        }
    }
    fn verify(&self, claim: &Self::Statement, proof: &Self::Proof) -> bool {
        // if proof.final_round_proof.coeff.degree() + 1 > self.config.stopping_degree { // TODO: fix this
        //     return false;
        // }

        // Step 1: Verify merkle paths for all rounds
        if !proof.verify_challenge_answers() {
            return false;
        }

        // Step 2: Recompute
        let mut state = STIRVerifierState::new(
            self.config.clone(),
            claim.commitment_digest(),
            proof.clone(),
        );
        for round_proof in &proof.rounds {
            if !round_proof.is_final_round {
                let round_result = self.round(proof, round_proof, state);
                if round_result.is_none() {
                    return false;
                }
                state = round_result.unwrap();
            }
        }

        // Step 3: Randomness indices
        let randomness_indices: Vec<usize> = state.randomness_indices();

        // Step 4: Proof of work
        if !state.verify_proof_of_work(proof) {
            return false;
        }

        // Step 4: Folded answers
        let oracle_answers = proof.rounds.last().unwrap().challenge_values.clone();
        let folded_answers = state.folded_evaluations(randomness_indices, oracle_answers);
        folded_answers
            .into_iter()
            .all(|(point, value)| proof.rounds.last().unwrap().coeff.evaluate(&point) == value)
    }
}

impl<F, M, S, W> STIRVerifier<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M> + Clone,
    W::ChallengeAnswers: Clone,
{
    fn round(
        &self,
        proof: &STIRProof<F, M, S>,
        round_proof: &STIRProofRound<F, M, S>,
        mut state: STIRVerifierState<F, M, S>,
    ) -> Option<STIRVerifierState<F, M, S>> {
        // Step 1: handle randomness
        let (out_of_domain_randomness, comb_randomness, folding_randomness, randomness_indices) =
            state.randomness(
                round_proof.commitment_digest.clone(),
                round_proof.out_of_domain_evaluations.clone(),
            );

        // Step 2: proof of work
        if !state.verify_proof_of_work(proof) {
            return None;
        }

        // Step 3: more randomness
        let shake_randomness = state.sponge_squeeze();

        // Step 4: quotient answers
        let quotient_answers: Vec<(F, F)> = state.quotient_answers(
            &round_proof.challenge_values,
            &out_of_domain_randomness,
            &round_proof.out_of_domain_evaluations,
            &randomness_indices,
        );

        // Step 5: verify quotient answers
        if !round_proof.verify_quotient_answers(&quotient_answers, &shake_randomness) {
            return None;
        }

        Some(STIRVerifierState {
            comb_randomness: comb_randomness.clone(),
            config: self.config.clone(),
            domain_gen: state.domain_gen * state.domain_gen,
            domain_offset: state.domain_offset * state.domain_offset * state.root_of_unity,
            domain_size: state.domain_size / 2,
            folding_randomness: folding_randomness,
            interpolating_coeff: round_proof.coeff.clone(),
            proof: state.proof,
            quotient_set: quotient_answers.into_iter().map(|(x, _)| x).collect(),
            root_of_unity: state.root_of_unity,
            round_num: state.round_num + 1,
            sponge: state.sponge,
        })
    }
}
