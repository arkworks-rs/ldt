use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{batch_inversion, FftField, PrimeField};
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
    utils::{dedup, proof_of_work_verify, squeeze_integer},
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
        let mut state = STIRVerifierState::new(self.config.clone(), claim.commitment_digest());

        for round_proof in &proof.rounds {
            if !round_proof.is_final_round {
                let round_result = self.round(round_proof, state);
                if round_result.is_none() {
                    return false;
                }
                state = round_result.unwrap();
            }
        }

        // Now, we sample the last points that we want to check consisntency at
        let final_repetitions = self.config.num_repetitions[self.config.num_rounds];
        let scaling_factor = state.domain_size / self.config.folding_factor;
        let final_randomness_indexes = dedup(
            (0..final_repetitions).map(|_| squeeze_integer(&mut state.sponge, scaling_factor)),
        );

        if !proof_of_work_verify(
            &mut state.sponge,
            self.config.num_proof_of_work_bits[self.config.num_rounds],
            proof.rounds.last().unwrap().proof_of_work_nonce,
        ) {
            return false;
        }

        // First, we want to query back the last oracle at this point, which is, again, just a
        // lookup
        let oracle_answers = proof.rounds.last().unwrap().challenge_values.clone();

        let folded_answers = state.folded_evaluations(final_randomness_indexes, oracle_answers);

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
        round_proof: &STIRProofRound<F, M, S>,
        mut state: STIRVerifierState<F, M, S>,
    ) -> Option<STIRVerifierState<F, M, S>> {
        // Redo FS
        state.sponge_absorb(&round_proof.commitment_digest);
        let ood_randomness = state.sponge_squeeze_multiple(self.config.num_out_of_domain_samples);
        state.sponge_absorb(&round_proof.out_of_domain_evaluations);
        let comb_randomness = state.sponge_squeeze();
        let new_folding_randomness = state.sponge_squeeze();
        let scaling_factor = state.domain_size / self.config.folding_factor;

        let num_repetitions = self.config.num_repetitions[state.round_num];
        let stir_randomness_indexes =
            dedup((0..num_repetitions).map(|_| squeeze_integer(&mut state.sponge, scaling_factor)));

        // PoW verification
        if !proof_of_work_verify(
            &mut state.sponge,
            self.config.num_proof_of_work_bits[state.round_num],
            round_proof.proof_of_work_nonce,
        ) {
            return None;
        }

        let shake_randomness = state.sponge_squeeze();

        // Now, we are starting to define the next function.
        // First, we need to query the previous oracle (which is either f_0 or g_i)
        // At the indexes B_i for i in stir_randomness_indexes
        // Since we previously verified the Merkle paths, this is easy
        // TODO: We should probably check the indexes
        let oracle_answers = round_proof.challenge_values.clone();

        // Now, for each of the selected random points, we need to compute the folding of the
        // previous oracle
        let folded_answers = state.folded_evaluations(stir_randomness_indexes, oracle_answers);

        // The quotient definining the function
        let quotient_answers: Vec<_> = ood_randomness
            .into_iter()
            .zip(&round_proof.out_of_domain_evaluations)
            .map(|(alpha, beta)| (alpha, *beta))
            .chain(folded_answers)
            .collect();
        let interpolating_polynomial = round_proof.coeff.clone();

        let ans_eval = interpolating_polynomial.evaluate(&shake_randomness);
        let shake_eval = round_proof.shake_coeff.evaluate(&shake_randomness);

        let mut denoms: Vec<_> = quotient_answers
            .iter()
            .map(|(x, _)| shake_randomness - x)
            .collect();

        batch_inversion(&mut denoms);
        // TODO: This maybe should be better
        if shake_eval
            != quotient_answers
                .iter()
                .zip(denoms)
                .map(|((_, y), d)| (ans_eval - y) * d)
                .sum()
        {
            return None;
        }

        let quotient_set = quotient_answers
            .into_iter()
            .map(|(x, _)| x)
            .collect::<Vec<_>>();

        Some(STIRVerifierState {
            comb_randomness: comb_randomness.clone(),
            config: self.config.clone(),
            domain_gen: state.domain_gen * state.domain_gen,
            domain_offset: state.domain_offset * state.domain_offset * state.root_of_unity,
            domain_size: state.domain_size / 2,
            folding_randomness: new_folding_randomness,
            interpolating_polynomial: interpolating_polynomial.clone(),
            quotient_set,
            root_of_unity: state.root_of_unity,
            round_num: state.round_num + 1,
            sponge: state.sponge,
        })
    }
}
