use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree, Path},
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_std::marker::PhantomData;

use crate::{
    ldt::{
        stir::{
            config::STIRConfig,
            proof::{STIRProof, STIRRoundProof},
            state::STIRRoundState,
        },
        Prover,
    },
    poly_utils,
    utils::{dedup, proof_of_work, squeeze_integer},
    witness::Witness,
};

pub struct STIRProver<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    W: Witness<F, M, MerkleConfig = M>,
{
    config: STIRConfig<W::MerkleConfig, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<W::MerkleConfig>,
    _sponge_config: PhantomData<S>,
}

impl<F, M, S, W> Prover<F> for STIRProver<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    S::Config: Clone,
    W: Witness<F, M, MerkleConfig = M, Commitment = MerkleTree<M>, CommittedValues = Vec<Vec<F>>>
        + Clone,
    W::ChallengeAnswers: Clone,
{
    type Witness = W;
    type ProverConfig = STIRConfig<M, S>;
    type Proof = STIRProof<F, M>;

    fn new(config: STIRConfig<M, S>) -> Self {
        Self {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }

    fn prove(&self, witness: &W) -> Self::Proof {
        assert!(witness.coeff().degree() < self.config.starting_degree);

        // Step 1: initial round state from witness
        let mut round_state = STIRRoundState::new(
            witness.domain(),
            witness.commitment(),
            witness.committed_values(),
            self.config.folding_factor,
            self.config.merkle_leaf_hash_param.clone(),
            self.config.merkle_two_to_one_param.clone(),
            self.config.num_out_of_domain_samples,
            self.config.proof_of_work_bits.clone(),
            self.config.repetitions.clone(),
            self.config.sponge_config.clone(),
            witness.coeff(),
        );

        // Step 2: compute inner rounds
        let mut round_proofs = Vec::with_capacity(self.config.num_rounds);
        for _round in 0..self.config.num_rounds {
            round_state.next_state();
            // let new_round_state = Self::compute_inner_round(&self.config, round_state);
            // round_state = new_round_state;
            round_proofs.push(round_state.round_proof());
        }

        // Step 3: compute final round (v similar but fewer things)
        let final_round_proof = Self::compute_final_round(&self.config, round_state);
        round_proofs.push(final_round_proof.clone());

        // Boom.
        STIRProof::<F, M> { round_proofs }
    }
}

impl<F, M, S, W> STIRProver<F, M, S, W>
where
    F: FftField + PrimeField + Absorb,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    W: Witness<F, M, MerkleConfig = M>,
{
    fn compute_final_round(
        config: &STIRConfig<W::MerkleConfig, S>,
        mut round_state: STIRRoundState<F, W::MerkleConfig, S>,
    ) -> STIRRoundProof<F, W::MerkleConfig> {
        // Step 1: Perfom fold operation
        let coeff = poly_utils::folding::poly_fold(
            &round_state.witness_coeff,
            config.folding_factor,
            round_state.folding_randomness,
        );

        // Step 2: Generate challenges and answers
        // round_state.update_challenges();
        let tmp_round_num = round_state.round_num;
        let tmp_domain = round_state.domain.clone();
        let challenges = Self::challenges(
            &mut round_state,
            tmp_domain.size() / config.folding_factor,
            config.repetitions[tmp_round_num],
        );
        let (committed_values, challenge_answers) = Self::challenge_answers(
            challenges,
            round_state.commitment.clone(),
            round_state.committed_values,
        );

        // Step 3: Proof of work
        let proof_of_work_nonce = proof_of_work(
            &mut round_state.sponge,
            config.proof_of_work_bits[config.num_rounds],
        );

        // Step 4: Return
        STIRRoundProof::new(
            coeff,
            challenge_answers,
            committed_values,
            true,
            vec![],
            round_state.commitment.root(),
            proof_of_work_nonce,
            DensePolynomial::<F>::from_coefficients_vec(Vec::new()),
        )
    }
    fn challenges(
        round_state: &mut STIRRoundState<F, M, S>,
        scaling_factor: usize,
        num_repetitions: usize,
    ) -> Vec<usize> {
        dedup(
            (0..num_repetitions).map(|_| squeeze_integer(&mut round_state.sponge, scaling_factor)),
        )
    }
    fn challenge_answers(
        challenges: Vec<usize>,
        last_round_commitment: MerkleTree<M>,
        last_round_committed_values: Vec<Vec<F>>,
    ) -> (Vec<Vec<F>>, Vec<Path<M>>) {
        let challenge_values: Vec<Vec<F>> = challenges
            .iter()
            .map(|index| last_round_committed_values[*index].clone())
            .collect();
        let mut challenge_answers: Vec<Path<M>> = Vec::with_capacity(challenge_values.len());
        for challenge in challenges {
            challenge_answers.push(last_round_commitment.generate_proof(challenge).unwrap());
        }
        (challenge_values, challenge_answers)
    }
}
