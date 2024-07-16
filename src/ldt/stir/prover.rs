use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, MerkleTree},
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

        // Step 1: initial state
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
            self.config.num_rounds,
            self.config.sponge_config.clone(),
            witness.coeff(),
        );

        // Step 2: compute inner rounds
        let mut round_proofs = Vec::with_capacity(self.config.num_rounds);
        for _round in 0..self.config.num_rounds {
            round_state.next();
            // let new_round_state = Self::compute_inner_round(&self.config, round_state);
            // round_state = new_round_state;
            round_proofs.push(round_state.round_proof());
        }

        // Step 3: compute final round (v similar but fewer things)
        if round_state.is_final_round() {
            let final_round_proof = Self::compute_final_round(&self.config, round_state);
            round_proofs.push(final_round_proof.clone());
        }

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
        round_state.fold();

        // Step 2: Generate challenges and answers
        round_state.update_challenges();

        // Step 3: Proof of work
        round_state.update_proof_of_work();

        // Step 4: Return
        round_state.round_proof()
    }
}
