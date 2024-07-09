use std::collections::BTreeMap;

use ark_crypto_primitives::{
    merkle_tree::Config as MerkleConfig,
    sponge::{Absorb, CryptographicSponge},
};
use ark_ff::{batch_inversion, FftField, PrimeField};
use ark_poly::{EvaluationDomain, Polynomial};
use ark_std::marker::PhantomData;

use crate::{
    commitment::Witness,
    domain::Domain,
    fri::{config::FRIConfig, proof::FRIProof},
    ldt::Verifier,
    poly_utils,
    utils::{dedup, proof_of_work_verify, squeeze_integer},
};
pub struct FRIVerifier<F, M, S, W>
where
    F: FftField,
    M: MerkleConfig,
    S: CryptographicSponge,
    W: Witness<F, M>,
{
    verifier_config: FRIConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<W::MerkleConfig>,
    _sponge_config: PhantomData<S>,
}

impl<F, M, S, W> Verifier<F> for FRIVerifier<F, M, S, W>
where
    F: FftField + PrimeField,
    M: MerkleConfig<Leaf = Vec<F>>,
    M::InnerDigest: Absorb,
    S: CryptographicSponge,
    W: Witness<F, M>,
{
    type VerifierConfig = FRIConfig<M, S>;
    type Proof = FRIProof<F, M>;
    fn new(verifier_config: FRIConfig<M, S>) -> Self {
        Self {
            verifier_config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<W::MerkleConfig>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn verify(&self, proof: &Self::Proof) -> bool {
        // TODO fix this
        // if proof.final_polynomial.degree() + 1 > self.parameters.stopping_degree {
        //     return false;
        // }

        // We do FS
        let mut sponge = S::new(&self.verifier_config.sponge_config);
        sponge.absorb(&proof.initial_p_commitment_root);

        let mut folding_randomnessness: Vec<F> = vec![];
        folding_randomnessness.push(sponge.squeeze_field_elements(1)[0]);
        // Absorb the roots
        for commitment in &proof.commitments {
            sponge.absorb(&commitment);
            folding_randomnessness.push(sponge.squeeze_field_elements(1)[0]);
        }

        // We adjoin the initial commitment
        let commitments: Vec<_> = std::iter::once(proof.initial_p_commitment_root.clone())
            .chain(proof.commitments.iter().cloned())
            .collect();

        // Verify merkle commitments
        for num_round in 0..=self.verifier_config.num_rounds {
            let (answers, proofs) = &proof.round_proofs[num_round].queries_to_prev;
            for (i, proof) in proofs.iter().enumerate() {
                let answer = answers[i].clone();
                if !proof
                    .verify(
                        &self.verifier_config.merkle_leaf_hash_param,
                        &self.verifier_config.merkle_two_to_one_param,
                        &commitments[num_round],
                        answer,
                    )
                    .unwrap()
                {
                    return false;
                }
            }
        }

        let mut g_domain = Domain::<F>::new(
            self.verifier_config.starting_degree,
            self.verifier_config.starting_rate,
        )
        .unwrap();

        let mut folded_evals_len = g_domain.size() / self.verifier_config.folding_factor;
        let query_indexes = dedup(
            (0..self.verifier_config.repetitions)
                .map(|_| squeeze_integer(&mut sponge, folded_evals_len)),
        );

        // Precomputation
        let (generators, coset_offsets) = {
            let mut folded_evals_len = folded_evals_len;
            let mut query_indexes = query_indexes.clone();
            let mut generators = vec![];
            let mut coset_offsets = vec![];
            for _ in 0..=self.verifier_config.num_rounds {
                let generator =
                    g_domain.element(g_domain.size() / self.verifier_config.folding_factor);

                generators.push(generator);

                let round_offsets: Vec<_> =
                    query_indexes.iter().map(|i| g_domain.element(*i)).collect();
                coset_offsets.push(round_offsets);

                g_domain = g_domain.scale(self.verifier_config.folding_factor);
                folded_evals_len = folded_evals_len / self.verifier_config.folding_factor;
                query_indexes = dedup(query_indexes.into_iter().map(|i| i % folded_evals_len));
            }

            (generators, coset_offsets)
        };

        let size = F::from(self.verifier_config.folding_factor as u64);

        let mut to_invert: Vec<F> = vec![];

        for co in &coset_offsets {
            to_invert.extend(co);
        }
        to_invert.extend(&generators);
        to_invert.push(size);
        batch_inversion(&mut to_invert);
        let size_inv = to_invert.pop().unwrap();
        let generators_inv = to_invert.split_off(to_invert.len() - generators.len());
        let mut coset_offsets_inv = vec![];
        for co in coset_offsets.iter().rev() {
            let co_inv = to_invert.split_off(to_invert.len() - co.len());
            coset_offsets_inv.push(co_inv);
        }
        coset_offsets_inv.reverse();

        let mut query_indexes: Vec<_> = query_indexes.into_iter().map(|i| (i, 0)).collect();
        let mut folded_answers: Option<Vec<F>> = None;

        for num_round in 0..=self.verifier_config.num_rounds {
            let folding_randomness = folding_randomnessness[num_round];
            let answers: Vec<_> = query_indexes
                .iter()
                .zip(proof.round_proofs[num_round].queries_to_prev.0.clone())
                .map(|(index, answer)| (index, answer.clone()))
                .collect();

            if let Some(folded_answers) = &folded_answers {
                if !folded_answers.iter().zip(answers.iter()).all(
                    |(folded_answer, ((_, checking_index), answer))| {
                        answer[*checking_index] == *folded_answer
                    },
                ) {
                    return false;
                }
            }

            let generator = generators[num_round];
            let generator_inv = generators_inv[num_round];

            let unordeded_folded_answers: Vec<_> = answers
                .iter()
                .zip(coset_offsets[num_round].iter())
                .zip(coset_offsets_inv[num_round].iter())
                .map(|(((_, answer), coset_offset), coset_offset_inv)| {
                    let folded_answer = poly_utils::interpolation::fft_interpolate(
                        generator,
                        *coset_offset,
                        generator_inv,
                        *coset_offset_inv,
                        size_inv,
                        answer,
                    )
                    .evaluate(&folding_randomness);

                    folded_answer
                })
                .collect();

            folded_evals_len = folded_evals_len / self.verifier_config.folding_factor;

            // Now we need to sort and dedup
            let query_answers: BTreeMap<_, _> = query_indexes
                .into_iter()
                .zip(unordeded_folded_answers)
                .map(|((i, _), a)| (i % folded_evals_len, (i / folded_evals_len, a)))
                .collect();

            folded_answers = Some(query_answers.values().map(|(_, a)| *a).collect());

            query_indexes = query_answers
                .into_iter()
                .map(|(i, (c, _))| (i, c))
                .collect();
        }

        let folded_answers = folded_answers.unwrap();

        let answers: Vec<_> = query_indexes
            .into_iter()
            .map(|(index, checking_index)| {
                proof
                    .polynomial
                    .evaluate(&g_domain.element(index + checking_index * folded_evals_len))
            })
            .collect();
        if !folded_answers
            .iter()
            .zip(answers)
            .all(|(folded_answer, poly_answer)| poly_answer == *folded_answer)
        {
            return false;
        }

        // Proof of work
        proof_of_work_verify(
            &mut sponge,
            self.verifier_config.proof_of_work_bits,
            proof.proof_of_work_nonce,
        )
    }
}
