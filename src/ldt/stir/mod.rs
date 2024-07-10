pub mod config;
pub mod ldt;
pub mod proof;
pub mod prover;
pub mod verifier;

#[cfg(test)]
mod tests {
    use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_std::test_rng;

    use crate::{
        crypto::{fields::Field256, fs, merkle_tree},
        domain::Domain,
        ldt::{
            stir::{config::STIRConfig, ldt::STIR},
            LowDegreeTest, Prover, Verifier,
        },
        witness::{
            single::{SingleWitness, SingleWitnessArgument},
            Witness,
        },
    };

    type TestField = Field256;
    type TestMerkleConfig = merkle_tree::poseidon::MerkleTreeParams<TestField>;
    type TestSpongeConfig = PoseidonSponge<TestField>;
    type TestWitness = SingleWitness<TestField, TestMerkleConfig, TestSpongeConfig>;

    #[test]
    fn test_stir_ldt() {
        // get ready
        let mut rng = test_rng();
        let (merkle_leaf_hash_param, merkle_two_to_one_param) =
            merkle_tree::poseidon::default_config::<Field256>(&mut rng, 2);
        let config: STIRConfig<TestMerkleConfig, TestSpongeConfig> = STIRConfig {
            folding_factor: 16,
            num_rounds: 4,
            merkle_leaf_hash_param: merkle_leaf_hash_param.clone(),
            merkle_two_to_one_param: merkle_two_to_one_param.clone(),
            num_out_of_domain_samples: 2,
            proof_of_work_bits: vec![2, 2, 2, 2, 2],
            repetitions: vec![2, 2, 2, 2, 2],
            sponge_config: fs::poseidon::default_fs_config::<Field256>(),
            starting_degree: 16,
            starting_rate: 8,
            stopping_degree: 8,
        };
        let (prover, verifier) =
            STIR::<TestField, TestMerkleConfig, TestSpongeConfig, TestWitness>::new(config.clone());

        let witness: SingleWitness<TestField, TestMerkleConfig, TestSpongeConfig> =
            SingleWitness::new(SingleWitnessArgument {
                coeff: DensePolynomial::<Field256>::rand(config.starting_degree, &mut rng),
                domain: Domain::<TestField>::new(config.starting_degree, config.starting_rate)
                    .unwrap(),
                folding_factor: 16,
                merkle_leaf_hash_param,
                merkle_two_to_one_param,
                sponge_config: config.sponge_config,
            });

        // prove
        let stir_proof = prover.prove(&witness);

        // verify
        assert_eq!(verifier.verify(&witness.claim(), &stir_proof), true);
    }
}
