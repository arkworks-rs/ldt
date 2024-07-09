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
        commitment::{
            single::{SingleWitness, SingleWitnessArgument},
            Witness,
        },
        crypto::{fields::Field256, fs, merkle_tree},
        domain::Domain,
        fri::{config::FRIConfig, ldt::FRI},
        ldt::{LowDegreeTest, Prover, Verifier},
    };

    type TestField = Field256;
    type TestMerkleConfig = merkle_tree::poseidon::MerkleTreeParams<TestField>;
    type TestSpongeConfig = PoseidonSponge<TestField>;
    type TestWitness = SingleWitness<TestField, TestMerkleConfig, TestSpongeConfig>;

    #[test]
    fn test_fri_ldt() {
        // get ready
        let mut rng = test_rng();
        let (merkle_leaf_hash_param, merkle_two_to_one_param) =
            merkle_tree::poseidon::default_config::<Field256>(&mut rng, 2);
        let config: FRIConfig<TestMerkleConfig, TestSpongeConfig> = FRIConfig {
            folding_factor: 2,
            num_rounds: 4,
            num_queries: 8,
            merkle_leaf_hash_param: merkle_leaf_hash_param.clone(),
            merkle_two_to_one_param: merkle_two_to_one_param.clone(),
            proof_of_work_bits: 8,
            repetitions: 4,
            sponge_config: fs::poseidon::default_fs_config::<Field256>(),
            starting_degree: 22,
            starting_rate: 8,
        };
        let (prover, verifier) =
            FRI::<TestField, TestMerkleConfig, TestSpongeConfig, TestWitness>::new(config.clone());

        // generate witness
        let witness: SingleWitness<TestField, TestMerkleConfig, TestSpongeConfig> =
            SingleWitness::new(SingleWitnessArgument {
                coeff: DensePolynomial::<Field256>::rand(config.starting_degree, &mut rng),
                domain: Domain::<TestField>::new(config.starting_degree, config.starting_rate)
                    .unwrap(),
                folding_factor: config.folding_factor,
                merkle_leaf_hash_param,
                merkle_two_to_one_param,
                sponge_config: config.sponge_config,
            });

        // prove
        let fri_proof = prover.prove(&witness);

        // verify
        assert_eq!(verifier.verify(&fri_proof), true);
    }
}
