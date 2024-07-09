pub mod config;
pub mod ldt;
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
        direct::{config::DirectConfig, ldt::DirectLDT},
        domain::Domain,
        ldt::{LowDegreeTest, Prover, Verifier},
    };

    type TestField = Field256;
    type TestMerkleConfig = merkle_tree::poseidon::MerkleTreeParams<TestField>;
    type TestSpongeConfig = PoseidonSponge<TestField>;
    type TestWitness = SingleWitness<TestField, TestMerkleConfig, TestSpongeConfig>;

    #[test]
    fn test_direct_ldt() {
        // get ready
        let mut rng = test_rng();
        let (merkle_leaf_hash_param, merkle_two_to_one_param) =
            merkle_tree::poseidon::default_config::<Field256>(&mut rng, 2);
        let config: DirectConfig<TestMerkleConfig, TestSpongeConfig> = DirectConfig {
            degree: 22,
            num_challenges: 2,
            merkle_leaf_hash_param: merkle_leaf_hash_param.clone(),
            merkle_two_to_one_param: merkle_two_to_one_param.clone(),
            sponge_config: fs::poseidon::default_fs_config::<Field256>(),
        };
        let (prover, verifier) =
            DirectLDT::<TestField, TestMerkleConfig, TestSpongeConfig, TestWitness>::new(
                config.clone(),
            );

        // generate witness
        let witness: SingleWitness<TestField, TestMerkleConfig, TestSpongeConfig> =
            SingleWitness::new(SingleWitnessArgument {
                coeff: DensePolynomial::<Field256>::rand(config.degree, &mut rng),
                domain: Domain::<TestField>::new(config.degree, 0).unwrap(),
                folding_factor: 1,
                merkle_leaf_hash_param,
                merkle_two_to_one_param,
                sponge_config: config.sponge_config,
            });

        // prove
        let direct_proof = prover.prove(&witness);

        // verify
        assert_eq!(verifier.verify(&direct_proof), true);
    }
}
