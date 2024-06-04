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
        commitment::Commitment,
        crypto::{fields::Field256, fs, merkle_tree},
        direct::{config::DirectConfig, ldt::DirectLDT},
        ldt::{LowDegreeTest, Prover, Verifier},
    };

    type TestField = Field256;
    type TestMerkleConfig = merkle_tree::poseidon::MerkleTreeParams<TestField>;
    type TestSpongeConfig = PoseidonSponge<TestField>;

    #[test]
    fn test_direct_ldt() {
        // get ready
        let mut rng = test_rng();
        let (merkle_leaf_hash_param, merkle_two_to_one_param) =
            merkle_tree::poseidon::default_config::<Field256>(&mut rng, 2);
        let config: DirectConfig<TestMerkleConfig, TestSpongeConfig> = DirectConfig {
            degree: 22,
            num_queries: 2,
            merkle_leaf_hash_param,
            merkle_two_to_one_param,
            sponge_config: fs::poseidon::default_fs_config::<Field256>(),
        };
        let (prover, verifier) =
            DirectLDT::<TestField, TestMerkleConfig, TestSpongeConfig>::new(config.clone());

        // generate random witness
        let polynomial = DensePolynomial::<Field256>::rand(config.degree, &mut rng);

        // commit
        let commitment = Commitment::<TestField, TestMerkleConfig>::new(
            config.degree.clone(),
            0,
            1,
            config.merkle_leaf_hash_param.clone(),
            config.merkle_two_to_one_param.clone(),
            vec![polynomial],
        );

        // prove
        let direct_proof = prover.prove(&commitment);

        // verify
        assert_eq!(verifier.verify(&commitment, &direct_proof), true);
    }
}
