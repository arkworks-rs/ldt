pub mod config;
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
        fri::{prover::FRIProver, verifier::FRIVerifier, config::FRIConfig},
        ldt::{Prover, Verifier},
    };

    type TestField = Field256;
    type TestMerkleConfig = merkle_tree::poseidon::MerkleTreeParams<Field256>;
    type TestSpongeConfig = PoseidonSponge<Field256>;

    #[test]
    fn test_fri_ldt() {
        // get ready
        let mut rng = test_rng();
        let mt_config = merkle_tree::poseidon::default_config::<Field256>(&mut rng, 2);
        let fs_config = fs::poseidon::default_fs_config::<Field256>();
        let config: FRIConfig<TestMerkleConfig, TestSpongeConfig> = FRIConfig::new(
            2,
            8,
            4,
            mt_config.0.clone(),
            mt_config.1.clone(),
            8,
            4,
            fs_config.clone(),
            22,
            8,
        );

        // generate random witness
        let polynomial = DensePolynomial::<Field256>::rand(config.starting_degree, &mut rng);

        // commit
        let commitment = Commitment::<Field256, TestMerkleConfig>::new(
            config.starting_degree,
            config.starting_rate,
            config.folding_factor,
            config.merkle_leaf_hash_param.clone(),
            config.merkle_two_to_one_param.clone(),
            vec![polynomial],
        );

        // prove
        let prover: FRIProver<TestField, TestMerkleConfig, TestSpongeConfig> = FRIProver::new(config.clone());
        let fri_proof = prover.prove(&commitment);

        // verify
        let verifier: FRIVerifier<
            Field256,
            merkle_tree::poseidon::MerkleTreeParams<Field256>,
            PoseidonSponge<Field256>,
        > = FRIVerifier::new(config);
        assert_eq!(verifier.verify(&commitment, &fri_proof), true);
    }
}
