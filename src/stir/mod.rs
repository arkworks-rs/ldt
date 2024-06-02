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
        ldt::{LowDegreeTest, Prover, Verifier},
        stir::{config::STIRConfig, ldt::STIR},
    };

    type TestField = Field256;
    type TestMerkleConfig = merkle_tree::poseidon::MerkleTreeParams<Field256>;
    type TestSpongeConfig = PoseidonSponge<Field256>;

    #[test]
    fn test_stir_ldt() {
        // get ready
        let mut rng = test_rng();
        let mt_config = merkle_tree::poseidon::default_config::<Field256>(&mut rng, 2);
        let fs_config = fs::poseidon::default_fs_config::<Field256>();
        let config: STIRConfig<TestMerkleConfig, TestSpongeConfig> = STIRConfig::new(
            2,
            4,
            mt_config.0.clone(),
            mt_config.1.clone(),
            2,
            vec![2, 2, 2, 2, 2], // TODO: these are getting 1-indexed for some reason
            vec![2, 2, 2, 2, 2],
            fs_config.clone(),
            22,
            8,
        );
        let (prover, verifier) =
            STIR::<TestField, TestMerkleConfig, TestSpongeConfig>::new(config.clone());

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
        let stir_proof = prover.prove(&commitment);

        // verify
        assert_eq!(verifier.verify(&commitment, &stir_proof), true);
    }
}
