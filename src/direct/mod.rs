pub mod config;
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
        direct::{config::DirectConfig, prover::DirectProver, verifier::DirectVerifier},
        ldt::Prover,
    };

    #[test]
    fn test_direct_ldt() {
        // config
        let mut rng = test_rng();
        let mt_config = merkle_tree::poseidon::default_config::<Field256>(&mut rng, 2);
        let fs_config = fs::poseidon::default_fs_config::<Field256>();
        let config = DirectConfig::new(
            22,
            8,
            mt_config.0.clone(),
            mt_config.1.clone(),
            fs_config.clone(),
        );

        // witness
        let polynomial = DensePolynomial::<Field256>::rand(config.degree, &mut rng);

        // commit
        let prover: DirectProver<
            Field256,
            merkle_tree::poseidon::MerkleTreeParams<Field256>,
            PoseidonSponge<Field256>,
        > = DirectProver::new(config.clone());
        let commitment = prover.commit(polynomial);

        // prove
        let direct_proof = prover.prove(commitment);

        // verify
        let verifier: DirectVerifier<
            Field256,
            merkle_tree::poseidon::MerkleTreeParams<Field256>,
            PoseidonSponge<Field256>,
        > = DirectVerifier::new(config);
        assert_eq!(verifier.verify(&direct_proof), true);
    }
}
