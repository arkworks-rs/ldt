use ark_crypto_primitives::{merkle_tree::{Config, MerkleTree, Path}, sponge::poseidon::PoseidonConfig};
use ark_ff::FftField;
use ark_poly::{EvaluationDomain, univariate::DensePolynomial};
use ark_std::{marker::PhantomData, test_rng};

use crate::{crypto::merkle_tree, domain::Domain, ldt::{LDTConfig, LowDegreeTest, Prover, Verifier}};

pub struct DirectProof<F: FftField, MerkleConfig: Config> {
    polynomial: Vec<F>,
    openings: Vec<Path<MerkleConfig>>
}

// Config
pub struct DirectConfig<F: FftField, MerkleConfig: Config, SpongeConfig> {
    pub ldt_config: LDTConfig<F>,
    pub log_root_of_unity: usize,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<MerkleConfig>,
    _sponge_config: PhantomData<SpongeConfig>,
}
impl<F: FftField, MerkleConfig: Config, SpongeConfig> DirectConfig<F, MerkleConfig, SpongeConfig> {}

// Prover
struct DirectProver<F: FftField, MerkleConfig: Config, SpongeConfig> {
    config: DirectConfig<F, MerkleConfig, SpongeConfig>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<MerkleConfig>,
    _sponge_config: PhantomData<SpongeConfig>,
}
impl<F: FftField + ark_ff::PrimeField + ark_crypto_primitives::sponge::Absorb, MerkleConfig: Config, SpongeConfig> Prover<F> for DirectProver<F, MerkleConfig, SpongeConfig> {
    fn commit(&self, witness_polynomial: DensePolynomial<F>) { // TODO -> MerkleConfig::InnerDigest {
        // config
        let degree = 22;
        let log_root_of_unity = 0;
        let num_queries = 8;
        let mut rng = test_rng();
        let domain = Domain::<F>::new(degree, log_root_of_unity).unwrap();
        let evals: Vec<Vec<F>> = witness_polynomial
            .evaluate_over_domain_by_ref(domain.backing_domain)
            .evals
            .iter()
            .map(|f| -> Vec<F> { vec![*f] })
            .collect();
        // commit
        let (leaf_hash_params, two_to_one_params) = merkle_tree::poseidon::default_config(&mut rng, 2);
        let mt = MerkleTree::<MerkleConfig>::new(
            &leaf_hash_params,
            &two_to_one_params,
            &evals,
        )
        .unwrap();
        let commitment = mt.root(); // TODO
    }
    fn prove(&self, witness_polynomial: DensePolynomial<F>) {
        // let domain =
        //     Domain::<F>::new(self.config.ldt_config.degree, self.config.log_root_of_unity).unwrap();
        // let evals = domain.fft(&witness_polynomial);
        // let mut rng = test_rng(); // TODO: where should this randomness live? Config?
        // let (leaf_hash_params, two_to_one_params) =
        //     MerkleConfig::default_config(&mut rng, 2);
        // let merkle_tree: MerkleTree::<MerkleConfig> =
        // MerkleTree::<MerkleConfig>::new(&leaf_hash_params, &two_to_one_params, &evals).unwrap();
        // let commitment = merkle_tree.root();
    }
}

// Verifier
struct DirectVerifier<F, MerkleConfig: Config, SpongeConfig> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<MerkleConfig>,
    _sponge_config: PhantomData<SpongeConfig>,
}
impl<F: FftField, MerkleConfig: Config, SpongeConfig> Verifier<F>
    for DirectVerifier<F, MerkleConfig, SpongeConfig>
{
}

// LDT
struct DirectLDT<F, MerkleConfig, SpongeConfig> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<MerkleConfig>,
    _sponge_config: PhantomData<SpongeConfig>,
}
impl<F: FftField, MerkleConfig: Config, SpongeConfig> LowDegreeTest<F> for DirectLDT<F, MerkleConfig, SpongeConfig> {
    type Proof = DirectProof<F, MerkleConfig>;
    type Config = DirectConfig<F, MerkleConfig, SpongeConfig>;
    type Prover = DirectProver<F, MerkleConfig, SpongeConfig>;
    type Verifier = DirectVerifier<F, MerkleConfig, SpongeConfig>;
    fn config(ldt_config: LDTConfig<F>) -> Self::Config {
        DirectConfig {
            ldt_config,
            log_root_of_unity: 0, // TODO
            _field: PhantomData,
            _merkle_config: PhantomData,
            _sponge_config: PhantomData,
        }
    }
    fn prover(config: Self::Config) -> Self::Prover {
        Self::Prover {
            config,
            _field: PhantomData,
            _merkle_config: PhantomData,
            _sponge_config: PhantomData,
        }
    }
    fn verifier(config: Self::Config) -> Self::Verifier {
        Self::Verifier {
            _field: PhantomData,
            _merkle_config: PhantomData,
            _sponge_config: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::crypto;
    use crate::crypto::merkle_tree::blake3::{Blake3LeafHash, Blake3TwoToOneCRHScheme};
    use crate::crypto::{fields::Field256, fs, merkle_tree};
    use crate::direct::{DirectConfig, DirectLDT, DirectProver, DirectVerifier};
    use crate::domain::Domain;
    use crate::ldt::{LDTConfig, LowDegreeTest};
    use crate::utils::squeeze_integer;
    use ark_crypto_primitives::merkle_tree::{MerkleTree, Path};
    use ark_crypto_primitives::sponge::poseidon::{PoseidonConfig, PoseidonSponge};
    use ark_crypto_primitives::sponge::CryptographicSponge;
    use ark_ff::Field;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_poly::EvaluationDomain;
    use ark_std::test_rng;

    #[test]
    fn test_direct_ldt() {
        // config
        let degree = 22;
        // let rate = 4;
        let log_root_of_unity = 0;
        let num_queries = 8;

        // witness
        let mut rng = test_rng();
        let witness_polynomial = DensePolynomial::<Field256>::rand(degree, &mut rng);
        let domain = Domain::<Field256>::new(degree, log_root_of_unity).unwrap();
        let evals: Vec<Vec<Field256>> = witness_polynomial
            .evaluate_over_domain_by_ref(domain.backing_domain)
            .evals
            .iter()
            .map(|f| -> Vec<Field256> { vec![*f] })
            .collect();

        // commit
        let (leaf_hash_params, two_to_one_params): (
            PoseidonConfig<Field256>,
            PoseidonConfig<Field256>,
        ) = merkle_tree::poseidon::default_config(&mut rng, 2);
        let mt = MerkleTree::<merkle_tree::poseidon::MerkleTreeParams<Field256>>::new(
            &leaf_hash_params,
            &two_to_one_params,
            &evals,
        )
        .unwrap();
        let commitment = mt.root();

        // prove
        let mut fs1: PoseidonSponge<Field256> = fs::poseidon::Sponge::new(&leaf_hash_params);
        fs1.absorb(&commitment);
        let mut queries1: Vec<usize> = Vec::with_capacity(num_queries);
        for _ in 0..num_queries {
            queries1.push(squeeze_integer(&mut fs1, 32));
        }
        let mut auth: Vec<Path<merkle_tree::poseidon::MerkleTreeParams<Field256>>> =
            Vec::with_capacity(num_queries);
        for query in queries1 {
            auth.push(mt.generate_proof(query).unwrap());
        }

        // verify
        let mut fs2: PoseidonSponge<Field256> = fs::poseidon::Sponge::new(&leaf_hash_params);
        fs2.absorb(&commitment);
        for i in 0..num_queries {
            let query = squeeze_integer(&mut fs2, 32);
            // is correct query index
            assert_eq!(auth[i].leaf_index, query);
            // is correct path
            assert_eq!(
                auth[i]
                    .verify(
                        &leaf_hash_params,
                        &two_to_one_params,
                        &commitment,
                        evals[query].clone()
                    )
                    .unwrap(),
                true
            );
        }

        // SKELETON
        // let config = DirectLDT::config(LDTConfig::new(degree, rate));

        // // witness
        // let mut rng = test_rng();
        // let witness_polynomial = DensePolynomial::<BN254>::rand(degree, &mut rng);

        // // commit and prove
        // let prover = DirectLDT::prover(config);
        // let (commitment, witness) = prover.commit(witness_polynomial);
        // let proof = prover.prove(witness);

        // // verify
        // let verifier = DirectLDT::verifier(config);
        // verifier.verify(commitment, proof);
    }
}
