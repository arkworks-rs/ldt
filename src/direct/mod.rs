use ark_crypto_primitives::merkle_tree::{MerkleTree, Config};
use ark_ff::FftField;
use ark_poly::{EvaluationDomain, univariate::DensePolynomial};
use ark_std::{marker::PhantomData, test_rng};

use crate::{domain::Domain, ldt::{LDTConfig, LowDegreeTest, Prover, Verifier}};

pub struct DirectProof {}

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
impl<F: FftField, MerkleConfig: Config, SpongeConfig> Prover<F> for DirectProver<F, MerkleConfig, SpongeConfig> {
    fn prove(&self, witness_polynomial: DensePolynomial<F>) {
        let domain =
            Domain::<F>::new(self.config.ldt_config.degree, self.config.log_root_of_unity).unwrap();
        let evals = domain.fft(&witness_polynomial);
        let mut rng = test_rng(); // TODO: where should this randomness live? Config?
        let (leaf_hash_params, two_to_one_params) =
            MerkleConfig::default_config(&mut rng, 2);
        let merkle_tree: MerkleTree::<MerkleConfig> =
        MerkleTree::<MerkleConfig>::new(&leaf_hash_params, &two_to_one_params, &evals).unwrap();
        let commitment = merkle_tree.root();
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
    type Proof = DirectProof;
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
    // use ark_crypto_primitives::{
    //     merkle_tree::Config as MerkleConfig,
    //     sponge::CryptographicSponge as Sponge
    // };
    // use ark_ff::FftField;
    use crate::crypto::{fields::Field256, fs, merkle_tree};
    use crate::direct::{DirectConfig, DirectLDT, DirectProver, DirectVerifier};
    use crate::domain::Domain;
    use crate::ldt::{LDTConfig, LowDegreeTest};
    use ark_crypto_primitives::merkle_tree::MerkleTree;
    use ark_crypto_primitives::sponge::poseidon::PoseidonConfig;
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
        let rate = 4;
        let ldt_config = LDTConfig::new(degree, rate);
        let config: DirectConfig<Field256, MerkleTree<crypto::merkle_tree::poseidon::MerkleTreeParams<Field256>>, PoseidonConfig<Field256>> =
            DirectLDT::config(ldt_config);

        // commit and prove
        let prover: DirectProver<Field256, MerkleTree<merkle_tree::poseidon::MerkleTreeParams<Field256>>, PoseidonConfig<Field256>> =
            DirectLDT::prover(config);
        let mut rng = test_rng();
        let witness_polynomial = DensePolynomial::<Field256>::rand(degree, &mut rng);

        
        let (commitment, witness) = prover.commit(witness_polynomial);
        let proof = prover.prove(witness);

        let num_queries = 8;
        let domain =
            Domain::<F>::new(self.config.ldt_config.degree, self.config.log_root_of_unity).unwrap();
        let evals = domain.fft(&witness_polynomial);
        let mut rng = test_rng(); // TODO: where should this randomness live? Config?
        let (leaf_hash_params, two_to_one_params) =
            MerkleConfig::default_config(&mut rng, 2);
        let merkle_tree: MerkleTree::<MerkleConfig> =
        MerkleTree::<MerkleConfig>::new(&leaf_hash_params, &two_to_one_params, &evals).unwrap();
        let commitment = merkle_tree.root();
        // witness
        let mut rng = test_rng();
        let witness_polynomial = DensePolynomial::<Field256>::rand(degree, &mut rng);
        let evals = domain.fft(&witness_polynomial);
        let (leaf_hash_params, two_to_one_params) =
            merkle_tree::poseidon::default_config(&mut rng, 2);
        let merkle_tree: MerkleTree<merkle_tree::poseidon::MerkleTreeParams<Field256>> =
            MerkleTree::new(&leaf_hash_params, &two_to_one_params, &evals).unwrap();
        let commitment = merkle_tree.root();

        let fs_config: PoseidonConfig<Field256> = fs::poseidon::default_fs_config();
        let fs = fs::poseidon::Sponge::new(&fs_config);
        for eval in evals {
            // Absorb the leaf data into the sponge
            fs.absorb(&eval);
            // Store the query
        }

        // Squeeze out the query result
        let query: Vec<Field256> = fs.squeeze_field_elements(num_queries);

        let mut openings = Vec::new();
        for leaf_index in 0..evals.len() {
            // Get the leaf value
            let leaf_value = merkle_tree.leaf_at(leaf_index).unwrap();
            // Generate the authentication path (opening) for the leaf
            let opening = merkle_tree.generate_proof(leaf_index).unwrap();
            // Store the leaf value and its corresponding opening
            openings.push((leaf_value, opening));
        }
        let auth = merkle_tree.multiopen(query);

        // commit and prove
        let prover = DirectLDT::prover(config);
        let (commitment, witness) = prover.commit(witness_polynomial);
        let proof = prover.prove(witness);

        // verify
        let verifier = DirectLDT::verifier(config);
        verifier.verify(commitment, proof);
    }
}
