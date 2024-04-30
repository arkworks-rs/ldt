use ark_std::marker::PhantomData;

use crate::ldt::{LDTConfig, LowDegreeTest, Prover, Verifier};

// Config
pub struct DirectConfig<F, MerkleConfig, SpongeConfig> {
    pub ldt_config: LDTConfig,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<MerkleConfig>,
    _sponge_config: PhantomData<SpongeConfig>,
}
impl<Field, MerkleConfig, SpongeConfig> Prover<Field>
    for DirectConfig<Field, MerkleConfig, SpongeConfig>
{
}

// Prover
struct DirectProver<F, MerkleConfig, SpongeConfig> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<MerkleConfig>,
    _sponge_config: PhantomData<SpongeConfig>,
}
impl<Field, MerkleConfig, SpongeConfig> Prover<Field>
    for DirectProver<Field, MerkleConfig, SpongeConfig>
{
}

// Verifier
struct DirectVerifier<F, MerkleConfig, SpongeConfig> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<MerkleConfig>,
    _sponge_config: PhantomData<SpongeConfig>,
}
impl<Field, MerkleConfig, SpongeConfig> Verifier<Field>
    for DirectVerifier<Field, MerkleConfig, SpongeConfig>
{
}

// LDT
struct DirectLDT<F, MerkleConfig, SpongeConfig> {
    _field: PhantomData<F>,
    _merkle_config: PhantomData<MerkleConfig>,
    _sponge_config: PhantomData<SpongeConfig>,
}
impl<F, MerkleConfig, SpongeConfig> LowDegreeTest<F>
    for DirectLDT<F, MerkleConfig, SpongeConfig>
{
    type Config = DirectConfig<F, MerkleConfig, SpongeConfig>;
    type Prover = DirectProver<F, MerkleConfig, SpongeConfig>;
    type Verifier = DirectVerifier<F, MerkleConfig, SpongeConfig>;
    fn config(ldt_config: LDTConfig) -> Self::Config {
        DirectConfig {
            ldt_config,
            _field: PhantomData,
            _merkle_config: PhantomData,
            _sponge_config: PhantomData
        }
    }
    fn prover(config: Self::Config) -> Self::Prover {
        Self::Prover {
            _field: PhantomData,
            _merkle_config: PhantomData,
            _sponge_config: PhantomData
        }
    }
    fn verifier(config: Self::Config) -> Self::Verifier {
        Self::Verifier {
            _field: PhantomData,
            _merkle_config: PhantomData,
            _sponge_config: PhantomData
        }
    }
}

#[cfg(test)]
mod tests {
    // use ark_crypto_primitives::{
    //     merkle_tree::Config as MerkleConfig,
    //     sponge::CryptographicSponge as Sponge
    // };
    // use ark_ff::FftField;
    use crate::direct::{DirectProver, DirectVerifier, DirectLDT};
    use crate::ldt::{LDTConfig, LowDegreeTest};
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_std::test_rng;
    use ark_bn254::Fr as BN254;

    #[test]
    fn test_direct_ldt() {
        // config
        let degree = 22;
        let rate = 4;
        let config = DirectLDT::config(LDTConfig::new(degree, rate));
        
        // witness
        let mut rng = test_rng();
        let witness_polynomial = DensePolynomial::<BN254>::rand(degree, &mut rng);

        // commit and prove
        let prover = DirectLDT::prover(config);
        let (commitment, witness) = prover.commit(witness_polynomial);
        let proof = prover.prove(witness);

        // verify
        let verifier = DirectLDT::verifier(config);
        verifier.verify(commitment, proof);
    }
}