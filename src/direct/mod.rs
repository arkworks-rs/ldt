use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, MerkleTree, Path, TwoToOneParam}, sponge::{Absorb, CryptographicSponge}
};
use ark_ff::FftField;
use ark_poly::univariate::DensePolynomial;
use ark_std::marker::PhantomData;

use crate::{domain::Domain, utils::squeeze_integer};

pub struct DirectConfig<M: MerkleConfig, S: CryptographicSponge> {
    pub degree: usize,
    pub num_queries: usize,
    pub merkle_leaf_hash_param: LeafParam<M>,
    pub merkle_two_to_one_param: TwoToOneParam<M>,
    pub sponge_config: S::Config,
}

// pub struct DirectProof<F: FftField, M: MerkleConfig> {
//     polynomial: Vec<F>,
//     openings: Vec<Path<M>>,
// }

// Prover
struct DirectProver<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    config: DirectConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> DirectProver<F, M, S>
where
    M::InnerDigest: Absorb,
 {
    fn new(config: DirectConfig<M, S>) -> Self {
        DirectProver {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn commit(&self, polynomial: DensePolynomial<F>) -> (Vec<Vec<F>>, MerkleTree<M>) {
        // get evaluations over a domain
        let domain = Domain::<F>::new(self.config.degree, 0).unwrap();
        let evals: Vec<Vec<F>> = polynomial
            .evaluate_over_domain_by_ref(domain.backing_domain)
            .evals
            .iter()
            .map(|f| -> Vec<F> { vec![*f] })
            .collect();
        // generate the committment
        let mt = MerkleTree::<M>::new(
            &self.config.merkle_leaf_hash_param,
            &self.config.merkle_two_to_one_param,
            evals.clone(),
        )
        .unwrap();
        (evals, mt)
    }
    fn prove(&self, mt: &MerkleTree<M>) -> Vec<Path<M>> {
        // absorb committment
        let commitment: M::InnerDigest = mt.root();
        let mut sponge = S::new(&self.config.sponge_config);
        sponge.absorb(&commitment);
        // squeeze out queries
        let mut queries: Vec<usize> = Vec::with_capacity(self.config.num_queries);
        for _ in 0..self.config.num_queries {
            queries.push(squeeze_integer(&mut sponge,  32));
        }
        // get the openings
        let mut openings: Vec<Path<M>> = Vec::with_capacity(self.config.num_queries);
        for query in queries {
            openings.push(mt.generate_proof(query).unwrap());
        }
        openings
    }
}

// Verifier
struct DirectVerifier<F: FftField, M: MerkleConfig, S: CryptographicSponge> {
    config: DirectConfig<M, S>,
    _field: PhantomData<F>,
    _merkle_config: PhantomData<M>,
    _sponge_config: PhantomData<S>,
}
impl<F: FftField, M: MerkleConfig<Leaf = Vec<F>>, S: CryptographicSponge> DirectVerifier<F, M, S>
where
    M::InnerDigest: Absorb
{
    fn new(config: DirectConfig<M, S>) -> Self {
        DirectVerifier {
            config,
            _field: PhantomData::<F>,
            _merkle_config: PhantomData::<M>,
            _sponge_config: PhantomData::<S>,
        }
    }
    fn verify(&self, evals: Vec<Vec<F>>, mt: &MerkleTree<M>, openings: Vec<Path<M>>) -> bool {
        // derive queries from committment to validate
        let commitment: M::InnerDigest = mt.root();
        let mut sponge = S::new(&self.config.sponge_config);
        sponge.absorb(&commitment);
        let mut sponge = S::new(&self.config.sponge_config);
        sponge.absorb(&commitment);
        // squeeze out queries
        for i in 0_usize..self.config.num_queries {
            let query = squeeze_integer(&mut sponge, 32);
            // is correct query index
            if openings[i].leaf_index != query {
                return false;
            }
            // is correct path
            if openings[i]
                .verify(
                    &self.config.merkle_leaf_hash_param,
                    &self.config.merkle_two_to_one_param,
                    &commitment,
                    evals[query].clone(),
                )
                .unwrap()
                != true
            {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_std::test_rng;

    use crate::{
        crypto::{fields::Field256, fs, merkle_tree},
        direct::DirectVerifier,
    };

    use super::{DirectConfig, DirectProver};

    #[test]
    fn test_direct_ldt() {
        // config
        let mut rng = test_rng();
        let mt_config = merkle_tree::poseidon::default_config::<Field256>(&mut rng, 2);
        let fs_config = fs::poseidon::default_fs_config::<Field256>();
        let config = DirectConfig {
            degree: 22,
            num_queries: 8,
            merkle_leaf_hash_param: mt_config.0.clone(),
            merkle_two_to_one_param: mt_config.1.clone(),
            sponge_config: fs_config.clone(),
        };

        // witness
        let polynomial = DensePolynomial::<Field256>::rand(config.degree, &mut rng);

        // commit
        let prover: DirectProver<
            Field256,
            merkle_tree::poseidon::MerkleTreeParams<Field256>,
            PoseidonSponge<Field256>,
        > = DirectProver::new(config);
        let (evals, mt) = prover.commit(polynomial);

        // prove
        let openings = prover.prove(&mt);

        // verify
        // TODO implement clone for DirectConfig
        let config2: DirectConfig<merkle_tree::poseidon::MerkleTreeParams<Field256>, PoseidonSponge<Field256>> = DirectConfig {
            degree: 22,
            num_queries: 8,
            merkle_leaf_hash_param: mt_config.0,
            merkle_two_to_one_param: mt_config.1,
            sponge_config: fs_config,
        };
        let verifier: DirectVerifier<
            Field256,
            merkle_tree::poseidon::MerkleTreeParams<Field256>,
            PoseidonSponge<Field256>,
        > = DirectVerifier::new(config2);
        assert_eq!(verifier.verify(evals, &mt, openings), true);
    }
}
