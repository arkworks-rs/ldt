use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, TwoToOneParam},
    sponge::CryptographicSponge,
};

pub struct DirectConfig<M: MerkleConfig, S: CryptographicSponge> {
    pub degree: usize,
    pub num_challenges: usize,
    pub merkle_leaf_hash_param: LeafParam<M>,
    pub merkle_two_to_one_param: TwoToOneParam<M>,
    pub sponge_config: S::Config,
}

impl<M: MerkleConfig, S: CryptographicSponge> DirectConfig<M, S> {
    pub fn new(
        degree: usize,
        num_challenges: usize,
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        sponge_config: S::Config,
    ) -> Self {
        DirectConfig {
            degree,
            num_challenges,
            merkle_leaf_hash_param,
            merkle_two_to_one_param,
            sponge_config,
        }
    }
}

impl<M: MerkleConfig, S: CryptographicSponge> Clone for DirectConfig<M, S>
where
    LeafParam<M>: Clone,
    TwoToOneParam<M>: Clone,
    S::Config: Clone,
{
    fn clone(&self) -> Self {
        DirectConfig {
            degree: self.degree,
            num_challenges: self.num_challenges,
            merkle_leaf_hash_param: self.merkle_leaf_hash_param.clone(),
            merkle_two_to_one_param: self.merkle_two_to_one_param.clone(),
            sponge_config: self.sponge_config.clone(),
        }
    }
}
