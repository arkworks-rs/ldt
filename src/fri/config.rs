use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, TwoToOneParam},
    sponge::CryptographicSponge,
};

#[derive(Clone)]
pub struct FRIConfig<M: MerkleConfig, S: CryptographicSponge> {
    pub folding_factor: usize,
    pub num_queries: usize,
    pub num_rounds: usize,
    pub merkle_leaf_hash_param: LeafParam<M>,
    pub merkle_two_to_one_param: TwoToOneParam<M>,
    pub proof_of_work_bits: usize,
    pub repetitions: usize,
    pub sponge_config: S::Config,
    pub starting_degree: usize,
    pub starting_rate: usize,
}

impl<M: MerkleConfig, S: CryptographicSponge> FRIConfig<M, S> {
    pub fn new(
        folding_factor: usize,
        num_queries: usize,
        num_rounds: usize,
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        proof_of_work_bits: usize,
        repetitions: usize,
        sponge_config: S::Config,
        starting_degree: usize,
        starting_rate: usize,
    ) -> Self {
        Self {
            folding_factor,
            num_queries,
            num_rounds,
            merkle_leaf_hash_param,
            merkle_two_to_one_param,
            proof_of_work_bits,
            repetitions,
            sponge_config,
            starting_degree,
            starting_rate,
        }
    }
}
