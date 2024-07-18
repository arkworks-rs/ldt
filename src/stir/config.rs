use ark_crypto_primitives::{
    merkle_tree::{Config as MerkleConfig, LeafParam, TwoToOneParam},
    sponge::CryptographicSponge,
};

pub struct STIRConfig<M: MerkleConfig, S: CryptographicSponge> {
    pub folding_factor: usize,
    pub num_rounds: usize,
    pub merkle_leaf_hash_param: LeafParam<M>,
    pub merkle_two_to_one_param: TwoToOneParam<M>,
    pub num_out_of_domain_samples: usize,
    pub num_proof_of_work_bits: Vec<usize>,
    pub num_repetitions: Vec<usize>,
    pub sponge_config: S::Config,
    pub starting_degree: usize,
    pub starting_rate: usize,
    pub stopping_degree: usize,
}

impl<M: MerkleConfig, S: CryptographicSponge> STIRConfig<M, S> {
    pub fn new(
        folding_factor: usize,
        num_rounds: usize,
        merkle_leaf_hash_param: LeafParam<M>,
        merkle_two_to_one_param: TwoToOneParam<M>,
        num_out_of_domain_samples: usize,
        num_proof_of_work_bits: Vec<usize>,
        num_repetitions: Vec<usize>,
        sponge_config: S::Config,
        starting_degree: usize,
        starting_rate: usize,
        stopping_degree: usize,
    ) -> Self {
        Self {
            folding_factor,
            num_rounds,
            merkle_leaf_hash_param,
            merkle_two_to_one_param,
            num_out_of_domain_samples,
            num_proof_of_work_bits,
            num_repetitions,
            sponge_config,
            starting_degree,
            starting_rate,
            stopping_degree,
        }
    }
}

impl<M: MerkleConfig, S: CryptographicSponge> Clone for STIRConfig<M, S>
where
    LeafParam<M>: Clone,
    TwoToOneParam<M>: Clone,
    S::Config: Clone,
{
    fn clone(&self) -> Self {
        Self {
            folding_factor: self.folding_factor,
            num_rounds: self.num_rounds,
            num_out_of_domain_samples: self.num_out_of_domain_samples,
            merkle_leaf_hash_param: self.merkle_leaf_hash_param.clone(),
            merkle_two_to_one_param: self.merkle_two_to_one_param.clone(),
            num_proof_of_work_bits: self.num_proof_of_work_bits.clone(),
            num_repetitions: self.num_repetitions.clone(),
            sponge_config: self.sponge_config.clone(),
            starting_degree: self.starting_degree,
            starting_rate: self.starting_rate,
            stopping_degree: self.stopping_degree,
        }
    }
}
