use ark_crypto_primitives::merkle_tree::Config as MerkleConfig;

pub struct SingleClaim<M>
where
    M: MerkleConfig,
{
    commitment_digest: M::InnerDigest,
}

impl<M> SingleClaim<M>
where
    M: MerkleConfig,
{
    pub fn new(commitment_digest: M::InnerDigest) -> Self {
        Self { commitment_digest }
    }
    pub fn commitment_digest(&self) -> M::InnerDigest {
        self.commitment_digest.clone()
    }
}
