#![forbid(unsafe_code)]

use ark_ff::PrimeField;
use std::marker::PhantomData;

pub mod prover;
pub mod verifier;
#[cfg(test)]
mod test;

pub struct FRI<F: PrimeField>{
    _protocol: PhantomData<F>
}
