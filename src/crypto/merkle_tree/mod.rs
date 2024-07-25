use core::sync::atomic::{AtomicUsize, Ordering};
use spin::{Mutex, Once};

pub mod mock;
pub mod poseidon;

use ark_std::{borrow::Borrow, marker::PhantomData, vec::Vec};

use ark_crypto_primitives::crh::CRHScheme;
use ark_serialize::CanonicalSerialize;
use ark_std::{rand::RngCore, vec};

#[derive(Debug, Default)]
pub struct HashCounter {
    counter: AtomicUsize,
}

static INIT: Once = Once::new();
static mut HASH_COUNTER: Option<Mutex<HashCounter>> = None;

impl HashCounter {
    fn get_instance() -> &'static Mutex<HashCounter> {
        unsafe {
            INIT.call_once(|| {
                HASH_COUNTER = Some(Mutex::new(HashCounter::default()));
            });
            HASH_COUNTER.as_ref().unwrap()
        }
    }

    pub(crate) fn add() -> usize {
        let counter = Self::get_instance().lock();
        counter.counter.fetch_add(1, Ordering::SeqCst)
    }

    pub fn reset() {
        let counter = Self::get_instance().lock();
        counter.counter.store(0, Ordering::SeqCst)
    }

    pub fn get() -> usize {
        let counter = Self::get_instance().lock();
        counter.counter.load(Ordering::SeqCst)
    }
}

#[derive(Debug, Default)]
pub struct LeafIdentityHasher<F>(PhantomData<F>);

impl<F: CanonicalSerialize + Send> CRHScheme for LeafIdentityHasher<F> {
    type Input = F;
    type Output = Vec<u8>;
    type Parameters = ();

    fn setup<R: RngCore>(_: &mut R) -> Result<Self::Parameters, ark_crypto_primitives::Error> {
        Ok(())
    }

    fn evaluate<T: Borrow<Self::Input>>(
        _: &Self::Parameters,
        input: T,
    ) -> Result<Self::Output, ark_crypto_primitives::Error> {
        let mut buf = vec![];
        CanonicalSerialize::serialize_compressed(input.borrow(), &mut buf)?;
        Ok(buf)
    }
}
