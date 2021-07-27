#![cfg_attr(not(feature = "std"), no_std)]

//! A crate for low-degree tests.
#![deny(
    const_err,
    future_incompatible,
    missing_docs,
    non_shorthand_field_patterns,
    renamed_and_removed_lints,
    rust_2018_idioms,
    stable_features,
    trivial_casts,
    trivial_numeric_casts,
    unused,
    variant_size_differences,
    warnings
)]
#![forbid(unsafe_code)]

/// Direct low-degree tests
pub mod direct;

/// Domain represented as coset.
pub mod domain;
/// Implementations for FRI Protocol
pub mod fri;
