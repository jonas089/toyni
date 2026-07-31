//! # toyni
//!
//! A modular STARK library with two pluggable proving engines over one
//! generic [`Air`] frontend:
//!
//! - the **classical** roots-of-unity engine over BabyBear
//!   ([`ClassicalEngine`]) — NTT, multiplicative-coset domains, 2-adic FRI;
//! - the **circle** STARK engine over M31 ([`CircleEngine`], `circle`
//!   feature) — circle FFT, twin-coset domains, circle FRI, optional Metal
//!   GPU backend.
//!
//! The same AIR definition, example programs, Merkle commitments and
//! Fiat-Shamir transcript run on either engine; only the field/geometry/FFT/
//! FRI differ, behind the [`Engine`] trait.
//!
//! ```
//! use toyni::{examples::FibonacciAir, ClassicalEngine, ProofOptions,
//!             StarkProver, StarkVerifier};
//! use toyni::field::BabyBear;
//!
//! let (air, trace) = FibonacciAir::<BabyBear>::with_trace(8);
//! let options = ProofOptions::default();
//! let proof = StarkProver::<ClassicalEngine, _>::new(&air, options).prove(&trace);
//! assert!(StarkVerifier::<ClassicalEngine, _>::new(&air, options).verify(&proof));
//! ```

// `x = x * y` is idiomatic in the vendored NTT / circle-FFT / field math;
// `a / b := a * b.inverse()` trips the Div-impl lint. Both are stylistic.
#![allow(clippy::assign_op_pattern, clippy::suspicious_arithmetic_impl)]

use sha2::{Digest, Sha256};

pub mod air;
pub mod engine;
pub mod examples;
pub mod field;
pub mod mask;
pub mod merkle;
pub mod polynomial;
pub mod proof;
pub mod prover;
pub mod transcript;
pub mod verifier;

pub use air::{Air, Boundary, TraceTable};
pub use engine::classical::ClassicalEngine;
pub use engine::{Engine, Mode};
pub use field::{BabyBear, BabyBearExt, Field};
pub use proof::{ProofOptions, StarkProof};
pub use prover::StarkProver;
pub use verifier::StarkVerifier;

#[cfg(feature = "circle")]
pub use engine::circle::CircleEngine;
#[cfg(feature = "circle")]
pub use field::{CM31, M31, QM31};

/// SHA-256 of a byte slice (kept for downstream compatibility).
pub fn digest_sha2(data: &[u8]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(data);
    hasher.finalize().into()
}

/// True if the circle engine's Metal GPU backend is compiled in and a device
/// is available.
pub fn metal_available() -> bool {
    #[cfg(all(feature = "circle", feature = "metal"))]
    {
        engine::circle::metal::available()
    }
    #[cfg(not(all(feature = "circle", feature = "metal")))]
    {
        false
    }
}
