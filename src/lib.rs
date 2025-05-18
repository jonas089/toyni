//! Toyni Stark - A Zero-Knowledge Virtual Machine Implementation
//!
//! This crate provides a toy implementation of a Stark proving system that can be extended
//! for use in an experimental ZKVM or other proving context like circuit arithmetic.
//!
//! # Warning
//! This project is not ready for production and has not been audited.
//! Use at own risk.
//!
//! # Modules
//!
//! * `math` - Mathematical utilities for polynomial operations and FRI protocol
//! * `vm` - Virtual machine implementation with execution tracing

use sha2::{Digest, Sha256};

pub mod math;
pub mod merkle;
pub mod vm;
pub mod prover;
pub mod verifier;

pub fn digest_sha2(data: &[u8]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(data);
    hasher.finalize().into()
}
