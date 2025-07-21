use sha2::{Digest, Sha256};

pub mod math;
pub mod merkle;
pub mod program;
pub mod prover;
pub mod verifier;

pub fn digest_sha2(data: &[u8]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(data);
    hasher.finalize().into()
}
