use crate::prover::StarkProof;

pub struct StarkVerifier;

impl StarkVerifier {
    pub fn verify(&self, proof: &StarkProof) -> bool {
        true
    }
}
