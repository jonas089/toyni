use crate::babybear::BabyBear;
use crate::digest_sha2;

/// Fiat-Shamir transcript for deriving verifier challenges deterministically.
/// Both prover and verifier build identical transcripts to get the same challenges.
pub struct FiatShamirTranscript {
    state: Vec<u8>,
}

impl FiatShamirTranscript {
    pub fn new() -> Self {
        Self {
            state: b"toyni-stark-v1".to_vec(),
        }
    }

    /// Absorb raw bytes into the transcript.
    pub fn absorb(&mut self, data: &[u8]) {
        self.state.extend_from_slice(data);
    }

    /// Absorb a BabyBear field element.
    pub fn absorb_field(&mut self, val: BabyBear) {
        self.absorb(&val.to_bytes());
    }

    /// Absorb a Merkle root (32-byte hash).
    pub fn absorb_commitment(&mut self, root: &[u8]) {
        self.absorb(root);
    }

    /// Squeeze a BabyBear challenge from the transcript.
    pub fn squeeze_challenge(&mut self) -> BabyBear {
        let hash = digest_sha2(&self.state);
        // Feed the hash back into state so subsequent squeezes differ
        self.state = hash.to_vec();
        BabyBear::from_bytes_mod_order(&hash)
    }

    /// Squeeze `count` distinct query indices in [0, max).
    pub fn squeeze_indices(&mut self, count: usize, max: usize) -> Vec<usize> {
        let mut indices = Vec::with_capacity(count);
        let mut seen = std::collections::HashSet::new();
        while indices.len() < count {
            let hash = digest_sha2(&self.state);
            self.state = hash.to_vec();
            // Use first 8 bytes as u64
            let val = u64::from_le_bytes(hash[..8].try_into().unwrap());
            let idx = (val as usize) % max;
            if seen.insert(idx) {
                indices.push(idx);
            }
        }
        indices
    }
}
