//! Fiat-Shamir transcript over SHA-256, shared by both engines.
//!
//! A fixed-width sponge: the 32-byte state absorbs data as
//! `state = H(0x02 || state || len || data)` and squeezes as
//! `out = H(0x03 || state || ctr)`, ratcheting the state after every squeeze
//! so challenges drawn after a commitment is absorbed are binding on it.
//!
//! The transcript itself is field-agnostic (a byte oracle). Base-field draws
//! go through [`Transcript::draw_field`]; extension-field and geometric draws
//! (circle points) are built on top by each engine, which owns how many base
//! limbs a challenge consumes.

use sha2::{Digest, Sha256};

use crate::field::Field;
use crate::merkle::Hash;

const ABSORB_TAG: u8 = 0x02;
const SQUEEZE_TAG: u8 = 0x03;

pub struct Transcript {
    state: Hash,
    counter: u64,
}

impl Transcript {
    pub fn new(protocol_label: &[u8]) -> Self {
        let mut h = Sha256::new();
        h.update(b"toyni-stark-v2/");
        h.update(protocol_label);
        Self {
            state: h.finalize().into(),
            counter: 0,
        }
    }

    pub fn absorb(&mut self, data: &[u8]) {
        let mut h = Sha256::new();
        h.update([ABSORB_TAG]);
        h.update(self.state);
        h.update((data.len() as u64).to_le_bytes());
        h.update(data);
        self.state = h.finalize().into();
        self.counter = 0;
    }

    pub fn absorb_commitment(&mut self, root: &Hash) {
        self.absorb(root);
    }

    pub fn absorb_field<F: Field>(&mut self, v: F) {
        self.absorb(&v.to_bytes());
    }

    /// Squeeze 32 raw bytes and advance the squeeze counter (state unchanged
    /// so successive squeezes without a ratchet differ only by counter).
    pub fn squeeze_bytes(&mut self) -> Hash {
        let mut h = Sha256::new();
        h.update([SQUEEZE_TAG]);
        h.update(self.state);
        h.update(self.counter.to_le_bytes());
        let out: Hash = h.finalize().into();
        self.counter += 1;
        out
    }

    /// Fold a squeeze back into the state so later absorbs commit to it.
    pub fn ratchet(&mut self) {
        let out = self.squeeze_bytes();
        self.state = out;
        self.counter = 0;
    }

    /// Draw a uniform base-field element.
    pub fn draw_field<F: Field>(&mut self) -> F {
        let bytes = self.squeeze_bytes();
        self.ratchet();
        F::from_bytes_reduce(&bytes)
    }

    /// Draw `count` distinct indices in `[0, max)`; `max` must be a power of
    /// two (masking introduces no modular bias).
    pub fn draw_indices(&mut self, count: usize, max: usize) -> Vec<usize> {
        assert!(max.is_power_of_two());
        assert!(count <= max);
        let mask = (max - 1) as u64;
        let mut seen = std::collections::HashSet::new();
        let mut out = Vec::with_capacity(count);
        while out.len() < count {
            let bytes = self.squeeze_bytes();
            for chunk in bytes.chunks(8) {
                if out.len() >= count {
                    break;
                }
                let v = u64::from_le_bytes(chunk.try_into().unwrap()) & mask;
                let idx = v as usize;
                if seen.insert(idx) {
                    out.push(idx);
                }
            }
        }
        self.ratchet();
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::BabyBear;

    #[test]
    fn deterministic_replay() {
        let mut a = Transcript::new(b"t");
        let mut b = Transcript::new(b"t");
        a.absorb(b"hello");
        b.absorb(b"hello");
        assert_eq!(a.draw_field::<BabyBear>(), b.draw_field::<BabyBear>());
        assert_eq!(a.draw_indices(5, 64), b.draw_indices(5, 64));
    }

    #[test]
    fn diverges_on_input() {
        let mut a = Transcript::new(b"t");
        let mut b = Transcript::new(b"t");
        a.absorb(b"x");
        b.absorb(b"y");
        assert_ne!(a.draw_field::<BabyBear>(), b.draw_field::<BabyBear>());
    }
}
