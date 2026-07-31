//! A Poseidon-style permutation-chain AIR - the heavy benchmark workload.
//! Generic over the base field. Each row applies one x^5 S-box round with
//! intermediate power columns keeping every constraint degree <= 2.

use crate::air::{Air, Boundary, TraceTable};
use crate::field::{ConstraintField, Field};

pub const STATE_WIDTH: usize = 8;
/// state, squares, fourth powers.
pub const NUM_COLS: usize = 3 * STATE_WIDTH;

/// Circulant mixing coefficients (first row of `M`); `M_ij = MIX[(j - i) mod 8]`.
const MIX: [u64; STATE_WIDTH] = [2, 1, 1, 3, 1, 1, 1, 1];
/// Per-lane round constants.
const RC: [u64; STATE_WIDTH] = [
    0x1a2b3c4, 0x5d6e7f8, 0x9abcdef, 0x1234567, 0x89abcde, 0xf013579, 0x2468ace, 0x7fffffe,
];

#[derive(Debug, Clone)]
pub struct HashChainAir<B: Field> {
    pub input: [B; STATE_WIDTH],
    /// Claimed state at row `N-2` (public output).
    pub output: [B; STATE_WIDTH],
}

/// One permutation round applied to a state.
pub fn round<B: Field>(state: &[B; STATE_WIDTH]) -> [B; STATE_WIDTH] {
    let sbox: Vec<B> = state
        .iter()
        .map(|&s| {
            let q = s * s;
            q * q * s
        })
        .collect();
    let mut out = [B::ZERO; STATE_WIDTH];
    for (i, o) in out.iter_mut().enumerate() {
        let mut acc = B::from_u64(RC[i]);
        for (j, &sb) in sbox.iter().enumerate() {
            acc += B::from_u64(MIX[(j + STATE_WIDTH - i) % STATE_WIDTH]) * sb;
        }
        *o = acc;
    }
    out
}

impl<B: Field> HashChainAir<B> {
    /// Build the AIR and its honest trace: `2^log_len` rounds from a seed state.
    pub fn with_trace(log_len: u32, seed: u64) -> (Self, TraceTable<B>) {
        let n = 1usize << log_len;
        let mut input = [B::ZERO; STATE_WIDTH];
        for (i, v) in input.iter_mut().enumerate() {
            *v = B::from_u64(seed.wrapping_mul(0x9e3779b97f4a7c15).wrapping_add(i as u64));
        }

        let mut states = Vec::with_capacity(n);
        let mut s = input;
        for _ in 0..n {
            states.push(s);
            s = round(&s);
        }

        let mut cols: Vec<Vec<B>> = (0..NUM_COLS).map(|_| Vec::with_capacity(n)).collect();
        for st in &states {
            for i in 0..STATE_WIDTH {
                let q = st[i] * st[i];
                cols[i].push(st[i]);
                cols[STATE_WIDTH + i].push(q);
                cols[2 * STATE_WIDTH + i].push(q * q);
            }
        }
        let output = states[n - 2];
        (Self { input, output }, TraceTable::new(cols))
    }
}

impl<B: Field> Air<B> for HashChainAir<B> {
    fn num_columns(&self) -> usize {
        NUM_COLS
    }
    fn mask_offsets(&self) -> &'static [usize] {
        &[0, 1]
    }
    fn num_transition_constraints(&self) -> usize {
        3 * STATE_WIDTH
    }
    fn eval_transitions<F: ConstraintField<B>>(&self, mask: &[Vec<F>], out: &mut [F]) {
        let cur = &mask[0];
        let next = &mask[1];
        for i in 0..STATE_WIDTH {
            let s = cur[i];
            let q = cur[STATE_WIDTH + i];
            let f = cur[2 * STATE_WIDTH + i];
            out[i] = q - s * s;
            out[STATE_WIDTH + i] = f - q * q;
            let mut acc: F = F::from_u64(RC[i]);
            for j in 0..STATE_WIDTH {
                let m: F = F::from_u64(MIX[(j + STATE_WIDTH - i) % STATE_WIDTH]);
                acc = acc + m * cur[2 * STATE_WIDTH + j] * cur[j];
            }
            out[2 * STATE_WIDTH + i] = next[i] - acc;
        }
    }
    fn boundary_constraints(&self, trace_len: usize) -> Vec<Boundary<B>> {
        let mut out = Vec::with_capacity(2 * STATE_WIDTH);
        for i in 0..STATE_WIDTH {
            out.push(Boundary { column: i, row: 0, value: self.input[i] });
        }
        for i in 0..STATE_WIDTH {
            out.push(Boundary { column: i, row: trace_len - 2, value: self.output[i] });
        }
        out
    }
    fn label(&self) -> &'static str {
        "hash-chain"
    }
    fn metal_transitions(&self) -> Option<String> {
        let mix = MIX.map(|v| v.to_string()).join(", ");
        let rc = RC.map(|v| format!("{v}u")).join(", ");
        Some(format!(
            r#"
constant uint HC_MIX[8] = {{ {mix} }};
constant uint HC_RC[8] = {{ {rc} }};

inline void air_transitions(thread const uint m[2][24], thread uint* out) {{
    for (uint i = 0; i < 8; i++) {{
        uint s = m[0][i];
        uint q = m[0][8 + i];
        uint f = m[0][16 + i];
        out[i] = m31_sub(q, m31_mul(s, s));
        out[8 + i] = m31_sub(f, m31_mul(q, q));
        uint acc = HC_RC[i];
        for (uint j = 0; j < 8; j++) {{
            uint mij = HC_MIX[(j + 8 - i) % 8];
            acc = m31_add(acc, m31_mul(mij, m31_mul(m[0][16 + j], m[0][j])));
        }}
        out[16 + i] = m31_sub(m[1][i], acc);
    }}
}}
"#
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::BabyBear;

    #[test]
    fn honest_and_corrupt() {
        let (air, trace) = HashChainAir::<BabyBear>::with_trace(5, 42);
        assert!(trace.satisfies(&air));
        let mut bad = trace.clone();
        bad.columns[3][7] += BabyBear::ONE;
        assert!(!bad.satisfies(&air));
    }
}
