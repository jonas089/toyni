//! A multi-column range-check AIR: bytes `v < 256` with a running sum.
//! Generic over the base field.

use crate::air::{Air, Boundary, TraceTable};
use crate::field::{ConstraintField, Field};

pub const NUM_BITS: usize = 8;
/// Columns: value, sum, 8 bit columns.
pub const NUM_COLS: usize = 2 + NUM_BITS;

#[derive(Debug, Clone)]
pub struct RangeCheckAir<B: Field> {
    /// Claimed running sum at row `N-2` (public output).
    pub claimed_sum: B,
}

impl<B: Field> RangeCheckAir<B> {
    /// Build the AIR and an honest trace from a seeded byte stream.
    pub fn with_trace(log_len: u32, seed: u64) -> (Self, TraceTable<B>) {
        let n = 1usize << log_len;
        let mut state = seed | 1;
        let mut v = Vec::with_capacity(n);
        for _ in 0..n {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            v.push((state & 0xff) as u32);
        }
        let trace = range_check_trace::<B>(&v);
        let claimed_sum = trace.columns[1][n - 2];
        (Self { claimed_sum }, trace)
    }
}

/// Build the honest trace from byte values (`< 256`).
pub fn range_check_trace<B: Field>(values: &[u32]) -> TraceTable<B> {
    let n = values.len();
    assert!(n.is_power_of_two());
    assert!(values.iter().all(|&v| v < 256));

    let v: Vec<B> = values.iter().map(|&x| B::from_u64(x as u64)).collect();
    let mut s = Vec::with_capacity(n);
    let mut acc = B::ZERO;
    for &val in &v {
        s.push(acc);
        acc += val;
    }
    let mut cols = vec![v, s];
    for j in 0..NUM_BITS {
        cols.push(
            values
                .iter()
                .map(|&x| B::from_u64(((x >> j) & 1) as u64))
                .collect(),
        );
    }
    TraceTable::new(cols)
}

impl<B: Field> Air<B> for RangeCheckAir<B> {
    fn num_columns(&self) -> usize {
        NUM_COLS
    }
    fn mask_offsets(&self) -> &'static [usize] {
        &[0, 1]
    }
    fn num_transition_constraints(&self) -> usize {
        NUM_BITS + 2
    }
    fn eval_transitions<F: ConstraintField<B>>(&self, mask: &[Vec<F>], out: &mut [F]) {
        let cur = &mask[0];
        let next = &mask[1];
        for j in 0..NUM_BITS {
            let b = cur[2 + j];
            out[j] = b * b - b;
        }
        let mut acc = F::ZERO;
        let mut pow = F::ONE;
        let two = F::ONE + F::ONE;
        for j in 0..NUM_BITS {
            acc = acc + pow * cur[2 + j];
            pow = pow * two;
        }
        out[NUM_BITS] = cur[0] - acc;
        out[NUM_BITS + 1] = next[1] - cur[1] - cur[0];
    }
    fn boundary_constraints(&self, trace_len: usize) -> Vec<Boundary<B>> {
        vec![
            Boundary { column: 1, row: 0, value: B::ZERO },
            Boundary { column: 1, row: trace_len - 2, value: self.claimed_sum },
        ]
    }
    fn label(&self) -> &'static str {
        "range-check"
    }
    fn metal_transitions(&self) -> Option<String> {
        Some(
            r#"
inline void air_transitions(thread const uint m[2][10], thread uint* out) {
    uint acc = 0;
    uint pow = 1;
    for (uint j = 0; j < 8; j++) {
        uint b = m[0][2 + j];
        out[j] = m31_sub(m31_mul(b, b), b);
        acc = m31_add(acc, m31_mul(pow, b));
        pow = m31_add(pow, pow);
    }
    out[8] = m31_sub(m[0][0], acc);
    out[9] = m31_sub(m31_sub(m[1][1], m[0][1]), m[0][0]);
}
"#
            .to_string(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::BabyBear;

    #[test]
    fn honest_and_corrupt() {
        let (air, trace) = RangeCheckAir::<BabyBear>::with_trace(6, 0xDEAD);
        assert!(trace.satisfies(&air));
        let mut bad = trace.clone();
        bad.columns[2][9] = BabyBear::from_u64(2);
        assert!(!bad.satisfies(&air));
    }
}
