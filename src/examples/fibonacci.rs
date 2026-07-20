//! The Fibonacci AIR: one column, `t[i+2] = t[i+1] + t[i]`. Generic over the
//! base field, so it runs on either engine unchanged.

use crate::air::{Air, Boundary, TraceTable};
use crate::field::{ConstraintField, Field};

#[derive(Debug, Clone)]
pub struct FibonacciAir<B: Field> {
    /// Claimed value of the last row (public output).
    pub claimed_output: B,
}

impl<B: Field> FibonacciAir<B> {
    /// Build the AIR together with its honest trace of length `2^log_len`.
    pub fn with_trace(log_len: u32) -> (Self, TraceTable<B>) {
        let trace = fibonacci_trace::<B>(log_len);
        let claimed_output = *trace.columns[0].last().unwrap();
        (Self { claimed_output }, trace)
    }
}

/// The honest Fibonacci trace of length `2^log_len` (values mod the field).
pub fn fibonacci_trace<B: Field>(log_len: u32) -> TraceTable<B> {
    let n = 1usize << log_len;
    let mut col = Vec::with_capacity(n);
    let (mut a, mut b) = (B::ONE, B::ONE);
    for _ in 0..n {
        col.push(a);
        let next = a + b;
        a = b;
        b = next;
    }
    TraceTable::new(vec![col])
}

impl<B: Field> Air<B> for FibonacciAir<B> {
    fn num_columns(&self) -> usize {
        1
    }
    fn mask_offsets(&self) -> &'static [usize] {
        &[0, 1, 2]
    }
    fn num_transition_constraints(&self) -> usize {
        1
    }
    fn eval_transitions<F: ConstraintField<B>>(&self, mask: &[Vec<F>], out: &mut [F]) {
        out[0] = mask[2][0] - mask[1][0] - mask[0][0];
    }
    fn boundary_constraints(&self, trace_len: usize) -> Vec<Boundary<B>> {
        vec![
            Boundary { column: 0, row: 0, value: B::ONE },
            Boundary { column: 0, row: 1, value: B::ONE },
            Boundary { column: 0, row: trace_len - 1, value: self.claimed_output },
        ]
    }
    fn label(&self) -> &'static str {
        "fibonacci"
    }
    fn metal_transitions(&self) -> Option<String> {
        Some(
            r#"
inline void air_transitions(thread const uint m[3][1], thread uint* out) {
    out[0] = m31_sub(m[2][0], m31_add(m[1][0], m[0][0]));
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
        let (air, trace) = FibonacciAir::<BabyBear>::with_trace(6);
        assert!(trace.satisfies(&air));
        let mut bad = trace.clone();
        bad.columns[0][17] += BabyBear::ONE;
        assert!(!bad.satisfies(&air));
    }
}
