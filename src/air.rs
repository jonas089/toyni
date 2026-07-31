//! Algebraic Intermediate Representation: the constraint-system interface,
//! shared by both engines.
//!
//! An [`Air`] is written once, generic over a base field `B`, and runs on the
//! classical (BabyBear) or circle (M31) engine unchanged. Constraint
//! evaluators are generic over [`ConstraintField`], so the same code serves
//! the base-field low-degree extension (bulk work) and the out-of-domain check
//! over the challenge field.
//!
//! - **Transition constraints**: polynomial relations of total degree <= 2
//!   between trace cells on a window of consecutive rows (the *mask*),
//!   enforced on every row except the last [`Air::num_excluded_rows`].
//! - **Boundary constraints**: public assertions `column[row] = value`.

use crate::field::{ConstraintField, Field};

/// A boundary constraint: `trace[column][row] = value` (public input).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Boundary<B: Field> {
    pub column: usize,
    pub row: usize,
    pub value: B,
}

/// A constraint system over a power-of-two-length execution trace, generic
/// over the base field `B`.
pub trait Air<B: Field>: Sync {
    /// Number of trace columns.
    fn num_columns(&self) -> usize;

    /// Row offsets referenced by the transition constraints, contiguous from
    /// 0 (e.g. `[0, 1, 2]` for a second-order recurrence).
    fn mask_offsets(&self) -> &'static [usize];

    /// Number of transition constraints.
    fn num_transition_constraints(&self) -> usize;

    /// Evaluate every transition constraint on one mask window.
    /// `mask[k][c]` is column `c` at row offset `mask_offsets()[k]`.
    /// Constraints must have total degree <= 2 in the mask entries.
    fn eval_transitions<F: ConstraintField<B>>(&self, mask: &[Vec<F>], out: &mut [F]);

    /// Boundary constraints for a trace of length `trace_len`.
    fn boundary_constraints(&self, trace_len: usize) -> Vec<Boundary<B>>;

    /// Trailing rows on which transition constraints are *not* enforced. The
    /// exclusion set must have even size so both engines can build a
    /// low-degree selector for it (the circle engine needs a line through two
    /// points; the classical engine a product of two linear factors).
    fn num_excluded_rows(&self) -> usize {
        if self.mask_offsets().len() > 1 {
            2
        } else {
            0
        }
    }

    /// A short label mixed into the Fiat-Shamir transcript.
    fn label(&self) -> &'static str;

    /// Optional Metal Shading Language body of [`Air::eval_transitions`],
    /// enabling GPU constraint evaluation under the circle engine's `metal`
    /// feature. See `engine::circle::metal`.
    fn metal_transitions(&self) -> Option<String> {
        None
    }
}

/// A column-major execution trace over base field `B`.
#[derive(Debug, Clone)]
pub struct TraceTable<B: Field> {
    pub columns: Vec<Vec<B>>,
}

impl<B: Field> TraceTable<B> {
    pub fn new(columns: Vec<Vec<B>>) -> Self {
        assert!(!columns.is_empty(), "trace must have at least one column");
        let len = columns[0].len();
        assert!(len.is_power_of_two(), "trace length must be a power of two");
        assert!(columns.iter().all(|c| c.len() == len), "ragged trace");
        Self { columns }
    }

    pub fn len(&self) -> usize {
        self.columns[0].len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn log_len(&self) -> u32 {
        self.len().trailing_zeros()
    }

    pub fn num_columns(&self) -> usize {
        self.columns.len()
    }

    /// Check every constraint directly (debug/test helper - O(N*C)).
    pub fn satisfies<A: Air<B>>(&self, air: &A) -> bool {
        let n = self.len();
        let offsets = air.mask_offsets();
        let excluded = air.num_excluded_rows();
        let mut out = vec![B::ZERO; air.num_transition_constraints()];
        for row in 0..n.saturating_sub(excluded) {
            let mask: Vec<Vec<B>> = offsets
                .iter()
                .map(|&k| self.columns.iter().map(|c| c[(row + k) % n]).collect())
                .collect();
            air.eval_transitions(&mask, &mut out);
            if out.iter().any(|v| !v.is_zero()) {
                return false;
            }
        }
        air.boundary_constraints(n)
            .iter()
            .all(|b| self.columns[b.column][b.row] == b.value)
    }
}
