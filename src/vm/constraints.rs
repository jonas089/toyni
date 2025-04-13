//! Constraint system for STARK proofs.
//!
//! Defines and evaluates constraints over execution traces, including transition
//! constraints between consecutive rows and boundary constraints at specific rows.

use ark_bls12_381::Fr;
use ark_ff::{AdditiveGroup, Zero};
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
use std::collections::HashMap;

use crate::math::polynomial::Polynomial as ToyniPolynomial;
use crate::vm::trace::{ExecutionTrace, ProgramVariable};

/// Type alias for transition constraint evaluation function
type TransitionEvaluator =
    Box<dyn Fn(&HashMap<ProgramVariable, u64>, &HashMap<ProgramVariable, u64>) -> Fr>;

/// Type alias for boundary constraint evaluation function
type BoundaryEvaluator = Box<dyn Fn(&HashMap<ProgramVariable, u64>) -> Fr>;

/// Constraint between consecutive execution trace rows.
pub struct TransitionConstraint {
    /// Constraint name for debugging
    pub name: String,
    /// Variables used in constraint
    pub variables: Vec<ProgramVariable>,
    /// Function evaluating constraint
    pub evaluate: TransitionEvaluator,
}

/// Constraint at specific execution trace row.
pub struct BoundaryConstraint {
    /// Constraint name for debugging
    pub name: String,
    /// Row where constraint must hold
    pub row: u64,
    /// Variables used in constraint
    pub variables: Vec<ProgramVariable>,
    /// Function evaluating constraint
    pub evaluate: BoundaryEvaluator,
}

/// System holding all program constraints.
#[derive(Default)]
pub struct ConstraintSystem {
    /// Constraints between consecutive rows
    pub transition_constraints: Vec<TransitionConstraint>,
    /// Constraints at specific rows
    pub boundary_constraints: Vec<BoundaryConstraint>,
}

impl ConstraintSystem {
    /// Adds transition constraint to system.
    #[allow(clippy::type_complexity)]
    pub fn add_transition_constraint(
        &mut self,
        name: String,
        variables: Vec<ProgramVariable>,
        evaluate: TransitionEvaluator,
    ) {
        self.transition_constraints.push(TransitionConstraint {
            name,
            variables,
            evaluate,
        });
    }

    /// Adds boundary constraint to system.
    #[allow(clippy::type_complexity)]
    pub fn add_boundary_constraint(
        &mut self,
        name: String,
        row: u64,
        variables: Vec<ProgramVariable>,
        evaluate: BoundaryEvaluator,
    ) {
        self.boundary_constraints.push(BoundaryConstraint {
            name,
            row,
            variables,
            evaluate,
        });
    }

    /// Evaluates all constraints on trace.
    pub fn evaluate(&self, trace: &ExecutionTrace) -> Vec<Fr> {
        let mut evaluations = Vec::new();

        // Evaluate transition constraints
        for i in 0..trace.height - 1 {
            let current_row = trace.get_column(i);
            let next_row = trace.get_column(i + 1);

            for constraint in &self.transition_constraints {
                let eval = (constraint.evaluate)(current_row, next_row);
                evaluations.push(eval);
            }
        }

        // Evaluate boundary constraints
        for constraint in &self.boundary_constraints {
            let row = trace.get_column(constraint.row);
            let eval = (constraint.evaluate)(row);
            evaluations.push(eval);
        }

        evaluations
    }

    /// Checks if all constraints are satisfied.
    pub fn is_satisfied(&self, trace: &ExecutionTrace) -> bool {
        self.evaluate(trace).iter().all(|&x| x == Fr::ZERO)
    }

    /// Interpolates transition constraint as polynomial.
    pub fn interpolate_transition_constraint(
        &self,
        trace: &ExecutionTrace,
        constraint: &TransitionConstraint,
    ) -> ToyniPolynomial {
        let domain = GeneralEvaluationDomain::<Fr>::new(trace.height as usize)
            .expect("Trace height must be a power of 2");

        let mut evaluations = vec![Fr::zero(); trace.height as usize];
        for i in 0..trace.height - 1 {
            let current_row = trace.get_column(i);
            let next_row = trace.get_column(i + 1);
            evaluations[i as usize] = (constraint.evaluate)(current_row, next_row);
        }

        let evals = Evaluations::from_vec_and_domain(evaluations, domain);
        ToyniPolynomial::from_dense_poly(evals.interpolate())
    }

    /// Interpolates boundary constraint as polynomial.
    pub fn interpolate_boundary_constraint(
        &self,
        trace: &ExecutionTrace,
        constraint: &BoundaryConstraint,
    ) -> ToyniPolynomial {
        let domain = GeneralEvaluationDomain::<Fr>::new(trace.height as usize)
            .expect("Trace height must be a power of 2");

        let mut evaluations = vec![Fr::zero(); trace.height as usize];
        evaluations[constraint.row as usize] =
            (constraint.evaluate)(trace.get_column(constraint.row));

        let evals = Evaluations::from_vec_and_domain(evaluations, domain);
        ToyniPolynomial::from_dense_poly(evals.interpolate())
    }

    /// Interpolates all constraints as polynomials.
    pub fn interpolate_all_constraints(&self, trace: &ExecutionTrace) -> Vec<ToyniPolynomial> {
        let mut polys = Vec::new();

        for constraint in &self.transition_constraints {
            polys.push(self.interpolate_transition_constraint(trace, constraint));
        }

        for constraint in &self.boundary_constraints {
            polys.push(self.interpolate_boundary_constraint(trace, constraint));
        }

        polys
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vm::trace::ExecutionTrace;

    fn create_test_trace() -> ExecutionTrace {
        let mut trace = ExecutionTrace::new(3, 2);
        for i in 0..3 {
            let mut column = HashMap::new();
            column.insert("x".to_string(), i);
            column.insert("y".to_string(), i * 2);
            trace.insert_column(column);
        }
        trace
    }

    #[test]
    fn test_transition_constraint() {
        let mut system = ConstraintSystem::default();

        // Add a constraint: y[n] = 2 * x[n]
        system.add_transition_constraint(
            "y_equals_2x".to_string(),
            vec!["x".to_string(), "y".to_string()],
            Box::new(|current, _| {
                let x = Fr::from(*current.get("x").unwrap());
                let y = Fr::from(*current.get("y").unwrap());
                y - Fr::from(2u64) * x
            }),
        );

        let trace = create_test_trace();
        assert!(system.is_satisfied(&trace));
    }

    #[test]
    fn test_boundary_constraint() {
        let mut system = ConstraintSystem::default();

        // Add a constraint: x[0] = 0
        system.add_boundary_constraint(
            "x_starts_at_zero".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        let trace = create_test_trace();
        assert!(system.is_satisfied(&trace));
    }

    #[test]
    fn test_unsatisfied_constraint() {
        let mut system = ConstraintSystem::default();

        // Add a constraint: x[n] = x[n-1] + 1
        system.add_transition_constraint(
            "x_increments".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_current = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - (x_current + Fr::from(1u64))
            }),
        );

        // Create a trace that doesn't satisfy the constraint
        let mut trace = ExecutionTrace::new(3, 1);
        for i in 0..3 {
            let mut column = HashMap::new();
            column.insert("x".to_string(), i * 2); // x[n] = 2n instead of n
            trace.insert_column(column);
        }

        assert!(!system.is_satisfied(&trace));
    }

    #[test]
    fn test_constraint_interpolation() {
        let mut system = ConstraintSystem::default();

        // Add a transition constraint: x[n] = x[n-1] + 1
        system.add_transition_constraint(
            "increment".to_string(),
            vec!["x".to_string()],
            Box::new(|current, next| {
                let x_current = Fr::from(*current.get("x").unwrap());
                let x_next = Fr::from(*next.get("x").unwrap());
                x_next - (x_current + Fr::from(1u64))
            }),
        );

        // Add a boundary constraint: x[0] = 0
        system.add_boundary_constraint(
            "start".to_string(),
            0,
            vec!["x".to_string()],
            Box::new(|row| Fr::from(*row.get("x").unwrap())),
        );

        // Create a trace that satisfies the constraints
        let mut trace = ExecutionTrace::new(4, 1);
        for i in 0..4 {
            let mut column = HashMap::new();
            column.insert("x".to_string(), i);
            trace.insert_column(column);
        }

        // Interpolate all constraints
        let polynomials = system.interpolate_all_constraints(&trace);

        // Verify that the interpolated polynomials evaluate to zero at the trace points
        for (i, poly) in polynomials.iter().enumerate() {
            for j in 0..trace.height {
                let eval = poly.evaluate(Fr::from(j));
                assert!(
                    eval.is_zero(),
                    "Polynomial {} should evaluate to zero at x={}",
                    i,
                    j
                );
            }
        }
    }
}
