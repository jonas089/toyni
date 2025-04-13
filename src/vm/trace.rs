//! Execution trace recording for virtual machine.
//!
//! Records program execution as a matrix where columns are variables and rows are execution steps.

use std::collections::HashMap;

/// Program variable name type.
pub type ProgramVariable = String;

/// Execution trace storing program state changes.
pub struct ExecutionTrace {
    /// Number of execution steps
    pub height: u64,
    /// Number of program variables
    pub width: u64,
    /// Trace data as vector of variable-value maps
    pub trace: Vec<HashMap<ProgramVariable, u64>>,
}

impl ExecutionTrace {
    /// Creates empty trace with given dimensions.
    pub fn new(height: u64, width: u64) -> Self {
        Self {
            height,
            width,
            trace: Vec::new(),
        }
    }

    /// Adds new execution step to trace.
    pub fn insert_column(&mut self, column: HashMap<ProgramVariable, u64>) {
        assert!(column.len() == self.width as usize);
        assert!(self.trace.len() < self.height as usize);
        self.trace.push(column);
    }

    /// Gets execution step by index.
    pub fn get_column(&self, index: u64) -> &HashMap<ProgramVariable, u64> {
        &self.trace[index as usize]
    }

    /// Prints trace in tabular format.
    pub fn print_trace(&self, variables: Vec<ProgramVariable>) {
        for i in 0..self.height {
            let column = self.get_column(i);
            for var in &variables {
                print!("{} |", column.get(var).unwrap_or(&0));
            }
            println!();
        }
    }

    /// Interpolates variable value between two steps.
    pub fn interpolate(&self, variable: &ProgramVariable, step1: u64, step2: u64, t: u8) -> u64 {
        assert!(
            step1 < self.height && step2 < self.height,
            "Step indices out of bounds"
        );
        assert!(
            t <= 100,
            "Interpolation parameter must be between 0 and 100"
        );

        let col1 = self.get_column(step1);
        let col2 = self.get_column(step2);

        let val1 = *col1
            .get(variable)
            .expect("Variable not found in first step");
        let val2 = *col2
            .get(variable)
            .expect("Variable not found in second step");

        let diff = val2 - val1;
        let scaled_diff = (diff * (t as u64)) / 100;
        val1 + scaled_diff
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn generate_test_trace() -> ExecutionTrace {
        let mut execution_trace = ExecutionTrace::new(5, 5);
        for i in 0..execution_trace.height {
            let mut column = HashMap::new();
            column.insert("a".to_string(), i);
            column.insert("b".to_string(), i + 1);
            column.insert("c".to_string(), i + 2);
            column.insert("d".to_string(), i + 3);
            column.insert("e".to_string(), i + 4);
            execution_trace.insert_column(column);
        }
        execution_trace
    }

    #[test]
    fn print_test_trace() {
        let execution_trace = generate_test_trace();
        execution_trace.print_trace(vec![
            "a".to_string(),
            "b".to_string(),
            "c".to_string(),
            "d".to_string(),
            "e".to_string(),
        ]);
    }

    #[test]
    fn test_interpolation() {
        let execution_trace = generate_test_trace();

        // Test variable "a" interpolation between steps 0 and 1
        // step 0: 0, step 1: 1
        let interpolated = execution_trace.interpolate(&"a".to_string(), 0, 1, 50);
        assert_eq!(interpolated, 0); // At 50% between 0 and 1

        // Test variable "b" interpolation between steps 0 and 1
        // step 0: 1, step 1: 2
        let interpolated = execution_trace.interpolate(&"b".to_string(), 0, 1, 0);
        assert_eq!(interpolated, 1);

        // Test variable "c" interpolation between steps 0 and 1
        // step 0: 2, step 1: 3
        let interpolated = execution_trace.interpolate(&"c".to_string(), 0, 1, 100);
        assert_eq!(interpolated, 3);

        // Additional test cases
        let interpolated = execution_trace.interpolate(&"a".to_string(), 0, 1, 100);
        assert_eq!(interpolated, 1);

        let interpolated = execution_trace.interpolate(&"b".to_string(), 0, 1, 100);
        assert_eq!(interpolated, 2);
    }
}
