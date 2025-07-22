use std::collections::HashMap;
pub type ProgramVariable = String;

#[derive(Clone)]
pub struct ExecutionTrace {
    pub trace: Vec<Vec<u64>>,
}

impl ExecutionTrace {
    pub fn new() -> Self {
        Self { trace: Vec::new() }
    }

    /// Treats input as a column: inserts variable values over time.
    pub fn insert_column(&mut self, column: Vec<u64>) {
        if self.trace.is_empty() {
            self.trace = column.into_iter().map(|val| vec![val]).collect();
        } else {
            assert_eq!(self.trace.len(), column.len(), "Column length mismatch");
            for (row, val) in self.trace.iter_mut().zip(column.into_iter()) {
                row.push(val);
            }
        }
    }

    pub fn get_column(&self, index: u64) -> &Vec<u64> {
        &self.trace[index as usize]
    }

    pub fn interpolate_column(&self, domain: &[f64], column_idx: usize) -> impl Fn(f64) -> f64 {
        assert_eq!(
            domain.len(),
            self.trace.len(),
            "Domain length must match trace height"
        );

        let xs = domain.to_vec();
        let ys: Vec<f64> = self
            .trace
            .iter()
            .map(|row| row[column_idx] as f64)
            .collect();

        move |x: f64| {
            let mut result = 0.0;

            for (i, (&xi, &yi)) in xs.iter().zip(ys.iter()).enumerate() {
                let mut li = 1.0;
                for (j, &xj) in xs.iter().enumerate() {
                    if i != j {
                        li *= (x - xj) / (xi - xj);
                    }
                }
                result += yi * li;
            }

            result
        }
    }
}

#[test]
fn test_interpolate_column_exact_points() {
    let mut trace = ExecutionTrace::new();
    trace.insert_column(vec![1, 2, 3]); // Adds values for column 0
    trace.insert_column(vec![4, 5, 6]); // Adds values for column 1
    let poly = trace.interpolate_column(&[2.0, 3.0, 4.0], 0);
    assert_eq!(poly(2.0), 1.0);
    assert_eq!(poly(3.0), 2.0);
    assert_eq!(poly(4.0), 3.0);

    let poly = &trace.interpolate_column(&[5.0, 6.0, 7.0], 1);
    assert_eq!(poly(5.0), 4.0);
    assert_eq!(poly(6.0), 5.0);
    assert_eq!(poly(7.0), 6.0);
}
