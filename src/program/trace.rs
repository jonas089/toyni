use ark_bls12_381::Fr;
use ark_ff::Field;

use crate::math::polynomial::Polynomial;
pub type ProgramVariable = String;

#[derive(Clone)]
pub struct ExecutionTrace {
    pub trace: Vec<Vec<Fr>>,
}

impl ExecutionTrace {
    pub fn new() -> Self {
        Self { trace: Vec::new() }
    }

    /// Treats input as a column: inserts variable values over time.
    pub fn insert_column(&mut self, column: Vec<Fr>) {
        if self.trace.is_empty() {
            self.trace = column.into_iter().map(|val| vec![val]).collect();
        } else {
            assert_eq!(self.trace.len(), column.len(), "Column length mismatch");
            for (row, val) in self.trace.iter_mut().zip(column.into_iter()) {
                row.push(val);
            }
        }
    }

    pub fn get_column(&self, index: usize) -> Vec<Fr> {
        self.trace.iter().map(|row| row[index]).collect()
    }

    pub fn interpolate_column(&self, domain: &[Fr], column_idx: usize) -> Polynomial {
        assert_eq!(
            domain.len(),
            self.trace.len(),
            "Domain length must match trace height"
        );

        let xs = domain.to_vec();
        let ys: Vec<Fr> = self.trace.iter().map(|row| row[column_idx]).collect();

        let mut poly = Polynomial::zero();

        for (i, (xi, yi)) in xs.iter().zip(ys.iter()).enumerate() {
            let mut numerator = Polynomial::new(vec![Fr::ONE]);
            let mut denominator = Fr::ONE;

            for (j, xj) in xs.iter().enumerate() {
                if i != j {
                    numerator = numerator.mul(&Polynomial::new(vec![-*xj, Fr::ONE])); // (x - xj)
                    denominator *= *xi - *xj;
                }
            }

            let li = numerator.scale(denominator.inverse().unwrap());
            poly = poly.add(&li.scale(*yi));
        }

        poly
    }
}
