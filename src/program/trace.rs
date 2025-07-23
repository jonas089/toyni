use ark_bls12_381::Fr;
use ark_ff::{AdditiveGroup, Field};
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

    pub fn get_column(&self, index: u64) -> &Vec<Fr> {
        &self.trace[index as usize]
    }

    pub fn interpolate_column(&self, domain: &[Fr], column_idx: usize) -> impl Fn(Fr) -> Fr {
        assert_eq!(
            domain.len(),
            self.trace.len(),
            "Domain length must match trace height"
        );

        let xs = domain.to_vec();
        let ys: Vec<Fr> = self.trace.iter().map(|row| row[column_idx]).collect();

        move |x: Fr| {
            let mut result = Fr::ZERO;

            for (i, (xi, yi)) in xs.iter().zip(ys.iter()).enumerate() {
                let mut li = Fr::from(1);
                for (j, xj) in xs.iter().enumerate() {
                    if i != j {
                        let num = x - xj;
                        let denom = *xi - *xj;
                        li *= num * denom.inverse().expect("division by zero in interpolation");
                    }
                }
                result += *yi * li;
            }

            result
        }
    }
}

#[test]
fn test_interpolate_column_exact_points() {
    let mut trace = ExecutionTrace::new();
    trace.insert_column(vec![Fr::from(1), Fr::from(2), Fr::from(3)]); // Adds values for column 0
    trace.insert_column(vec![Fr::from(4), Fr::from(5), Fr::from(6)]); // Adds values for column 1
    let domain = &[Fr::from(2), Fr::from(3), Fr::from(4)];
    let poly = trace.interpolate_column(domain, 0);
    assert_eq!(poly(Fr::from(2)), Fr::from(1));
    assert_eq!(poly(Fr::from(3)), Fr::from(2));
    assert_eq!(poly(Fr::from(4)), Fr::from(3));

    let domain = &[Fr::from(5), Fr::from(6), Fr::from(7)];
    let poly = &trace.interpolate_column(domain, 1);
    assert_eq!(poly(Fr::from(5)), Fr::from(4));
    assert_eq!(poly(Fr::from(6)), Fr::from(5));
    assert_eq!(poly(Fr::from(7)), Fr::from(6));
}
