use crate::babybear::BabyBear;
use crate::math::polynomial::Polynomial;

pub type ProgramVariable = String;

#[derive(Clone)]
pub struct ExecutionTrace {
    pub trace: Vec<Vec<BabyBear>>,
}

impl ExecutionTrace {
    pub fn new() -> Self {
        Self { trace: Vec::new() }
    }

    /// Treats input as a column: inserts variable values over time.
    pub fn insert_column(&mut self, column: Vec<BabyBear>) {
        if self.trace.is_empty() {
            self.trace = column.into_iter().map(|val| vec![val]).collect();
        } else {
            assert_eq!(self.trace.len(), column.len(), "Column length mismatch");
            for (row, val) in self.trace.iter_mut().zip(column.into_iter()) {
                row.push(val);
            }
        }
    }

    pub fn interpolate_column(&self, domain: &[BabyBear], column_idx: usize) -> Polynomial {
        assert_eq!(
            domain.len(),
            self.trace.len(),
            "Domain length must match trace height"
        );

        let xs = domain.to_vec();
        let ys: Vec<BabyBear> = self.trace.iter().map(|row| row[column_idx]).collect();

        let mut poly = Polynomial::zero();

        for (i, (xi, yi)) in xs.iter().zip(ys.iter()).enumerate() {
            let mut numerator = Polynomial::new(vec![BabyBear::one()]);
            let mut denominator = BabyBear::one();

            for (j, xj) in xs.iter().enumerate() {
                if i != j {
                    numerator = numerator.mul(&Polynomial::new(vec![-*xj, BabyBear::one()]));
                    denominator = denominator * (*xi - *xj);
                }
            }

            let li = numerator.scale(denominator.inverse());
            poly = poly.add(&li.scale(*yi));
        }

        poly
    }
}
