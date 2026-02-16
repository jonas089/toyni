use crate::babybear::BabyBear;

pub fn fri_fold(evals: &[BabyBear], xs: &[BabyBear], beta: BabyBear) -> Vec<BabyBear> {
    assert!(evals.len() % 2 == 0, "Evaluations length must be even");
    let half = evals.len() / 2;
    let half_inv = BabyBear::new(2).inverse();
    let mut result = Vec::with_capacity(half);

    for i in 0..half {
        let a = evals[i]; // f(xi)
        let b = evals[i + half]; // f(-xi)
        let x = xs[i]; // xi

        // average part
        let avg = (a + b) * half_inv;
        // "odd" part, scaled and then divided by x to push onto x^2-domain
        let diff = (a - b) * half_inv;
        let folded = avg + diff * beta * x.inverse();

        result.push(folded);
    }

    result
}
