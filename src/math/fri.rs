use crate::babybear::BabyBear;
use crate::ext::Ext;

/// FRI fold of an extension-field codeword. The evaluation points `xs` stay in
/// the base field (the domain is base; squaring base points stays base), so only
/// the codeword values and the folding challenge `beta` are extension elements.
pub fn fri_fold_ext(evals: &[Ext], xs: &[BabyBear], beta: Ext) -> Vec<Ext> {
    assert!(evals.len() % 2 == 0, "Evaluations length must be even");
    let half = evals.len() / 2;
    let half_inv = BabyBear::new(2).inverse();
    let mut result = Vec::with_capacity(half);

    for i in 0..half {
        let a = evals[i];
        let b = evals[i + half];
        let x_inv = xs[i].inverse();

        let avg = (a + b).mul_base(half_inv);
        let diff = (a - b).mul_base(half_inv);
        let folded = avg + diff * beta * Ext::from_base(x_inv);
        result.push(folded);
    }

    result
}

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
