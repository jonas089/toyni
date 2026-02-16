// NTT implementation for BabyBear field
// Ported from joda

use crate::babybear::BabyBear;

#[inline]
fn bit_reverse(mut x: usize, log_n: usize) -> usize {
    let mut result = 0;
    for _ in 0..log_n {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    result
}

/// In-place Cooley-Tukey NTT
pub fn ntt(values: &mut [BabyBear], omega: BabyBear) {
    let n = values.len();
    assert!(n.is_power_of_two(), "NTT size must be power of 2");
    let log_n = n.trailing_zeros() as usize;

    for i in 0..n {
        let j = bit_reverse(i, log_n);
        if i < j {
            values.swap(i, j);
        }
    }

    let mut len = 2;
    while len <= n {
        let step = n / len;
        let w_len = omega.pow(step as u64);

        for i in (0..n).step_by(len) {
            let mut w = BabyBear::one();
            for j in 0..len / 2 {
                let u = values[i + j];
                let v = values[i + j + len / 2] * w;
                values[i + j] = u + v;
                values[i + j + len / 2] = u - v;
                w = w * w_len;
            }
        }
        len *= 2;
    }
}

/// In-place inverse NTT
pub fn intt(values: &mut [BabyBear], omega: BabyBear) {
    let n = values.len();

    let inv_omega = omega.pow(n as u64 - 1);
    ntt(values, inv_omega);

    let inv_n = BabyBear::new(n as u64).inverse();
    for v in values.iter_mut() {
        *v = *v * inv_n;
    }
}

/// Generate roots of unity domain
pub fn roots_of_unity_domain(n: usize) -> Vec<BabyBear> {
    assert!(n.is_power_of_two(), "Domain size must be power of 2");
    let log_n = n.trailing_zeros();
    let omega = BabyBear::get_root_of_unity(log_n);

    let mut domain = Vec::with_capacity(n);
    let mut cur = BabyBear::one();
    for _ in 0..n {
        domain.push(cur);
        cur = cur * omega;
    }
    domain
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_intt_roundtrip() {
        let n = 256usize;
        let omega = BabyBear::get_root_of_unity(n.trailing_zeros());

        let mut values: Vec<BabyBear> = (0..n).map(|i| BabyBear::new((i * 7 + 3) as u64)).collect();

        let original = values.clone();

        ntt(&mut values, omega);
        intt(&mut values, omega);

        for (a, b) in original.iter().zip(values.iter()) {
            assert_eq!(a.value, b.value, "NTT roundtrip failed");
        }
    }

    #[test]
    fn test_polynomial_evaluation() {
        let n = 8usize;
        let omega = BabyBear::get_root_of_unity(n.trailing_zeros());
        let domain = roots_of_unity_domain(n);

        let mut coeffs = vec![BabyBear::zero(); n];
        coeffs[0] = BabyBear::new(1);
        coeffs[1] = BabyBear::new(2);
        coeffs[2] = BabyBear::new(3);

        let mut evals = coeffs.clone();
        ntt(&mut evals, omega);

        assert_eq!(evals[0].value, 6);

        let x = domain[1];
        let expected = coeffs[0] + coeffs[1] * x + coeffs[2] * x * x;
        assert_eq!(evals[1].value, expected.value);
    }

    #[test]
    fn test_roots_of_unity() {
        let n = 16usize;
        let domain = roots_of_unity_domain(n);

        assert_eq!(domain[0].value, 1);

        let omega = domain[1];
        let result = omega.pow(n as u64);
        assert_eq!(result.value, 1);

        for i in 0..n {
            for j in i + 1..n {
                assert_ne!(
                    domain[i].value, domain[j].value,
                    "domain[{}] = {} equals domain[{}] = {}",
                    i, domain[i].value, j, domain[j].value
                );
            }
        }
    }
}
