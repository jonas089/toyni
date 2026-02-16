// BabyBear evaluation domain with NTT-based FFT/IFFT
// Replaces Arkworks GeneralEvaluationDomain

use crate::babybear::BabyBear;
use crate::ntt_babybear;

/// Evaluation domain over BabyBear field, supporting standard and coset domains.
#[derive(Debug, Clone)]
pub struct BabyBearDomain {
    pub size: usize,
    pub log_size: u32,
    pub omega: BabyBear,
    pub shift: BabyBear, // 1 for standard domain, h for coset {h*omega^i}
    pub use_gpu: bool,
}

impl BabyBearDomain {
    /// Create a standard evaluation domain of the given size (must be power of 2).
    pub fn new(size: usize) -> Self {
        assert!(size.is_power_of_two(), "Domain size must be power of 2");
        let log_size = size.trailing_zeros();
        let omega = BabyBear::get_root_of_unity(log_size);
        Self {
            size,
            log_size,
            omega,
            shift: BabyBear::one(),
            use_gpu: false,
        }
    }

    /// Create a coset domain: {shift * omega^i} for i = 0..size-1.
    pub fn get_coset(&self, shift: BabyBear) -> Self {
        Self {
            size: self.size,
            log_size: self.log_size,
            omega: self.omega,
            shift,
            use_gpu: self.use_gpu,
        }
    }

    /// Set whether to use GPU acceleration for NTT operations.
    pub fn with_gpu(mut self, use_gpu: bool) -> Self {
        self.use_gpu = use_gpu;
        self
    }

    /// The generator of this domain (primitive root of unity).
    pub fn group_gen(&self) -> BabyBear {
        self.omega
    }

    /// The domain size.
    pub fn size(&self) -> usize {
        self.size
    }

    /// Returns all domain elements: {shift * omega^i} for i = 0..size-1.
    pub fn elements(&self) -> Vec<BabyBear> {
        let mut result = Vec::with_capacity(self.size);
        let mut cur = self.shift;
        for _ in 0..self.size {
            result.push(cur);
            cur = cur * self.omega;
        }
        result
    }

    /// Returns the coefficients of the vanishing polynomial over this domain.
    /// Standard domain: x^n - 1
    /// Coset domain {h*omega^i}: x^n - h^n
    pub fn vanishing_poly_coeffs(&self) -> Vec<BabyBear> {
        let h_n = self.shift.pow(self.size as u64);
        let mut coeffs = vec![BabyBear::zero(); self.size + 1];
        coeffs[0] = -h_n;
        coeffs[self.size] = BabyBear::one();
        coeffs
    }

    /// IFFT: recover polynomial coefficients from evaluations on this domain.
    /// For standard domain: INTT(evals)
    /// For coset domain {h*omega^i}: INTT(evals), then divide coeff[i] by h^i.
    pub fn ifft(&self, evals: &[BabyBear]) -> Vec<BabyBear> {
        assert_eq!(evals.len(), self.size, "Evaluation count must match domain size");
        let mut values = evals.to_vec();

        // NTT/INTT (GPU or CPU)
        #[cfg(feature = "cuda")]
        if self.use_gpu {
            if crate::cuda_ntt::cuda_available() {
                crate::cuda_ntt::intt_cuda(&mut values).expect("CUDA INTT failed");
                self.undo_coset_shift(&mut values);
                return values;
            }
        }

        ntt_babybear::intt(&mut values, self.omega);
        self.undo_coset_shift(&mut values);
        values
    }

    /// FFT: evaluate polynomial at all domain points.
    /// For standard domain: NTT(coeffs)
    /// For coset domain {h*omega^i}: multiply coeffs[i] by h^i, then NTT.
    pub fn fft(&self, coeffs: &[BabyBear]) -> Vec<BabyBear> {
        let mut values = coeffs.to_vec();
        values.resize(self.size, BabyBear::zero());

        self.apply_coset_shift(&mut values);

        #[cfg(feature = "cuda")]
        if self.use_gpu {
            if crate::cuda_ntt::cuda_available() {
                crate::cuda_ntt::ntt_cuda(&mut values).expect("CUDA NTT failed");
                return values;
            }
        }

        ntt_babybear::ntt(&mut values, self.omega);
        values
    }

    /// Apply coset shift: multiply values[i] by shift^i (for coset FFT).
    fn apply_coset_shift(&self, values: &mut [BabyBear]) {
        if self.shift.value != 1 {
            let mut shift_pow = BabyBear::one();
            for v in values.iter_mut() {
                *v = *v * shift_pow;
                shift_pow = shift_pow * self.shift;
            }
        }
    }

    /// Undo coset shift: divide values[i] by shift^i (for coset IFFT).
    fn undo_coset_shift(&self, values: &mut [BabyBear]) {
        if self.shift.value != 1 {
            let shift_inv = self.shift.inverse();
            let mut shift_pow = BabyBear::one();
            for v in values.iter_mut() {
                *v = *v * shift_pow;
                shift_pow = shift_pow * shift_inv;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_domain_elements() {
        let domain = BabyBearDomain::new(8);
        let elements = domain.elements();
        assert_eq!(elements.len(), 8);
        assert_eq!(elements[0].value, 1); // omega^0 = 1

        // omega^8 should = 1
        let omega = domain.group_gen();
        assert_eq!(omega.pow(8).value, 1);
    }

    #[test]
    fn test_fft_ifft_roundtrip() {
        let domain = BabyBearDomain::new(8);
        let coeffs: Vec<BabyBear> = (0..8).map(|i| BabyBear::new(i * 3 + 1)).collect();

        let evals = domain.fft(&coeffs);
        let recovered = domain.ifft(&evals);

        for (a, b) in coeffs.iter().zip(recovered.iter()) {
            assert_eq!(a.value, b.value, "FFT/IFFT roundtrip failed");
        }
    }

    #[test]
    fn test_coset_fft_ifft_roundtrip() {
        let domain = BabyBearDomain::new(8);
        let coset = domain.get_coset(BabyBear::new(7));
        let coeffs: Vec<BabyBear> = (0..8).map(|i| BabyBear::new(i * 3 + 1)).collect();

        let evals = coset.fft(&coeffs);
        let recovered = coset.ifft(&evals);

        for (a, b) in coeffs.iter().zip(recovered.iter()) {
            assert_eq!(a.value, b.value, "Coset FFT/IFFT roundtrip failed");
        }
    }

    #[test]
    fn test_coset_evaluations_correct() {
        use crate::math::polynomial::Polynomial;

        let domain = BabyBearDomain::new(8);
        let coset = domain.get_coset(BabyBear::new(7));

        // Create polynomial: 1 + 2x + 3x^2
        let coeffs = vec![BabyBear::new(1), BabyBear::new(2), BabyBear::new(3)];
        let poly = Polynomial::new(coeffs.clone());

        let evals = coset.fft(&coeffs);
        let elements = coset.elements();

        for (i, &x) in elements.iter().enumerate() {
            let expected = poly.evaluate(x);
            assert_eq!(
                evals[i].value, expected.value,
                "Coset evaluation mismatch at index {}",
                i
            );
        }
    }

    #[test]
    fn test_vanishing_polynomial() {
        let domain = BabyBearDomain::new(8);
        let vp = domain.vanishing_poly_coeffs();
        let vpoly = crate::math::polynomial::Polynomial::new(vp);

        // The vanishing polynomial should be zero at all domain elements
        for x in domain.elements() {
            let val = vpoly.evaluate(x);
            assert_eq!(val.value, 0, "Vanishing poly should be 0 at domain element");
        }
    }

    #[test]
    fn test_extended_domain_contains_original() {
        let original = BabyBearDomain::new(4);
        let extended = BabyBearDomain::new(4 * 8);

        let orig_points = original.elements();
        let ext_points = extended.elements();

        for (i, &point) in orig_points.iter().enumerate() {
            assert_eq!(point.value, ext_points[i * 8].value);
        }
    }
}
