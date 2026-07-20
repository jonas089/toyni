//! The circle FFT (Algorithm 1, Appendix D of the paper).
//!
//! With the twin-coset ordering `P(i) = Q·G^i` (first half) /
//! `P(i) = J(P(i - half))` (second half), the circle FFT is an ordinary
//! radix-2 butterfly network. The first layer decomposes along the involution
//! `J` with y-coordinate twiddles; every later layer decomposes along the
//! squaring map `π(x) = 2x^2 - 1` with x-coordinate twiddles.
//!
//! Interpolation (values → coefficients) uses butterflies
//! `(a, b) -> (a + b, (a - b)·t^{-1})` and a final global scale by `2^{-n}`;
//! evaluation is the exact inverse network with butterflies
//! `(a, b) -> (a + t·b, a - t·b)`.
//!
//! Coefficients come out indexed so that position `i` (bits `i_{n-1}..i_0`)
//! corresponds to the FFT-basis polynomial
//! `y^{i_{n-1}} · v_1(x)^{i_{n-2}} · ... · v_{n-1}(x)^{i_0}`.
//! Consequently a coefficient vector of order `n` embeds into order `m > n`
//! by the index map `i -> i << (m - n)` (the low-degree extension embedding).

use crate::engine::circle::geometry::CircleDomain;
use crate::field::{batch_inverse, Field, M31, QM31};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Minimum slice length before we bother with rayon.
const PAR_THRESHOLD: usize = 1 << 13;

/// Forward twiddles of a domain: `twiddles[0]` are the y-coordinates of the
/// half coset (size `N/2`); `twiddles[1 + l]` are the x-coordinates of
/// `π^l(P(k))` for `k < N / 2^(l+2)`.
#[derive(Debug, Clone)]
pub struct Twiddles {
    pub log_size: u32,
    pub forward: Vec<Vec<M31>>,
    pub inverse: Vec<Vec<M31>>,
}

impl Twiddles {
    pub fn new(domain: &CircleDomain) -> Self {
        let n = domain.log_size();
        let half = domain.half_size();
        let mut forward: Vec<Vec<M31>> = Vec::with_capacity(n as usize);

        // Layer 0: y-coordinates of the half coset.
        let half_points = domain.half_coset.points();
        forward.push(half_points.iter().map(|p| p.y).collect());

        // Layer 1 + l: x-coordinates of the l-times squared half coset.
        let mut pts = half_points;
        let mut count = half / 2;
        for _ in 1..n {
            forward.push(pts[..count].iter().map(|p| p.x).collect());
            {
                let slice = &mut pts[..count];
                #[cfg(feature = "parallel")]
                if slice.len() >= PAR_THRESHOLD {
                    slice.par_iter_mut().for_each(|p| *p = p.double());
                } else {
                    slice.iter_mut().for_each(|p| *p = p.double());
                }
                #[cfg(not(feature = "parallel"))]
                slice.iter_mut().for_each(|p| *p = p.double());
            }
            count /= 2;
        }

        let inverse = forward.iter().map(|l| par_batch_inverse(l.as_slice())).collect();

        Self { log_size: n, forward, inverse }
    }
}

/// Batch inversion, chunked in parallel for large inputs.
fn par_batch_inverse(values: &[M31]) -> Vec<M31> {
    #[cfg(feature = "parallel")]
    if values.len() >= PAR_THRESHOLD {
        return values
            .par_chunks(1 << 12)
            .flat_map_iter(batch_inverse)
            .collect();
    }
    batch_inverse(values)
}

/// Value types the CFFT operates on: anything M31-linear.
pub trait CfftValue: Field + Send + Sync {
    fn mul_m31(self, s: M31) -> Self;
}

impl CfftValue for M31 {
    #[inline(always)]
    fn mul_m31(self, s: M31) -> Self {
        self * s
    }
}

impl CfftValue for QM31 {
    #[inline(always)]
    fn mul_m31(self, s: M31) -> Self {
        self.mul_base(s)
    }
}

fn for_each_block<V: Send>(
    data: &mut [V],
    block_size: usize,
    f: impl Fn(&mut [V]) + Sync + Send,
) {
    #[cfg(feature = "parallel")]
    if data.len() >= PAR_THRESHOLD && data.len() / block_size >= 2 {
        data.par_chunks_mut(block_size).for_each(f);
        return;
    }
    data.chunks_mut(block_size).for_each(f);
}

/// In-place circle FFT interpolation: domain-ordered values → FFT-basis
/// coefficients (see module docs for the coefficient order).
pub fn interpolate<V: CfftValue>(values: &mut [V], twiddles: &Twiddles) {
    let n = twiddles.log_size;
    assert_eq!(values.len(), 1 << n, "size mismatch");
    let size = values.len();

    // Layer 0 (J-fold, y-twiddles): pairs (k, k + N/2), shared across no blocks.
    let ty = &twiddles.inverse[0];
    let (lo, hi) = values.split_at_mut(size / 2);
    #[cfg(feature = "parallel")]
    if size >= PAR_THRESHOLD {
        lo.par_iter_mut()
            .zip(hi.par_iter_mut())
            .zip(ty.par_iter())
            .for_each(|((a, b), t)| butterfly_inv(a, b, *t));
    } else {
        for ((a, b), t) in lo.iter_mut().zip(hi.iter_mut()).zip(ty.iter()) {
            butterfly_inv(a, b, *t);
        }
    }
    #[cfg(not(feature = "parallel"))]
    for ((a, b), t) in lo.iter_mut().zip(hi.iter_mut()).zip(ty.iter()) {
        butterfly_inv(a, b, *t);
    }

    // Layers 1..n (π-folds, x-twiddles): step = N / 2^(l+1), blocks of 2·step.
    let mut step = size / 4;
    for l in 1..n as usize {
        let tw = &twiddles.inverse[l];
        for_each_block(values, 2 * step, |chunk| {
            let (lo, hi) = chunk.split_at_mut(step);
            for ((a, b), t) in lo.iter_mut().zip(hi.iter_mut()).zip(tw.iter()) {
                butterfly_inv(a, b, *t);
            }
        });
        step /= 2;
    }

    // Deferred scaling: each butterfly omitted a factor 1/2.
    let scale = M31::TWO.pow(n as u64).inverse();
    scale_slice(values, scale);
}

/// In-place circle FFT evaluation: FFT-basis coefficients → domain-ordered
/// values. Exact inverse of [`interpolate`].
pub fn evaluate<V: CfftValue>(coeffs: &mut [V], twiddles: &Twiddles) {
    let n = twiddles.log_size;
    assert_eq!(coeffs.len(), 1 << n, "size mismatch");
    let size = coeffs.len();

    // Reverse layer order: π-folds from smallest to largest ...
    let mut step = 1;
    for l in (1..n as usize).rev() {
        let tw = &twiddles.forward[l];
        for_each_block(coeffs, 2 * step, |chunk| {
            let (lo, hi) = chunk.split_at_mut(step);
            for ((a, b), t) in lo.iter_mut().zip(hi.iter_mut()).zip(tw.iter()) {
                butterfly_fwd(a, b, *t);
            }
        });
        step *= 2;
    }

    // ... then the J-fold with y-twiddles.
    let ty = &twiddles.forward[0];
    let (lo, hi) = coeffs.split_at_mut(size / 2);
    #[cfg(feature = "parallel")]
    if size >= PAR_THRESHOLD {
        lo.par_iter_mut()
            .zip(hi.par_iter_mut())
            .zip(ty.par_iter())
            .for_each(|((a, b), t)| butterfly_fwd(a, b, *t));
    } else {
        for ((a, b), t) in lo.iter_mut().zip(hi.iter_mut()).zip(ty.iter()) {
            butterfly_fwd(a, b, *t);
        }
    }
    #[cfg(not(feature = "parallel"))]
    for ((a, b), t) in lo.iter_mut().zip(hi.iter_mut()).zip(ty.iter()) {
        butterfly_fwd(a, b, *t);
    }
}

/// `(a, b) -> (a + b, (a - b)·t)` where `t` is an inverse twiddle.
#[inline(always)]
fn butterfly_inv<V: CfftValue>(a: &mut V, b: &mut V, t_inv: M31) {
    let sum = *a + *b;
    let diff = (*a - *b).mul_m31(t_inv);
    *a = sum;
    *b = diff;
}

/// `(a, b) -> (a + t·b, a - t·b)`.
#[inline(always)]
fn butterfly_fwd<V: CfftValue>(a: &mut V, b: &mut V, t: M31) {
    let tb = b.mul_m31(t);
    let (x, y) = (*a + tb, *a - tb);
    *a = x;
    *b = y;
}

fn scale_slice<V: CfftValue>(values: &mut [V], s: M31) {
    #[cfg(feature = "parallel")]
    if values.len() >= PAR_THRESHOLD {
        values.par_iter_mut().for_each(|v| *v = v.mul_m31(s));
        return;
    }
    for v in values.iter_mut() {
        *v = v.mul_m31(s);
    }
}

/// Embed a coefficient vector of order `log_from` into order `log_to ≥ log_from`
/// (zero-pad in the FFT basis): index map `i -> i << (log_to - log_from)`.
pub fn embed_coeffs<V: CfftValue>(coeffs: &[V], log_from: u32, log_to: u32) -> Vec<V> {
    assert_eq!(coeffs.len(), 1 << log_from);
    assert!(log_to >= log_from);
    let shift = (log_to - log_from) as usize;
    let mut out = vec![V::ZERO; 1 << log_to];
    for (i, &c) in coeffs.iter().enumerate() {
        out[i << shift] = c;
    }
    out
}

/// Evaluate an FFT-basis coefficient vector at an arbitrary circle point over
/// the challenge field. `O(N)` multiplications via the tensor-product
/// structure of the basis.
pub fn eval_at_point<V: CfftValue>(coeffs: &[V], point: crate::engine::circle::geometry::CirclePoint<QM31>) -> QM31
where
    QM31: From<V>,
{
    let n = coeffs.len().trailing_zeros();
    assert_eq!(coeffs.len(), 1 << n);

    // Basis values at `point`, position-ordered: build the tensor product so
    // that v_{n-1}(x) is bit 0, ..., v_1(x) is bit n-2, y is bit n-1.
    let mut factors: Vec<QM31> = Vec::with_capacity(n as usize);
    // v_l(x) for l = n-1 down to 1: v_l(x) applies x -> 2x^2 - 1 (l-1) times.
    // Compute iteratively: v_1 = x, v_{l+1} = π_x(v_l applied point)...
    // We need them individually; build v_1..v_{n-1} then push in reverse.
    let two = QM31::from_base(M31::TWO);
    let mut v = point.x; // v_1(x)
    let mut vs: Vec<QM31> = Vec::with_capacity((n as usize).saturating_sub(1));
    for _ in 1..n {
        vs.push(v);
        v = two * v * v - QM31::ONE;
    }
    for &vl in vs.iter().rev() {
        factors.push(vl);
    }
    factors.push(point.y);

    let mut basis: Vec<QM31> = vec![QM31::ONE];
    for f in factors {
        let mut next = Vec::with_capacity(basis.len() * 2);
        next.extend_from_slice(&basis);
        for b in &basis {
            next.push(*b * f);
        }
        basis = next;
    }

    let mut acc = QM31::ZERO;
    for (c, b) in coeffs.iter().zip(basis.iter()) {
        acc += *b * QM31::from(*c);
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::circle::geometry::{vanishing_eval, CircleDomain, CirclePoint};
    use rand::SeedableRng;

    fn rand_values(n: usize, seed: u64) -> Vec<M31> {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        (0..n).map(|_| M31::random(&mut rng)).collect()
    }

    #[test]
    fn interpolate_evaluate_roundtrip() {
        for log in [1u32, 2, 3, 6, 10] {
            let domain = CircleDomain::standard(log);
            let tw = Twiddles::new(&domain);
            let values = rand_values(domain.size(), log as u64);
            let mut work = values.clone();
            interpolate(&mut work, &tw);
            evaluate(&mut work, &tw);
            assert_eq!(work, values, "roundtrip failed at log={log}");
        }
    }

    #[test]
    fn qm31_roundtrip() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(9);
        let domain = CircleDomain::standard(7);
        let tw = Twiddles::new(&domain);
        let values: Vec<QM31> = (0..domain.size()).map(|_| QM31::random(&mut rng)).collect();
        let mut work = values.clone();
        interpolate(&mut work, &tw);
        evaluate(&mut work, &tw);
        assert_eq!(work, values);
    }

    /// The coefficients must match the claimed basis: check the interpolant
    /// reproduces the values via direct (slow) basis evaluation at each point.
    #[test]
    fn coefficients_match_basis_order() {
        let log = 4u32;
        let domain = CircleDomain::standard(log);
        let tw = Twiddles::new(&domain);
        let values = rand_values(domain.size(), 42);
        let mut coeffs = values.clone();
        interpolate(&mut coeffs, &tw);

        for (i, p) in domain.points().iter().enumerate() {
            let got = eval_at_point(&coeffs, p.into_qm31());
            assert_eq!(got, QM31::from(values[i]), "basis mismatch at point {i}");
        }
    }

    /// LDE: interpolate on a small domain, embed the coefficients, evaluate on
    /// a larger one; spot-check against direct basis evaluation.
    #[test]
    fn low_degree_extension_by_embedding() {
        let log_small = 4u32;
        let log_big = 7u32;
        let small = CircleDomain::standard(log_small);
        let big = CircleDomain::standard(log_big);
        let tw_small = Twiddles::new(&small);
        let tw_big = Twiddles::new(&big);

        let values = rand_values(small.size(), 77);
        let mut coeffs = values.clone();
        interpolate(&mut coeffs, &tw_small);

        let mut lde = embed_coeffs(&coeffs, log_small, log_big);
        evaluate(&mut lde, &tw_big);

        for (i, p) in big.points().iter().enumerate().step_by(11) {
            let direct = eval_at_point(&coeffs, p.into_qm31());
            assert_eq!(QM31::from(lde[i]), direct, "LDE mismatch at {i}");
        }
    }

    /// The FFT space L'_N contains no v_n component: interpolating the values
    /// of v_n over the standard domain of size N must NOT reproduce v_n
    /// off-domain (dimension gap), but interpolating any polynomial of the
    /// form p0(x) + y p1(x) with deg < N/2 must be exact everywhere.
    #[test]
    fn dimension_gap_is_real() {
        let log = 4u32;
        let domain = CircleDomain::standard(log);
        let tw = Twiddles::new(&domain);

        // v_n values over the domain: all zero! (v_n vanishes on its coset.)
        // So its interpolant is 0, which differs from v_n off-domain.
        let vn_vals: Vec<M31> = domain.points().iter().map(|p| vanishing_eval(log, p.x)).collect();
        assert!(vn_vals.iter().all(|v| v.is_zero()));

        // A low-degree y-linear function is interpolated exactly.
        let f = |p: CirclePoint<M31>| p.x * p.x + p.y * (p.x + M31::new(3)) + M31::new(5);
        let vals: Vec<M31> = domain.points().iter().map(|p| f(*p)).collect();
        let mut coeffs = vals.clone();
        interpolate(&mut coeffs, &tw);
        let other = CircleDomain::standard(log + 2);
        for p in other.points().iter().take(8) {
            let want = QM31::from(f(*p));
            assert_eq!(eval_at_point(&coeffs, p.into_qm31()), want);
        }
    }
}
