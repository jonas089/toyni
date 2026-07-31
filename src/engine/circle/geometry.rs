//! The circle curve `C(F): x^2 + y^2 = 1` as a cyclic group, plus the coset
//! and domain machinery of the circle FFT.
//!
//! Over M31 the curve has `p + 1 = 2^31` points and is cyclic; rotations
//! `T_P(Q) = P·Q` play the role that multiplication by roots of unity plays
//! in a classical STARK. Domains are *twin-cosets* `Q·G_{n-1} ∪ Q^{-1}·G_{n-1}`
//! (Definition 2 of the paper); the unique *standard position coset* of size
//! `2^n` is the special twin-coset that is itself a coset of `G_n`.
//!
//! Points of M31 domains are indexed by their exponent with respect to a fixed
//! generator of the full group, i.e. an integer mod `2^31`, which makes coset
//! arithmetic (rotation, squaring, conjugation) pure index arithmetic.

use crate::field::{Field, M31, QM31};

/// A point on the circle curve over the field `F`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CirclePoint<F> {
    pub x: F,
    pub y: F,
}

#[allow(clippy::should_implement_trait)] // `add`/`mul` mirror the group-law
// naming of the paper; operator traits are deliberately not implemented so
// call sites read as explicit group operations.
impl<F: Field> CirclePoint<F> {
    pub fn identity() -> Self {
        Self { x: F::ONE, y: F::ZERO }
    }

    /// The group law: `(x0,y0)·(x1,y1) = (x0x1 - y0y1, x0y1 + y0x1)`.
    #[inline]
    pub fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x * rhs.x - self.y * rhs.y,
            y: self.x * rhs.y + self.y * rhs.x,
        }
    }

    /// The group inverse / conjugation `J(x, y) = (x, -y)`.
    #[inline]
    pub fn conjugate(self) -> Self {
        Self { x: self.x, y: -self.y }
    }

    /// The squaring endomorphism `π(x, y) = (2x^2 - 1, 2xy)`.
    #[inline]
    pub fn double(self) -> Self {
        let two = F::ONE + F::ONE;
        Self {
            x: two * self.x * self.x - F::ONE,
            y: two * self.x * self.y,
        }
    }

    /// Group exponentiation by square-and-multiply.
    pub fn mul(self, mut exp: u128) -> Self {
        let mut base = self;
        let mut result = Self::identity();
        while exp > 0 {
            if exp & 1 == 1 {
                result = result.add(base);
            }
            base = base.double();
            exp >>= 1;
        }
        result
    }

    /// `self - rhs` in the group.
    #[inline]
    pub fn sub(self, rhs: Self) -> Self {
        self.add(rhs.conjugate())
    }
}

impl CirclePoint<M31> {
    /// Lift an M31 point into the challenge field.
    pub fn into_qm31(self) -> CirclePoint<QM31> {
        CirclePoint { x: self.x.into(), y: self.y.into() }
    }
}

/// log2 of the order of the M31 circle group: `p + 1 = 2^31`.
pub const M31_CIRCLE_LOG_ORDER: u32 = 31;

/// A generator of the full M31 circle group (order `2^31`).
pub const M31_CIRCLE_GEN: CirclePoint<M31> = CirclePoint { x: M31(2), y: M31(1268011823) };

/// A point index: the exponent of [`M31_CIRCLE_GEN`], an integer mod `2^31`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PointIndex(pub u32);

const INDEX_MASK: u32 = (1 << M31_CIRCLE_LOG_ORDER) - 1;

#[allow(clippy::should_implement_trait)] // index arithmetic mirrors the group ops
impl PointIndex {
    pub const ZERO: Self = Self(0);

    /// The index of a generator of the subgroup `G_n` of size `2^n`.
    pub fn subgroup_gen(log_size: u32) -> Self {
        assert!(log_size <= M31_CIRCLE_LOG_ORDER);
        Self(1 << (M31_CIRCLE_LOG_ORDER - log_size))
    }

    #[inline]
    pub fn add(self, rhs: Self) -> Self {
        Self((self.0.wrapping_add(rhs.0)) & INDEX_MASK)
    }

    #[inline]
    pub fn neg(self) -> Self {
        Self(self.0.wrapping_neg() & INDEX_MASK)
    }

    #[inline]
    pub fn mul(self, scalar: u32) -> Self {
        Self(self.0.wrapping_mul(scalar) & INDEX_MASK)
    }

    /// Materialize the point `GEN^index`.
    pub fn to_point(self) -> CirclePoint<M31> {
        M31_CIRCLE_GEN.mul(self.0 as u128)
    }
}

/// The coset `initial · <step>` of size `2^log_size`, as indices.
///
/// `elements` are `initial + i·step` for `i = 0..2^log_size`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Coset {
    pub initial: PointIndex,
    pub step: PointIndex,
    pub log_size: u32,
}

impl Coset {
    pub fn new(initial: PointIndex, log_size: u32) -> Self {
        Self { initial, step: PointIndex::subgroup_gen(log_size), log_size }
    }

    /// The subgroup `G_n` itself.
    pub fn subgroup(log_size: u32) -> Self {
        Self::new(PointIndex::ZERO, log_size)
    }

    /// The *standard position coset* of size `2^log_size`: `Q·G_n` with `Q` of
    /// order `2^(log_size+1)` (Proposition 1). It is unique, and standard
    /// position cosets of different sizes are disjoint (all points of `Q·G_n`
    /// have exact order `2^(log_size+1)`).
    pub fn standard(log_size: u32) -> Self {
        assert!(log_size < M31_CIRCLE_LOG_ORDER, "unsupported domain size");
        Self::new(PointIndex::subgroup_gen(log_size + 1), log_size)
    }

    pub fn size(&self) -> usize {
        1 << self.log_size
    }

    pub fn index_at(&self, i: usize) -> PointIndex {
        self.initial.add(self.step.mul(i as u32))
    }

    pub fn at(&self, i: usize) -> CirclePoint<M31> {
        self.index_at(i).to_point()
    }

    /// The coset of half the size `{initial + i·2·step}`.
    pub fn half(&self) -> Coset {
        assert!(self.log_size > 0);
        Coset {
            initial: self.initial,
            step: self.step.mul(2),
            log_size: self.log_size - 1,
        }
    }

    /// Image under the squaring map π: a coset of half the size.
    pub fn double(&self) -> Coset {
        assert!(self.log_size > 0);
        Coset {
            initial: self.initial.mul(2),
            step: self.step.mul(2),
            log_size: self.log_size - 1,
        }
    }

    /// All points of the coset, in index order. Parallelized by jump-ahead
    /// chunking for large cosets.
    pub fn points(&self) -> Vec<CirclePoint<M31>> {
        let size = self.size();
        let step = self.step.to_point();

        #[cfg(feature = "parallel")]
        if size >= 1 << 14 {
            use rayon::prelude::*;
            const CHUNK: usize = 1 << 12;
            return (0..size.div_ceil(CHUNK))
                .into_par_iter()
                .flat_map_iter(|chunk| {
                    let start = chunk * CHUNK;
                    let mut cur = self.index_at(start).to_point();
                    (start..(start + CHUNK).min(size)).map(move |_| {
                        let p = cur;
                        cur = cur.add(step);
                        p
                    })
                })
                .collect();
        }

        let mut cur = self.initial.to_point();
        let mut out = Vec::with_capacity(size);
        for _ in 0..size {
            out.push(cur);
            cur = cur.add(step);
        }
        out
    }
}

/// A circle-FFT evaluation domain: the twin-coset
/// `half_coset ∪ J(half_coset)` of size `2^log_size`.
///
/// The point ordering follows Appendix D of the paper:
/// `P(i) = Q·G^i` for `i < half`, and `P(i) = J(P(i - half))` for `i >= half`,
/// which makes the circle FFT a standard radix-2 butterfly network and makes
/// all fold "siblings" sit at `i ± half`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CircleDomain {
    pub half_coset: Coset,
}

impl CircleDomain {
    /// The standard position coset of size `2^log_size` as a domain.
    /// Its half coset is `Q·G_{n-1}` with `Q` of order `2^(n+1)`.
    pub fn standard(log_size: u32) -> Self {
        assert!(log_size >= 1);
        let full = Coset::standard(log_size);
        // Q·G_n = Q·G_{n-1} ∪ Q^{-1}·G_{n-1}; the half coset keeps the same
        // initial point but steps by a generator of G_{n-1}.
        Self { half_coset: full.half() }
    }

    pub fn log_size(&self) -> u32 {
        self.half_coset.log_size + 1
    }

    pub fn size(&self) -> usize {
        1 << self.log_size()
    }

    pub fn half_size(&self) -> usize {
        self.half_coset.size()
    }

    pub fn index_at(&self, i: usize) -> PointIndex {
        let half = self.half_size();
        if i < half {
            self.half_coset.index_at(i)
        } else {
            self.half_coset.index_at(i - half).neg()
        }
    }

    pub fn at(&self, i: usize) -> CirclePoint<M31> {
        self.index_at(i).to_point()
    }

    /// All domain points in FFT order.
    pub fn points(&self) -> Vec<CirclePoint<M31>> {
        let mut out = self.half_coset.points();
        out.extend(self.half_coset.points().iter().map(|p| p.conjugate()));
        out
    }

    /// The index of `J(P(i))` — the sibling in the first FFT/FRI fold.
    pub fn conjugate_index(&self, i: usize) -> usize {
        let half = self.half_size();
        if i < half { i + half } else { i - half }
    }

    /// The index of `T_r(P(i))` where `T_r` is rotation by `r` steps of the
    /// *trace* subgroup embedded in this domain. `steps` is measured in units
    /// of this domain's half-coset step.
    ///
    /// Rotation by `g^k` maps `P(i) -> P(i + k)` on the first half and
    /// `P(i) -> P(i - k)` (mod half) on the second half.
    pub fn rotate_index(&self, i: usize, steps: usize) -> usize {
        let half = self.half_size();
        if i < half {
            (i + steps) % half
        } else {
            half + (i - half + half - (steps % half)) % half
        }
    }
}

/// Evaluate the vanishing polynomial `v_n` of the standard position coset of
/// size `2^log_size` at x-coordinate `x`: apply `x -> 2x^2 - 1` `log_size - 1`
/// times (Section 3.3). Succinct: `O(log_size)` field operations.
pub fn vanishing_eval<F: Field>(log_size: u32, x: F) -> F {
    assert!(log_size >= 1);
    let two = F::ONE + F::ONE;
    let mut x = x;
    for _ in 0..(log_size - 1) {
        x = two * x * x - F::ONE;
    }
    x
}

/// Evaluate the line (pair-vanishing polynomial) through two distinct circle
/// points `p` and `q` at `at`: a degree-1 polynomial in `L_2` that vanishes
/// exactly at `p` and `q` on the curve.
pub fn line_eval<F: Field>(p: CirclePoint<F>, q: CirclePoint<F>, at: CirclePoint<F>) -> F {
    (at.x - p.x) * (q.y - p.y) - (at.y - p.y) * (q.x - p.x)
}

/// The single-point vanishing function `v_z` of Proposition 4:
/// `v_z(P) = 1 - ((P·z^{-1}).x + i·(P·z^{-1}).y)`,
/// an `F(i)`-rational function with a single simple zero at `z`. Since our
/// challenge field QM31 contains `i`, the value lives in QM31 directly.
pub fn point_vanishing_eval(z: CirclePoint<QM31>, at: CirclePoint<QM31>) -> QM31 {
    let rel = at.sub(z);
    QM31::ONE - rel.x - QM31::I * rel.y
}

/// Evaluate the interpolant used by boundary constraints: the constant `value`
/// viewed as the claimed evaluation at the boundary point. (For a single
/// boundary point the "interpolant" is just the constant.)
pub fn boundary_interpolant_eval(value: M31) -> QM31 {
    QM31::from_base(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generator_has_exact_order_2_31() {
        let g = M31_CIRCLE_GEN;
        assert_eq!(g.x * g.x + g.y * g.y, M31::ONE, "generator not on the curve");
        assert_eq!(g.mul(1u128 << 31), CirclePoint::identity());
        assert_ne!(g.mul(1u128 << 30), CirclePoint::identity());
    }

    #[test]
    fn group_axioms() {
        let a = M31_CIRCLE_GEN.mul(1234567);
        let b = M31_CIRCLE_GEN.mul(89101112);
        let c = M31_CIRCLE_GEN.mul(3141592653);
        assert_eq!(a.add(b), b.add(a));
        assert_eq!(a.add(b.add(c)), a.add(b).add(c));
        assert_eq!(a.add(a.conjugate()), CirclePoint::identity());
        assert_eq!(a.double(), a.add(a));
        // All points are on the curve.
        for p in [a, b, c, a.add(b)] {
            assert_eq!(p.x * p.x + p.y * p.y, M31::ONE);
        }
    }

    #[test]
    fn standard_coset_points_have_exact_order() {
        // All points of the standard position coset of size 2^n have exact
        // order 2^(n+1) — this is why standard cosets of different sizes are
        // disjoint (trace domain vs evaluation domain).
        let n = 4;
        let coset = Coset::standard(n);
        for p in coset.points() {
            assert_eq!(p.mul(1 << (n + 1)), CirclePoint::identity());
            assert_ne!(p.mul(1 << n), CirclePoint::identity());
        }
    }

    #[test]
    fn domain_ordering_and_conjugation() {
        let domain = CircleDomain::standard(5);
        let pts = domain.points();
        assert_eq!(pts.len(), 32);
        // Sibling structure: J(P(i)) = P(i + half).
        for i in 0..32 {
            let j = domain.conjugate_index(i);
            assert_eq!(pts[i].conjugate(), pts[j]);
            assert_eq!(domain.at(i), pts[i]);
        }
        // All points distinct.
        for i in 0..32 {
            for j in i + 1..32 {
                assert_ne!(pts[i], pts[j]);
            }
        }
    }

    #[test]
    fn rotate_index_matches_point_rotation() {
        let domain = CircleDomain::standard(6);
        let pts = domain.points();
        let step = domain.half_coset.step.to_point();
        for steps in [1usize, 2, 3, 7] {
            let rot = step.mul(steps as u128);
            for i in 0..pts.len() {
                let ri = domain.rotate_index(i, steps);
                assert_eq!(pts[i].add(rot), pts[ri], "i={i} steps={steps}");
            }
        }
    }

    #[test]
    fn domain_double_covers_squares() {
        let domain = CircleDomain::standard(5);
        let squared: Vec<_> = domain.points().iter().map(|p| p.double()).collect();
        let smaller = CircleDomain::standard(4);
        let target = smaller.points();
        for s in &squared {
            assert!(target.contains(s), "π(D_n) should be D_(n-1)");
        }
    }

    #[test]
    fn vanishing_poly_vanishes_exactly_on_its_coset() {
        for n in 1..=6u32 {
            let coset = Coset::standard(n);
            for p in coset.points() {
                assert_eq!(vanishing_eval(n, p.x), M31::ZERO, "v_{n} nonzero on coset");
            }
            // Nonzero on the standard coset of double size (the LDE domain).
            let bigger = Coset::standard(n + 1);
            for p in bigger.points() {
                assert_ne!(vanishing_eval(n, p.x), M31::ZERO, "v_{n} zero off coset");
            }
        }
    }

    #[test]
    fn line_vanishes_at_its_two_points_only() {
        let coset = Coset::standard(4);
        let pts = coset.points();
        let p = pts[3];
        let q = pts[7];
        for (i, r) in pts.iter().enumerate() {
            let v = line_eval(p, q, *r);
            if i == 3 || i == 7 {
                assert_eq!(v, M31::ZERO);
            } else {
                assert_ne!(v, M31::ZERO);
            }
        }
    }

    #[test]
    fn point_vanishing_zero_only_at_point() {
        let z = M31_CIRCLE_GEN.mul(987654321).into_qm31();
        assert_eq!(point_vanishing_eval(z, z), QM31::ZERO);
        let other = M31_CIRCLE_GEN.mul(123).into_qm31();
        assert_ne!(point_vanishing_eval(z, other), QM31::ZERO);
    }

    #[test]
    fn qm31_circle_points_via_stereographic_projection() {
        use rand::SeedableRng;
        let mut rng = rand::rngs::StdRng::seed_from_u64(5);
        for _ in 0..20 {
            let t = QM31::random(&mut rng);
            let denom = QM31::ONE + t * t;
            if denom.is_zero() {
                continue;
            }
            let inv = denom.inverse();
            let p = CirclePoint::<QM31> {
                x: (QM31::ONE - t * t) * inv,
                y: (t + t) * inv,
            };
            assert_eq!(p.x * p.x + p.y * p.y, QM31::ONE);
        }
    }
}
