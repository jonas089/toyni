//! The proving-engine abstraction: the geometry, transforms and low-degree
//! test that a STARK is built on, factored so the generic prover/verifier and
//! the example AIRs are written once and run on either engine.
//!
//! Two engines implement it:
//!
//! - [`classical::ClassicalEngine`] — the roots-of-unity STARK over BabyBear
//!   (NTT, multiplicative-coset domains, 2-adic FRI); the legacy toyni mode,
//!   now generic.
//! - [`circle::CircleEngine`] — the circle STARK over M31 (circle FFT,
//!   twin-coset domains, circle FRI with the dimension-gap decomposition);
//!   the `circle` feature.
//!
//! Both realize the same DEEP-ALI protocol shape (commit trace, commit the
//! constraint composition, out-of-domain sample, batch single-point DEEP
//! quotients, FRI to a constant final layer). Everything that differs between
//! roots-of-unity and circle geometry — domain points, FFT, vanishing
//! polynomials, single-point quotient denominators, the out-of-domain point
//! type, and the FRI fold — lives behind this trait.

pub mod classical;
#[cfg(feature = "circle")]
pub mod circle;

use crate::field::{ExtField, Field};
use crate::proof::ProofOptions;
use crate::transcript::Transcript;

/// Which engine a proof is produced/checked with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// Roots-of-unity STARK over BabyBear.
    Classical,
    /// Circle STARK over M31 (`circle` feature).
    Circle,
}

/// A pluggable STARK engine: field tower + domain geometry + FFT + FRI.
pub trait Engine: Sized + 'static {
    /// Base field the trace lives in.
    type Base: Field;
    /// Extension field the challenges live in (~124-bit soundness).
    type Ext: ExtField<Self::Base> + crate::field::ConstraintField<Self::Base>;
    /// An evaluation-domain point in the base field (a coset element `x`, or a
    /// circle point over M31).
    type Point: Copy + Send + Sync;
    /// An out-of-domain point in the extension field (an `Ext` scalar `z`, or
    /// a circle point over QM31).
    type Ood: Copy + Send + Sync;
    /// Precomputed transform data (twiddles) for a domain of a given size.
    type Transform: Send + Sync;
    /// FRI prover state carried between the commit and query phases.
    type FriProver;
    /// The serialized FRI proof embedded in a [`crate::proof::StarkProof`].
    type FriProof: Clone + std::fmt::Debug + Send + Sync;
    /// Per-query FRI openings.
    type FriQueryProof: Clone + std::fmt::Debug + Send + Sync;

    fn mode() -> Mode;
    fn label() -> &'static str;

    /// Commit to packed fixed-size leaves. Default is the CPU Merkle tree;
    /// the circle engine overrides this to use its GPU hasher when available.
    fn commit(packed: &[u8], leaf_len: usize) -> crate::merkle::MerkleTree {
        crate::merkle::MerkleTree::from_packed(packed, leaf_len)
    }

    // ── protocol sizes ─────────────────────────────────────────────────
    /// log2 of the evaluation-domain size.
    fn log_eval(log_trace: u32, options: &ProofOptions) -> u32;
    /// FRI degree-bound exponent (folds this many times to a constant layer).
    fn log_bound(log_trace: u32) -> u32;
    /// One trace-row rotation, in units of the evaluation domain's index step.
    fn rotation_step(log_trace: u32, log_eval: u32) -> usize;
    /// The index space queries are drawn from (`|D|/2` for both engines).
    fn query_space(log_eval: u32) -> usize {
        (1usize << log_eval) / 2
    }

    // ── transforms ─────────────────────────────────────────────────────
    fn transform(log_size: u32) -> Self::Transform;
    /// FFT-order trace values → basis coefficients, in place.
    fn interpolate(values: &mut [Self::Base], t: &Self::Transform);
    /// Basis coefficients (of a size-`2^log_trace` object) → low-degree
    /// extension values over the size-`2^log_eval` evaluation domain.
    fn evaluate_lde(
        coeffs: &[Self::Base],
        log_trace: u32,
        log_eval: u32,
        t_eval: &Self::Transform,
    ) -> Vec<Self::Base>;
    /// Evaluate a base-coefficient polynomial at an out-of-domain point.
    fn eval_at_ood(coeffs: &[Self::Base], ood: Self::Ood) -> Self::Ext;

    // Extension-field transforms (for auxiliary columns such as the zkVM's
    // permutation / lookup accumulators, which are committed over `Ext`).
    /// FFT-order extension values → basis coefficients, in place.
    fn interpolate_ext(values: &mut [Self::Ext], t: &Self::Transform);
    /// Extension coefficients → LDE values over the evaluation domain.
    fn evaluate_lde_ext(
        coeffs: &[Self::Ext],
        log_from: u32,
        log_eval: u32,
        t_eval: &Self::Transform,
    ) -> Vec<Self::Ext>;
    /// Evaluate an extension-coefficient polynomial at an out-of-domain point.
    fn eval_ext_at_ood(coeffs: &[Self::Ext], ood: Self::Ood) -> Self::Ext;

    /// Forward single-point vanishing function of a trace point, over the
    /// evaluation domain (a numerator factor for single-row exclusion). It is
    /// the reciprocal of [`Engine::boundary_denom_inv_over_eval`].
    fn point_vanishing_over_eval(
        trace_point: Self::Point,
        points: &[Self::Point],
    ) -> Vec<Self::Ext>;
    /// Forward single-point vanishing of a trace point at an out-of-domain point.
    fn point_vanishing_at_ood(trace_point: Self::Point, ood: Self::Ood) -> Self::Ext;

    // ── domain geometry ────────────────────────────────────────────────
    fn eval_points(log_eval: u32) -> Vec<Self::Point>;
    /// Position of trace row `row` in the interpolation (FFT) ordering.
    fn fft_index_of_row(log_trace: u32, row: usize) -> usize;
    /// The circle/coset point of trace row `row` (for boundary denominators).
    fn trace_point(log_trace: u32, row: usize) -> Self::Point;
    /// Index of `T^steps(P(i))` in the evaluation domain (mask rotation).
    fn rotation_index(log_eval: u32, i: usize, steps: usize) -> usize;
    /// The conjugate/sibling index `i ± |D|/2` (the first FRI fold pair).
    fn sibling_index(log_eval: u32, i: usize) -> usize {
        let half = 1usize << (log_eval - 1);
        if i < half {
            i + half
        } else {
            i - half
        }
    }

    // ── vanishing polynomials & exclusion selector ─────────────────────
    fn vanishing_over_eval(log_trace: u32, points: &[Self::Point]) -> Vec<Self::Base>;
    fn vanishing_at_ood(log_trace: u32, ood: Self::Ood) -> Self::Ext;
    /// Selector vanishing on the excluded rows, over the eval domain.
    /// `None` when no rows are excluded (single-row masks).
    fn selector_over_eval(
        excluded: &[Self::Point],
        points: &[Self::Point],
    ) -> Option<Vec<Self::Base>>;
    fn selector_at_ood(excluded: &[Self::Point], ood: Self::Ood) -> Self::Ext;

    // ── out-of-domain & DEEP quotients ─────────────────────────────────
    fn draw_ood(t: &mut Transcript) -> Self::Ood;
    fn draw_ext(t: &mut Transcript) -> Self::Ext;
    fn absorb_ext(t: &mut Transcript, v: Self::Ext);
    /// The mask point `ood` rotated by `offset` trace rows.
    fn mask_point(ood: Self::Ood, offset: usize, log_trace: u32) -> Self::Ood;
    /// Inverse single-point (DEEP) quotient denominators over the eval domain.
    fn deep_denom_inv_over_eval(ood: Self::Ood, points: &[Self::Point]) -> Vec<Self::Ext>;
    /// Inverse single-point (DEEP) quotient denominator at one eval point.
    fn deep_denom_inv_at(ood: Self::Ood, point: Self::Point) -> Self::Ext;
    /// Inverse boundary-quotient denominators over the eval domain.
    fn boundary_denom_inv_over_eval(
        boundary_point: Self::Point,
        points: &[Self::Point],
    ) -> Vec<Self::Ext>;
    /// Inverse boundary-quotient denominator at the out-of-domain point.
    fn boundary_denom_inv_at_ood(boundary_point: Self::Point, ood: Self::Ood) -> Self::Ext;

    // ── FRI (fully encapsulated, incl. the circle dimension-gap step) ──
    fn fri_commit(
        transcript: &mut Transcript,
        word: Vec<Self::Ext>,
        log_eval: u32,
        log_bound: u32,
        t_eval: &Self::Transform,
    ) -> Self::FriProver;
    fn fri_into_proof(prover: Self::FriProver, query_indices: &[usize]) -> Self::FriProof;
    fn fri_replay(
        transcript: &mut Transcript,
        proof: &Self::FriProof,
        log_bound: u32,
    ) -> Vec<Self::Ext>;
    fn fri_check_shape(
        proof: &Self::FriProof,
        log_eval: u32,
        log_bound: u32,
        num_queries: usize,
    ) -> bool;
    fn fri_query_proofs(proof: &Self::FriProof) -> &[Self::FriQueryProof];
    /// Verify one FRI query. `word_q` / `word_sib` are the DEEP-batched values
    /// at the query index and its sibling, recomputed by the verifier from its
    /// own Merkle-verified openings (this is what binds FRI to the trace and
    /// composition commitments). Any engine-specific pre-fold step (the circle
    /// dimension-gap removal) happens inside.
    #[allow(clippy::too_many_arguments)]
    fn fri_verify_query(
        proof: &Self::FriProof,
        log_eval: u32,
        log_bound: u32,
        challenges: &[Self::Ext],
        query_index: usize,
        word_q: Self::Ext,
        word_sib: Self::Ext,
        qp: &Self::FriQueryProof,
    ) -> bool;
}
