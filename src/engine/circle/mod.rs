//! The circle STARK engine over Mersenne-31 (the `circle` feature).
//!
//! Ports the circle-native geometry (circle group, twin-coset domains,
//! circle FFT, circle FRI with the dimension-gap decomposition) behind the
//! generic [`Engine`] interface, so it runs the same AIRs and shares the
//! prover/verifier/transcript/Merkle machinery with the classical engine.

pub mod backend;
pub mod cfft;
pub mod fri;
pub mod geometry;

#[cfg(feature = "metal")]
pub mod metal;

use cfft::{embed_coeffs, eval_at_point, Twiddles};
use geometry::{
    line_eval, point_vanishing_eval, vanishing_eval, CircleDomain, CirclePoint, PointIndex,
};

use crate::engine::{Engine, Mode};
use crate::field::{batch_inverse, M31, QM31};
use crate::merkle::MerkleTree;
use crate::proof::ProofOptions;
use crate::transcript::Transcript;

/// The circle M31 STARK engine.
#[derive(Clone, Copy)]
pub struct CircleEngine;

/// Draw a uniform QM31 challenge (one squeeze → four M31 limbs).
pub(crate) fn draw_ext_shared(t: &mut Transcript) -> QM31 {
    let b = t.squeeze_bytes();
    t.ratchet();
    QM31::from_m31_array([
        M31::from_bytes_mod_order(&b[0..4]),
        M31::from_bytes_mod_order(&b[4..8]),
        M31::from_bytes_mod_order(&b[8..12]),
        M31::from_bytes_mod_order(&b[12..16]),
    ])
}

impl Engine for CircleEngine {
    type Base = M31;
    type Ext = QM31;
    type Point = CirclePoint<M31>;
    type Ood = CirclePoint<QM31>;
    type Transform = Twiddles;
    type FriProver = fri::FriProver;
    type FriProof = fri::FriProof;
    type FriQueryProof = fri::FriQueryProof;

    fn mode() -> Mode {
        Mode::Circle
    }
    fn label() -> &'static str {
        "circle"
    }

    fn commit(packed: &[u8], leaf_len: usize) -> MerkleTree {
        backend::merkle_packed(packed, leaf_len)
    }

    fn log_eval(log_trace: u32, options: &ProofOptions) -> u32 {
        log_trace + 1 + options.log_blowup
    }
    fn log_bound(log_trace: u32) -> u32 {
        log_trace + 1
    }
    fn rotation_step(log_trace: u32, log_eval: u32) -> usize {
        1 << (log_eval - log_trace - 1)
    }

    fn transform(log_size: u32) -> Twiddles {
        Twiddles::new(&CircleDomain::standard(log_size))
    }

    fn interpolate(values: &mut [M31], t: &Twiddles) {
        cfft::interpolate(values, t);
    }

    fn evaluate_lde(coeffs: &[M31], log_from: u32, log_eval: u32, t_eval: &Twiddles) -> Vec<M31> {
        let mut lde = embed_coeffs(coeffs, log_from, log_eval);
        backend::evaluate_m31(&mut lde, t_eval);
        lde
    }

    fn eval_at_ood(coeffs: &[M31], ood: CirclePoint<QM31>) -> QM31 {
        eval_at_point(coeffs, ood)
    }

    fn eval_points(log_eval: u32) -> Vec<CirclePoint<M31>> {
        CircleDomain::standard(log_eval).points()
    }

    fn fft_index_of_row(log_trace: u32, row: usize) -> usize {
        let n = 1usize << log_trace;
        if row.is_multiple_of(2) {
            row / 2
        } else {
            n / 2 + (n - 1 - row) / 2
        }
    }

    fn trace_point(log_trace: u32, row: usize) -> CirclePoint<M31> {
        let q = PointIndex::subgroup_gen(log_trace + 1);
        q.add(PointIndex::subgroup_gen(log_trace).mul(row as u32))
            .to_point()
    }

    fn rotation_index(log_eval: u32, i: usize, steps: usize) -> usize {
        CircleDomain::standard(log_eval).rotate_index(i, steps)
    }

    fn vanishing_over_eval(log_trace: u32, points: &[CirclePoint<M31>]) -> Vec<M31> {
        points.iter().map(|p| vanishing_eval(log_trace, p.x)).collect()
    }

    fn vanishing_at_ood(log_trace: u32, ood: CirclePoint<QM31>) -> QM31 {
        vanishing_eval(log_trace, ood.x)
    }

    fn selector_over_eval(
        excluded: &[CirclePoint<M31>],
        points: &[CirclePoint<M31>],
    ) -> Option<Vec<M31>> {
        match excluded.len() {
            0 => None,
            2 => Some(
                points
                    .iter()
                    .map(|&p| line_eval(excluded[0], excluded[1], p))
                    .collect(),
            ),
            k => panic!("unsupported exclusion size {k}"),
        }
    }

    fn selector_at_ood(excluded: &[CirclePoint<M31>], ood: CirclePoint<QM31>) -> QM31 {
        line_eval(excluded[0].into_qm31(), excluded[1].into_qm31(), ood)
    }

    fn draw_ood(t: &mut Transcript) -> CirclePoint<QM31> {
        loop {
            let z = draw_ext_shared(t);
            let denom = QM31::ONE + z * z;
            if denom.is_zero() {
                continue;
            }
            let inv = denom.inverse();
            let p = CirclePoint {
                x: (QM31::ONE - z * z) * inv,
                y: (z + z) * inv,
            };
            if p.x.is_base() && p.y.is_base() {
                continue;
            }
            return p;
        }
    }
    fn draw_ext(t: &mut Transcript) -> QM31 {
        draw_ext_shared(t)
    }
    fn absorb_ext(t: &mut Transcript, v: QM31) {
        t.absorb_field(v);
    }

    fn mask_point(ood: CirclePoint<QM31>, offset: usize, log_trace: u32) -> CirclePoint<QM31> {
        let g = PointIndex::subgroup_gen(log_trace).to_point().into_qm31();
        let mut p = ood;
        for _ in 0..offset {
            p = p.add(g);
        }
        p
    }

    fn deep_denom_inv_over_eval(
        ood: CirclePoint<QM31>,
        points: &[CirclePoint<M31>],
    ) -> Vec<QM31> {
        let denom: Vec<QM31> = points
            .iter()
            .map(|&p| point_vanishing_eval(ood, p.into_qm31()))
            .collect();
        batch_inverse(&denom)
    }
    fn deep_denom_inv_at(ood: CirclePoint<QM31>, point: CirclePoint<M31>) -> QM31 {
        point_vanishing_eval(ood, point.into_qm31()).inverse()
    }

    fn boundary_denom_inv_over_eval(
        boundary_point: CirclePoint<M31>,
        points: &[CirclePoint<M31>],
    ) -> Vec<QM31> {
        let z = boundary_point.into_qm31();
        let denom: Vec<QM31> = points
            .iter()
            .map(|&p| point_vanishing_eval(z, p.into_qm31()))
            .collect();
        batch_inverse(&denom)
    }
    fn boundary_denom_inv_at_ood(
        boundary_point: CirclePoint<M31>,
        ood: CirclePoint<QM31>,
    ) -> QM31 {
        point_vanishing_eval(boundary_point.into_qm31(), ood).inverse()
    }

    fn fri_commit(
        transcript: &mut Transcript,
        word: Vec<QM31>,
        log_eval: u32,
        log_bound: u32,
        t_eval: &Twiddles,
    ) -> fri::FriProver {
        let domain = CircleDomain::standard(log_eval);
        fri::fri_commit(transcript, word, &domain, &t_eval.inverse, log_bound)
    }
    fn fri_into_proof(prover: fri::FriProver, query_indices: &[usize]) -> fri::FriProof {
        prover.into_proof(query_indices)
    }
    fn fri_replay(transcript: &mut Transcript, proof: &fri::FriProof, log_bound: u32) -> Vec<QM31> {
        fri::fri_replay_transcript(transcript, proof, log_bound)
    }
    fn fri_check_shape(
        proof: &fri::FriProof,
        log_eval: u32,
        log_bound: u32,
        num_queries: usize,
    ) -> bool {
        fri::fri_check_shape(proof, &CircleDomain::standard(log_eval), log_bound, num_queries)
    }
    fn fri_query_proofs(proof: &fri::FriProof) -> &[fri::FriQueryProof] {
        &proof.query_proofs
    }
    #[allow(clippy::too_many_arguments)]
    fn fri_verify_query(
        proof: &fri::FriProof,
        log_eval: u32,
        log_bound: u32,
        challenges: &[QM31],
        query_index: usize,
        word_q: QM31,
        word_sib: QM31,
        qp: &fri::FriQueryProof,
    ) -> bool {
        let domain = CircleDomain::standard(log_eval);
        // Remove the dimension-gap component before folding (v_b depends only
        // on x, shared by a point and its conjugate sibling).
        let vb = vanishing_eval(log_bound, domain.at(query_index).x);
        let g = word_q - proof.lambda.mul_base(vb);
        let g_sib = word_sib - proof.lambda.mul_base(vb);
        fri::fri_verify_query(proof, &domain, challenges, query_index, g, g_sib, qp)
    }
}
