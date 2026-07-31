//! The classical roots-of-unity engine over BabyBear.
//!
//! This is the legacy toyni STARK — NTT interpolation, multiplicative-coset
//! evaluation domains, 2-adic FRI — refactored to plug into the generic
//! [`Engine`] interface so it runs any [`crate::air::Air`] and shares the
//! prover/verifier/transcript/Merkle machinery with the circle engine.

pub mod domain;
pub mod fri;
pub mod fri_fold;
pub mod ntt;

use domain::{coset_elements, BabyBearDomain, COSET_SHIFT};

use crate::engine::{Engine, Mode};
use crate::field::{batch_inverse, BabyBear, BabyBearExt, ExtField, Field};
use crate::proof::ProofOptions;
use crate::transcript::Transcript;

type Ext = BabyBearExt;

/// The classical BabyBear STARK engine.
#[derive(Clone, Copy)]
pub struct ClassicalEngine;

/// Precomputed transform data for a domain of a given size.
pub struct Transform {
    pub log_size: u32,
    pub omega: BabyBear,
}

/// Draw a uniform BabyBear quartic-extension challenge (four base squeezes).
pub(crate) fn draw_ext(t: &mut Transcript) -> Ext {
    Ext::new([
        t.draw_field::<BabyBear>(),
        t.draw_field::<BabyBear>(),
        t.draw_field::<BabyBear>(),
        t.draw_field::<BabyBear>(),
    ])
}

impl Engine for ClassicalEngine {
    type Base = BabyBear;
    type Ext = Ext;
    type Point = BabyBear;
    type Ood = Ext;
    type Transform = Transform;
    type FriProver = fri::FriProver;
    type FriProof = fri::FriProof;
    type FriQueryProof = fri::FriQueryProof;

    fn mode() -> Mode {
        Mode::Classical
    }
    fn label() -> &'static str {
        "roots-of-unity"
    }

    fn log_eval(log_trace: u32, options: &ProofOptions) -> u32 {
        log_trace + 1 + options.log_blowup
    }
    fn log_bound(log_trace: u32) -> u32 {
        log_trace + 1
    }
    fn rotation_step(log_trace: u32, log_eval: u32) -> usize {
        1 << (log_eval - log_trace)
    }

    fn transform(log_size: u32) -> Transform {
        Transform {
            log_size,
            omega: BabyBear::get_root_of_unity(log_size),
        }
    }

    fn interpolate(values: &mut [BabyBear], t: &Transform) {
        ntt::intt(values, t.omega);
    }

    fn evaluate_lde(
        coeffs: &[BabyBear],
        _log_from: u32,
        log_eval: u32,
        _t_eval: &Transform,
    ) -> Vec<BabyBear> {
        let coset = BabyBearDomain::new(1usize << log_eval).get_coset(BabyBear::new(COSET_SHIFT));
        coset.fft(coeffs)
    }

    fn eval_at_ood(coeffs: &[BabyBear], ood: Ext) -> Ext {
        ood.eval_base_poly(coeffs)
    }

    fn interpolate_ext(values: &mut [Ext], t: &Transform) {
        let d = BabyBearDomain::new(1usize << t.log_size);
        let coeffs = d.ifft_ext(values);
        values.copy_from_slice(&coeffs);
    }
    fn evaluate_lde_ext(coeffs: &[Ext], _log_from: u32, log_eval: u32, _t: &Transform) -> Vec<Ext> {
        let coset = BabyBearDomain::new(1usize << log_eval).get_coset(BabyBear::new(COSET_SHIFT));
        let mut padded = coeffs.to_vec();
        padded.resize(1usize << log_eval, Ext::ZERO);
        coset.fft_ext(&padded)
    }
    fn eval_ext_at_ood(coeffs: &[Ext], ood: Ext) -> Ext {
        ood.eval_ext_poly(coeffs)
    }
    fn point_vanishing_over_eval(trace_point: BabyBear, points: &[BabyBear]) -> Vec<Ext> {
        points.iter().map(|&x| Ext::from(x - trace_point)).collect()
    }
    fn point_vanishing_at_ood(trace_point: BabyBear, ood: Ext) -> Ext {
        ood - Ext::from(trace_point)
    }

    fn eval_points(log_eval: u32) -> Vec<BabyBear> {
        coset_elements(log_eval)
    }

    fn fft_index_of_row(_log_trace: u32, row: usize) -> usize {
        row
    }

    fn trace_point(log_trace: u32, row: usize) -> BabyBear {
        BabyBear::get_root_of_unity(log_trace).pow(row as u64)
    }

    fn rotation_index(log_eval: u32, i: usize, steps: usize) -> usize {
        (i + steps) % (1usize << log_eval)
    }

    fn vanishing_over_eval(log_trace: u32, points: &[BabyBear]) -> Vec<BabyBear> {
        let n = 1u64 << log_trace;
        points.iter().map(|x| x.pow(n) - BabyBear::ONE).collect()
    }

    fn vanishing_at_ood(log_trace: u32, ood: Ext) -> Ext {
        ood.pow(1u64 << log_trace) - Ext::ONE
    }

    fn selector_over_eval(
        excluded: &[BabyBear],
        points: &[BabyBear],
    ) -> Option<Vec<BabyBear>> {
        if excluded.is_empty() {
            return None;
        }
        Some(
            points
                .iter()
                .map(|&x| excluded.iter().fold(BabyBear::ONE, |acc, &e| acc * (x - e)))
                .collect(),
        )
    }

    fn selector_at_ood(excluded: &[BabyBear], ood: Ext) -> Ext {
        excluded
            .iter()
            .fold(Ext::ONE, |acc, &e| acc * (ood - Ext::from(e)))
    }

    fn draw_ood(t: &mut Transcript) -> Ext {
        loop {
            let z = draw_ext(t);
            if !z.is_base() {
                return z;
            }
        }
    }
    fn draw_ext(t: &mut Transcript) -> Ext {
        draw_ext(t)
    }
    fn absorb_ext(t: &mut Transcript, v: Ext) {
        t.absorb_field(v);
    }

    fn mask_point(ood: Ext, offset: usize, log_trace: u32) -> Ext {
        let g = BabyBear::get_root_of_unity(log_trace).pow(offset as u64);
        ood * Ext::from(g)
    }

    fn deep_denom_inv_over_eval(ood: Ext, points: &[BabyBear]) -> Vec<Ext> {
        let denom: Vec<Ext> = points.iter().map(|&x| Ext::from(x) - ood).collect();
        batch_inverse(&denom)
    }
    fn deep_denom_inv_at(ood: Ext, point: BabyBear) -> Ext {
        (Ext::from(point) - ood).inverse()
    }

    fn boundary_denom_inv_over_eval(boundary_point: BabyBear, points: &[BabyBear]) -> Vec<Ext> {
        let denom: Vec<BabyBear> = points.iter().map(|&x| x - boundary_point).collect();
        batch_inverse(&denom).into_iter().map(Ext::from).collect()
    }
    fn boundary_denom_inv_at_ood(boundary_point: BabyBear, ood: Ext) -> Ext {
        (ood - Ext::from(boundary_point)).inverse()
    }

    fn fri_commit(
        transcript: &mut Transcript,
        word: Vec<Ext>,
        log_eval: u32,
        log_bound: u32,
        _t_eval: &Transform,
    ) -> fri::FriProver {
        fri::commit(transcript, word, log_eval, log_bound)
    }
    fn fri_into_proof(prover: fri::FriProver, query_indices: &[usize]) -> fri::FriProof {
        prover.into_proof(query_indices)
    }
    fn fri_replay(transcript: &mut Transcript, proof: &fri::FriProof, log_bound: u32) -> Vec<Ext> {
        fri::replay(transcript, proof, log_bound)
    }
    fn fri_check_shape(
        proof: &fri::FriProof,
        log_eval: u32,
        log_bound: u32,
        num_queries: usize,
    ) -> bool {
        fri::check_shape(proof, log_eval, log_bound, num_queries)
    }
    fn fri_query_proofs(proof: &fri::FriProof) -> &[fri::FriQueryProof] {
        &proof.query_proofs
    }
    fn fri_verify_query(
        proof: &fri::FriProof,
        log_eval: u32,
        _log_bound: u32,
        challenges: &[Ext],
        query_index: usize,
        word_q: Ext,
        word_sib: Ext,
        qp: &fri::FriQueryProof,
    ) -> bool {
        fri::verify_query(proof, log_eval, challenges, query_index, word_q, word_sib, qp)
    }
}
