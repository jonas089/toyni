//! The generic STARK verifier, parameterized by the [`Engine`].
//!
//! Replays the transcript, checks the out-of-domain composition identity at
//! `gamma` (binding the committed composition to the constraints), then
//! spot-checks the DEEP/FRI chain: every query recomputes the batched DEEP
//! value from its own Merkle-verified openings at the query position and its
//! sibling, so the FRI fold input is fully bound to the commitments.

use std::marker::PhantomData;

use crate::air::Air;
use crate::engine::Engine;
use crate::field::Field;
use crate::merkle::verify_path;
use crate::proof::{ProofOptions, StarkProof};
use crate::prover::{
    composition_leaf_bytes, eval_composition_at, init_transcript, powers, trace_leaf_bytes,
};

pub struct StarkVerifier<'a, E: Engine, A: Air<E::Base>> {
    air: &'a A,
    options: ProofOptions,
    _engine: PhantomData<E>,
}

impl<'a, E: Engine, A: Air<E::Base>> StarkVerifier<'a, E, A> {
    pub fn new(air: &'a A, options: ProofOptions) -> Self {
        Self {
            air,
            options,
            _engine: PhantomData,
        }
    }

    pub fn verify(&self, proof: &StarkProof<E>) -> bool {
        let air = self.air;
        let num_cols = air.num_columns();
        let num_offsets = air.mask_offsets().len();
        let log_trace = proof.log_trace_len;
        if log_trace < 2 {
            return false;
        }
        let log_eval = E::log_eval(log_trace, &self.options);
        let log_bound = E::log_bound(log_trace);
        let domain_size = 1usize << log_eval;
        let half = domain_size / 2;

        if proof.ood_trace.len() != num_offsets
            || proof.ood_trace.iter().any(|r| r.len() != num_cols)
            || proof.queries.len() != self.options.num_queries
        {
            return false;
        }

        // ── transcript replay ───────────────────────────────────────────
        let mut transcript = init_transcript::<E, A>(air, &self.options, log_trace, log_eval);
        transcript.absorb_commitment(&proof.trace_root);
        let beta = E::draw_ext(&mut transcript);
        transcript.absorb_commitment(&proof.composition_root);
        let gamma = E::draw_ood(&mut transcript);
        let mask_points: Vec<E::Ood> = air
            .mask_offsets()
            .iter()
            .map(|&k| E::mask_point(gamma, k, log_trace))
            .collect();

        for row in &proof.ood_trace {
            for &v in row {
                E::absorb_ext(&mut transcript, v);
            }
        }
        E::absorb_ext(&mut transcript, proof.ood_composition);

        // ── out-of-domain composition identity ──────────────────────────
        if proof.ood_composition
            != eval_composition_at::<E, A>(air, log_trace, beta, gamma, &proof.ood_trace)
        {
            return false;
        }

        let mu = E::draw_ext(&mut transcript);
        let fold_challenges = E::fri_replay(&mut transcript, &proof.fri, log_bound);
        let query_indices = transcript.draw_indices(self.options.num_queries, half);

        if !E::fri_check_shape(&proof.fri, log_eval, log_bound, self.options.num_queries) {
            return false;
        }

        let points = E::eval_points(log_eval);
        let fri_qps = E::fri_query_proofs(&proof.fri);

        for (qi, &q) in query_indices.iter().enumerate() {
            let openings = &proof.queries[qi];
            let q_sib = q + half;

            for (pos, t, c) in [
                (q, &openings.trace, &openings.composition),
                (q_sib, &openings.trace_sibling, &openings.composition_sibling),
            ] {
                if t.values.len() != num_cols {
                    return false;
                }
                let leaf = trace_leaf_bytes::<E>(&t.values, &t.salt);
                if !verify_path(&proof.trace_root, pos, &leaf, &t.path) {
                    return false;
                }
                let leaf = composition_leaf_bytes::<E>(c.value, &c.salt);
                if !verify_path(&proof.composition_root, pos, &leaf, &c.path) {
                    return false;
                }
            }

            let word_q = self.deep_at(
                q,
                &openings.trace.values,
                openings.composition.value,
                &points,
                &mask_points,
                &proof.ood_trace,
                proof.ood_composition,
                mu,
            );
            let word_sib = self.deep_at(
                q_sib,
                &openings.trace_sibling.values,
                openings.composition_sibling.value,
                &points,
                &mask_points,
                &proof.ood_trace,
                proof.ood_composition,
                mu,
            );

            if !E::fri_verify_query(
                &proof.fri,
                log_eval,
                log_bound,
                &fold_challenges,
                q,
                word_q,
                word_sib,
                &fri_qps[qi],
            ) {
                return false;
            }
        }
        true
    }

    #[allow(clippy::too_many_arguments)]
    fn deep_at(
        &self,
        i: usize,
        trace_values: &[E::Base],
        comp_value: E::Ext,
        points: &[E::Point],
        mask_points: &[E::Ood],
        ood_trace: &[Vec<E::Ext>],
        ood_composition: E::Ext,
        mu: E::Ext,
    ) -> E::Ext {
        let num_cols = trace_values.len();
        let mu_pows = powers(mu, mask_points.len() * num_cols + 1);
        let inv_vz: Vec<E::Ext> = mask_points
            .iter()
            .map(|&z| E::deep_denom_inv_at(z, points[i]))
            .collect();
        let mut acc = E::Ext::ZERO;
        let mut term = 0;
        for (k, ood_row) in ood_trace.iter().enumerate() {
            for (c, &ood) in ood_row.iter().enumerate() {
                acc = acc + mu_pows[term] * (E::Ext::from(trace_values[c]) - ood) * inv_vz[k];
                term += 1;
            }
        }
        acc + mu_pows[term] * (comp_value - ood_composition) * inv_vz[0]
    }
}
