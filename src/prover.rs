//! The generic STARK prover, parameterized by the [`Engine`].
//!
//! The protocol is identical for both engines — only the geometry, FFT and
//! FRI differ, and those live behind the engine trait:
//!
//! 1. Interpolate each trace column, blind it (zero-knowledge masking with the
//!    trace-domain vanishing polynomial), low-degree-extend and commit.
//! 2. Draw `beta`, fold all constraints into the composition polynomial
//!    `F = Sum beta^i * C_i * selector / vanishing + boundary quotients`, commit.
//! 3. Draw the out-of-domain point `gamma`, open the trace mask and `F` at it.
//! 4. Draw `mu`, batch all single-point (DEEP) quotients into one word `u`.
//! 5. Run FRI on `u` and answer spot-check queries.

use std::marker::PhantomData;

use rand::Rng;

use crate::air::Air;
use crate::air::TraceTable;
use crate::engine::Engine;
use crate::field::{batch_inverse, ExtField, Field};
use crate::proof::{
    CompositionOpening, ProofOptions, QueryOpenings, StarkProof, TraceOpening,
};
use crate::transcript::Transcript;

const SALT_LEN: usize = 16;

pub struct StarkProver<'a, E: Engine, A: Air<E::Base>> {
    air: &'a A,
    options: ProofOptions,
    _engine: PhantomData<E>,
}

impl<'a, E: Engine, A: Air<E::Base>> StarkProver<'a, E, A> {
    pub fn new(air: &'a A, options: ProofOptions) -> Self {
        Self {
            air,
            options,
            _engine: PhantomData,
        }
    }

    pub fn prove(&self, trace: &TraceTable<E::Base>) -> StarkProof<E> {
        let air = self.air;
        assert_eq!(trace.num_columns(), air.num_columns());
        assert!(trace.satisfies(air), "trace does not satisfy the AIR");

        let log_trace = trace.log_len();
        let log_eval = E::log_eval(log_trace, &self.options);
        let log_bound = E::log_bound(log_trace);
        let domain_size = 1usize << log_eval;
        let num_cols = air.num_columns();
        let mut rng = rand::thread_rng();

        let trace_tf = E::transform(log_trace);
        let eval_tf = E::transform(log_eval);
        let points = E::eval_points(log_eval);
        let rotation_step = E::rotation_step(log_trace, log_eval);
        let vn = E::vanishing_over_eval(log_trace, &points);

        // ── 1. trace: interpolate, mask, LDE, commit ────────────────────
        let mask_log = self.options.mask_log_size(log_trace);
        let mut trace_coeffs: Vec<Vec<E::Base>> = Vec::with_capacity(num_cols);
        let mut mask_coeffs: Vec<Vec<E::Base>> = Vec::with_capacity(num_cols);
        let mut trace_ldes: Vec<Vec<E::Base>> = Vec::with_capacity(num_cols);

        for col in &trace.columns {
            let mut values = vec![E::Base::ZERO; trace.len()];
            for (row, &v) in col.iter().enumerate() {
                values[E::fft_index_of_row(log_trace, row)] = v;
            }
            E::interpolate(&mut values, &trace_tf);
            let coeffs = values;

            let mut lde = E::evaluate_lde(&coeffs, log_trace, log_eval, &eval_tf);
            let r_coeffs: Vec<E::Base> = if self.options.zk {
                (0..1usize << mask_log)
                    .map(|_| E::Base::random(&mut rng))
                    .collect()
            } else {
                Vec::new()
            };
            if self.options.zk {
                let r_lde = E::evaluate_lde(&r_coeffs, mask_log, log_eval, &eval_tf);
                for i in 0..domain_size {
                    lde[i] += vn[i] * r_lde[i];
                }
            }
            trace_coeffs.push(coeffs);
            mask_coeffs.push(r_coeffs);
            trace_ldes.push(lde);
        }

        let salt_len = if self.options.zk { SALT_LEN } else { 0 };
        let trace_salts = salts(domain_size, self.options.zk, &mut rng);
        let trace_leaf_len = num_cols * E::Base::NUM_BYTES + salt_len;
        let mut trace_packed = vec![0u8; domain_size * trace_leaf_len];
        for i in 0..domain_size {
            let leaf = &mut trace_packed[i * trace_leaf_len..(i + 1) * trace_leaf_len];
            for (c, col) in trace_ldes.iter().enumerate() {
                let b = col[i].to_bytes();
                leaf[c * E::Base::NUM_BYTES..(c + 1) * E::Base::NUM_BYTES].copy_from_slice(&b);
            }
            leaf[num_cols * E::Base::NUM_BYTES..].copy_from_slice(&trace_salts[i]);
        }
        let trace_tree = E::commit(&trace_packed, trace_leaf_len);
        let trace_root = trace_tree.root();

        let mut transcript = init_transcript::<E, A>(air, &self.options, log_trace, log_eval);
        transcript.absorb_commitment(&trace_root);

        // ── 2. composition polynomial ───────────────────────────────────
        let beta = E::draw_ext(&mut transcript);
        let composition = self.compute_composition(
            log_trace,
            log_eval,
            rotation_step,
            &points,
            &vn,
            &trace_ldes,
            beta,
        );

        let comp_salts = salts(domain_size, self.options.zk, &mut rng);
        let comp_leaf_len = E::Ext::NUM_BYTES + salt_len;
        let mut comp_packed = vec![0u8; domain_size * comp_leaf_len];
        for i in 0..domain_size {
            let leaf = &mut comp_packed[i * comp_leaf_len..(i + 1) * comp_leaf_len];
            leaf[..E::Ext::NUM_BYTES].copy_from_slice(&composition[i].to_bytes());
            leaf[E::Ext::NUM_BYTES..].copy_from_slice(&comp_salts[i]);
        }
        let comp_tree = E::commit(&comp_packed, comp_leaf_len);
        let composition_root = comp_tree.root();
        transcript.absorb_commitment(&composition_root);

        // ── 3. out-of-domain point and openings ─────────────────────────
        let gamma = E::draw_ood(&mut transcript);
        let mask_points: Vec<E::Ood> = air
            .mask_offsets()
            .iter()
            .map(|&k| E::mask_point(gamma, k, log_trace))
            .collect();

        let ood_trace: Vec<Vec<E::Ext>> = mask_points
            .iter()
            .map(|&z| {
                (0..num_cols)
                    .map(|c| {
                        let mut v = E::eval_at_ood(&trace_coeffs[c], z);
                        if self.options.zk {
                            v = v + E::vanishing_at_ood(log_trace, z)
                                * E::eval_at_ood(&mask_coeffs[c], z);
                        }
                        v
                    })
                    .collect()
            })
            .collect();

        let ood_composition = eval_composition_at::<E, A>(air, log_trace, beta, gamma, &ood_trace);

        for row in &ood_trace {
            for &v in row {
                E::absorb_ext(&mut transcript, v);
            }
        }
        E::absorb_ext(&mut transcript, ood_composition);

        // ── 4. DEEP batching ────────────────────────────────────────────
        let mu = E::draw_ext(&mut transcript);
        let deep = self.compute_deep(
            &points,
            &trace_ldes,
            &composition,
            &mask_points,
            &ood_trace,
            ood_composition,
            mu,
        );

        // ── 5. FRI ──────────────────────────────────────────────────────
        let fri_prover = E::fri_commit(&mut transcript, deep, log_eval, log_bound, &eval_tf);

        // ── 6. queries ──────────────────────────────────────────────────
        let query_indices =
            transcript.draw_indices(self.options.num_queries, E::query_space(log_eval));
        let half = domain_size / 2;
        let queries: Vec<QueryOpenings<E>> = query_indices
            .iter()
            .map(|&q| {
                let open_trace = |i: usize| TraceOpening::<E> {
                    values: trace_ldes.iter().map(|c| c[i]).collect(),
                    salt: trace_salts[i].clone(),
                    path: trace_tree.prove(i),
                };
                let open_comp = |i: usize| CompositionOpening::<E> {
                    value: composition[i],
                    salt: comp_salts[i].clone(),
                    path: comp_tree.prove(i),
                };
                QueryOpenings {
                    trace: open_trace(q),
                    trace_sibling: open_trace(q + half),
                    composition: open_comp(q),
                    composition_sibling: open_comp(q + half),
                }
            })
            .collect();

        let fri = E::fri_into_proof(fri_prover, &query_indices);

        StarkProof {
            log_trace_len: log_trace,
            trace_root,
            composition_root,
            ood_trace,
            ood_composition,
            fri,
            queries,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn compute_composition(
        &self,
        log_trace: u32,
        log_eval: u32,
        rotation_step: usize,
        points: &[E::Point],
        vn: &[E::Base],
        trace_ldes: &[Vec<E::Base>],
        beta: E::Ext,
    ) -> Vec<E::Ext> {
        let air = self.air;
        let domain_size = points.len();
        let num_transitions = air.num_transition_constraints();
        let offsets = air.mask_offsets();

        let inv_vn = batch_inverse(vn);
        let excluded = excluded_points::<E, A>(air, log_trace);
        let selector = E::selector_over_eval(&excluded, points);

        let rotations: Vec<Vec<usize>> = offsets
            .iter()
            .map(|&k| {
                (0..domain_size)
                    .map(|i| E::rotation_index(log_eval, i, k * rotation_step))
                    .collect()
            })
            .collect();
        let beta_pows = powers(beta, num_transitions);

        let mut composition: Vec<E::Ext> = (0..domain_size)
            .map(|i| {
                let mask: Vec<Vec<E::Base>> = rotations
                    .iter()
                    .map(|rot| trace_ldes.iter().map(|c| c[rot[i]]).collect())
                    .collect();
                let mut c_out = vec![E::Base::ZERO; num_transitions];
                air.eval_transitions(&mask, &mut c_out);
                let scale = selector.as_ref().map_or(E::Base::ONE, |s| s[i]) * inv_vn[i];
                let mut acc = E::Ext::ZERO;
                for (t, &c) in c_out.iter().enumerate() {
                    acc = acc + beta_pows[t].mul_base(c * scale);
                }
                acc
            })
            .collect();

        let mut beta_pow = beta.pow(num_transitions as u64);
        for b in air.boundary_constraints(1usize << log_trace) {
            let bpoint = E::trace_point(log_trace, b.row);
            let inv_denom = E::boundary_denom_inv_over_eval(bpoint, points);
            let col = &trace_ldes[b.column];
            let value = E::Ext::from(b.value);
            for i in 0..domain_size {
                composition[i] = composition[i]
                    + beta_pow * (E::Ext::from(col[i]) - value) * inv_denom[i];
            }
            beta_pow *= beta;
        }
        composition
    }

    #[allow(clippy::too_many_arguments)]
    fn compute_deep(
        &self,
        points: &[E::Point],
        trace_ldes: &[Vec<E::Base>],
        composition: &[E::Ext],
        mask_points: &[E::Ood],
        ood_trace: &[Vec<E::Ext>],
        ood_composition: E::Ext,
        mu: E::Ext,
    ) -> Vec<E::Ext> {
        let domain_size = points.len();
        let num_cols = trace_ldes.len();
        let inv_vz: Vec<Vec<E::Ext>> = mask_points
            .iter()
            .map(|&z| E::deep_denom_inv_over_eval(z, points))
            .collect();
        let mu_pows = powers(mu, mask_points.len() * num_cols + 1);

        (0..domain_size)
            .map(|i| {
                let mut acc = E::Ext::ZERO;
                let mut term = 0;
                for (k, ood_row) in ood_trace.iter().enumerate() {
                    for (c, &ood) in ood_row.iter().enumerate() {
                        let diff = E::Ext::from(trace_ldes[c][i]) - ood;
                        acc = acc + mu_pows[term] * diff * inv_vz[k][i];
                        term += 1;
                    }
                }
                acc + mu_pows[term] * (composition[i] - ood_composition) * inv_vz[0][i]
            })
            .collect()
    }
}

/// The excluded trace points for transition constraints.
fn excluded_points<E: Engine, A: Air<E::Base>>(air: &A, log_trace: u32) -> Vec<E::Point> {
    let n = 1usize << log_trace;
    (0..air.num_excluded_rows())
        .map(|k| E::trace_point(log_trace, n - 1 - k))
        .collect()
}

/// Evaluate the composition polynomial `F` at an out-of-domain point from the
/// claimed mask values — the verifier's identity, and the prover's `F(gamma)`.
pub(crate) fn eval_composition_at<E: Engine, A: Air<E::Base>>(
    air: &A,
    log_trace: u32,
    beta: E::Ext,
    ood: E::Ood,
    mask: &[Vec<E::Ext>],
) -> E::Ext {
    let num_transitions = air.num_transition_constraints();
    let mut c_out = vec![E::Ext::ZERO; num_transitions];
    air.eval_transitions(mask, &mut c_out);

    let excluded = excluded_points::<E, A>(air, log_trace);
    let selector = if excluded.is_empty() {
        E::Ext::ONE
    } else {
        E::selector_at_ood(&excluded, ood)
    };
    let inv_vh = E::vanishing_at_ood(log_trace, ood).inverse();

    let mut acc = E::Ext::ZERO;
    let mut beta_pow = E::Ext::ONE;
    for c in c_out {
        acc = acc + beta_pow * c * selector * inv_vh;
        beta_pow *= beta;
    }
    for b in air.boundary_constraints(1usize << log_trace) {
        let bpoint = E::trace_point(log_trace, b.row);
        let quotient = (mask[0][b.column] - E::Ext::from(b.value))
            * E::boundary_denom_inv_at_ood(bpoint, ood);
        acc = acc + beta_pow * quotient;
        beta_pow *= beta;
    }
    acc
}

pub(crate) fn init_transcript<E: Engine, A: Air<E::Base>>(
    air: &A,
    options: &ProofOptions,
    log_trace: u32,
    log_eval: u32,
) -> Transcript {
    let mut t = Transcript::new(air.label().as_bytes());
    let mut params = Vec::new();
    params.extend_from_slice(&(log_trace as u64).to_le_bytes());
    params.extend_from_slice(&(log_eval as u64).to_le_bytes());
    params.extend_from_slice(&(air.num_columns() as u64).to_le_bytes());
    params.extend_from_slice(&(air.num_transition_constraints() as u64).to_le_bytes());
    params.extend_from_slice(&(options.num_queries as u64).to_le_bytes());
    params.push(options.zk as u8);
    params.extend_from_slice(E::label().as_bytes());
    for b in air.boundary_constraints(1usize << log_trace) {
        params.extend_from_slice(&(b.column as u64).to_le_bytes());
        params.extend_from_slice(&(b.row as u64).to_le_bytes());
        params.extend_from_slice(&b.value.to_bytes());
    }
    t.absorb(&params);
    t
}

pub(crate) fn powers<F: Field>(base: F, count: usize) -> Vec<F> {
    let mut out = Vec::with_capacity(count);
    let mut acc = F::ONE;
    for _ in 0..count {
        out.push(acc);
        acc *= base;
    }
    out
}

fn salts(count: usize, zk: bool, rng: &mut impl Rng) -> Vec<Vec<u8>> {
    if zk {
        (0..count)
            .map(|_| rng.r#gen::<[u8; SALT_LEN]>().to_vec())
            .collect()
    } else {
        vec![Vec::new(); count]
    }
}

pub(crate) fn trace_leaf_bytes<E: Engine>(values: &[E::Base], salt: &[u8]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(values.len() * E::Base::NUM_BYTES + salt.len());
    for v in values {
        bytes.extend_from_slice(&v.to_bytes());
    }
    bytes.extend_from_slice(salt);
    bytes
}

pub(crate) fn composition_leaf_bytes<E: Engine>(value: E::Ext, salt: &[u8]) -> Vec<u8> {
    let mut bytes = value.to_bytes();
    bytes.extend_from_slice(salt);
    bytes
}
