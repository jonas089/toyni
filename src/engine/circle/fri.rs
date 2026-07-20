//! Circle FRI: the low-degree test over the circle (Section 6, Protocol 1).
//!
//! Proves that a word `f ∈ QM31^D` over a standard position coset `D` of size
//! `2^m` is close to the circle code `C_{2^b}(QM31, D)`, i.e. agrees with a
//! polynomial from `L_{2^b}` on most of `D`.
//!
//! Structure:
//! 1. **Decomposition** (the dimension gap): `L_{2^b} = L'_{2^b} ⊕ <v_b>`,
//!    so the prover extracts `λ = <f, v_b>_D / <v_b, v_b>_D` (Lemma 7:
//!    `v_b ⊥ L'` over any `G`- and `J`-invariant domain) and sends it; the
//!    test proceeds on `g = f - λ·v_b ∈ L'_{2^b}`.
//! 2. **Folding**: the first fold is along the involution `J` with twiddle
//!    `y`; each later fold is along `π(x) = 2x^2 - 1` with twiddle `x`.
//!    `b` folds take `deg < 2^(b-1)` down to constants.
//! 3. **Query phase**: spot checks of every fold at random positions.
//!
//! The final (constant, for an honest prover) layer is sent in the clear.

use crate::engine::circle::geometry::{vanishing_eval, CircleDomain};
use crate::field::{M31, QM31};
use crate::merkle::{verify_path, Hash, MerkleTree};
use crate::transcript::Transcript;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// One committed FRI layer opened at a query: the value at the (canonical)
/// position and at its fold sibling, each with a Merkle path.
#[derive(Debug, Clone)]
pub struct FriLayerOpening {
    pub value: QM31,
    pub sibling: QM31,
    pub path: Vec<Hash>,
    pub sibling_path: Vec<Hash>,
}

/// Per-query FRI openings across all committed layers.
#[derive(Debug, Clone)]
pub struct FriQueryProof {
    pub layers: Vec<FriLayerOpening>,
}

#[derive(Debug, Clone)]
pub struct FriProof {
    /// Decomposition coefficient of the dimension gap.
    pub lambda: QM31,
    /// Roots of the committed folded layers (folds `1..b`, excluding final).
    pub layer_roots: Vec<Hash>,
    /// The final layer in the clear; must be constant.
    pub final_layer: Vec<QM31>,
    pub query_proofs: Vec<FriQueryProof>,
}

/// Prover state kept between the commit and query phases.
pub struct FriProver {
    layers: Vec<Vec<QM31>>,
    trees: Vec<MerkleTree>,
    pub lambda: QM31,
    pub layer_roots: Vec<Hash>,
    pub final_layer: Vec<QM31>,
}

/// Fold one layer: `g'[i] = (a + b)/2 + λ · (a - b)/2t` with `(a, b)` the
/// values at `i` and `i + half`, and `t` the layer twiddle at `i`
/// (y-coordinate for the first fold, x-coordinate iterates afterwards).
pub(crate) fn fold_layer_cpu(values: &[QM31], inv_twiddles: &[M31], lambda: QM31) -> Vec<QM31> {
    let half = values.len() / 2;
    assert_eq!(inv_twiddles.len(), half);
    let half_scalar = M31::half();
    let fold = |i: usize| {
        let a = values[i];
        let b = values[i + half];
        let avg = (a + b).mul_base(half_scalar);
        let diff = (a - b).mul_base(half_scalar * inv_twiddles[i]);
        avg + lambda * diff
    };
    #[cfg(feature = "parallel")]
    if half >= 1 << 13 {
        return (0..half).into_par_iter().map(fold).collect();
    }
    (0..half).map(fold).collect()
}

/// Compute the decomposition coefficient λ of `f = g + λ·v_b` over the
/// domain, using orthogonality of `v_b` to the FFT space (Lemma 7).
fn decompose_lambda(evals: &[QM31], vb: &[M31]) -> QM31 {
    let mut num = QM31::ZERO;
    let mut den = M31::ZERO;
    for (f, v) in evals.iter().zip(vb.iter()) {
        num += f.mul_base(*v);
        den += *v * *v;
    }
    assert!(!den.is_zero(), "degenerate FRI domain: <v_b, v_b> = 0");
    num.mul_base(den.inverse())
}

/// Evaluations of `v_b` (`b = log_bound`) over the domain, in domain order.
fn vanishing_on_domain(domain: &CircleDomain, log_bound: u32) -> Vec<M31> {
    let half: Vec<M31> = domain
        .half_coset
        .points()
        .iter()
        .map(|p| vanishing_eval(log_bound, p.x))
        .collect();
    // v_b depends only on x, and J preserves x: second half mirrors the first.
    let mut out = half.clone();
    out.extend_from_slice(&half);
    out
}

/// Commit phase. `evals` is the word over `domain` (size `2^m`); the test
/// enforces proximity to `L_{2^b}` with `b = log_bound`, folding `b` times.
///
/// `inv_twiddles` are the CFFT inverse twiddles of `domain` (shared with the
/// LDE machinery).
pub fn fri_commit(
    transcript: &mut Transcript,
    evals: Vec<QM31>,
    domain: &CircleDomain,
    inv_twiddles: &[Vec<M31>],
    log_bound: u32,
) -> FriProver {
    let m = domain.log_size();
    assert_eq!(evals.len(), domain.size());
    assert!(log_bound >= 1 && log_bound < m);

    // 1. Dimension-gap decomposition.
    let vb = vanishing_on_domain(domain, log_bound);
    let lambda = decompose_lambda(&evals, &vb);
    let g: Vec<QM31> = evals
        .iter()
        .zip(vb.iter())
        .map(|(f, v)| *f - lambda.mul_base(*v))
        .collect();
    transcript.absorb_field(lambda);

    // 2. Folding: b folds; commit every intermediate layer, send the last in
    //    the clear.
    let mut layers: Vec<Vec<QM31>> = vec![g];
    let mut trees: Vec<MerkleTree> = Vec::new();
    let mut layer_roots: Vec<Hash> = Vec::new();

    for fold in 0..log_bound {
        let cur = layers.last().unwrap();
        let challenge = super::draw_ext_shared(transcript);
        let folded = crate::engine::circle::backend::fri_fold(cur, &inv_twiddles[fold as usize], challenge);
        if fold + 1 < log_bound {
            let tree = crate::engine::circle::backend::merkle_qm31(&folded);
            let root = tree.root();
            transcript.absorb_commitment(&root);
            layer_roots.push(root);
            trees.push(tree);
            layers.push(folded);
        } else {
            // Final layer: sent in the clear, absorbed wholesale.
            for v in &folded {
                transcript.absorb_field(*v);
            }
            return FriProver {
                layers,
                trees,
                lambda,
                layer_roots,
                final_layer: folded,
            };
        }
    }
    unreachable!("log_bound >= 1")
}

impl FriProver {
    /// Open all committed layers at the positions a query at `index`
    /// (canonical, `< |D|/2`) traces through the fold chain.
    pub fn open_query(&self, index: usize) -> FriQueryProof {
        let mut openings = Vec::new();
        // layers[0] is the input layer (not committed here; the STARK derives
        // it from the trace openings). Committed layers start at fold 1.
        let mut pos = index;
        for (layer, tree) in self.layers[1..].iter().zip(self.trees.iter()) {
            let half = layer.len() / 2;
            pos %= half; // canonical position in this layer's pairing
            openings.push(FriLayerOpening {
                value: layer[pos],
                sibling: layer[pos + half],
                path: tree.prove(pos),
                sibling_path: tree.prove(pos + half),
            });
        }
        FriQueryProof { layers: openings }
    }

    pub fn into_proof(self, query_indices: &[usize]) -> FriProof {
        let query_proofs = query_indices.iter().map(|&q| self.open_query(q)).collect();
        FriProof {
            lambda: self.lambda,
            layer_roots: self.layer_roots,
            final_layer: self.final_layer,
            query_proofs,
        }
    }
}

/// Replay of the FRI transcript on the verifier side: absorb λ and the layer
/// roots, drawing the fold challenges in the same order as the prover.
pub fn fri_replay_transcript(
    transcript: &mut Transcript,
    proof: &FriProof,
    log_bound: u32,
) -> Vec<QM31> {
    transcript.absorb_field(proof.lambda);
    let mut challenges = Vec::with_capacity(log_bound as usize);
    for fold in 0..log_bound {
        challenges.push(super::draw_ext_shared(transcript));
        if fold + 1 < log_bound {
            transcript.absorb_commitment(&proof.layer_roots[fold as usize]);
        } else {
            for v in &proof.final_layer {
                transcript.absorb_field(*v);
            }
        }
    }
    challenges
}

/// Verify the shape and final layer of a FRI proof (everything that is not
/// per-query). Returns false on malformed proofs.
pub fn fri_check_shape(proof: &FriProof, domain: &CircleDomain, log_bound: u32, num_queries: usize) -> bool {
    if proof.layer_roots.len() != (log_bound as usize).saturating_sub(1) {
        return false;
    }
    let final_size = domain.size() >> log_bound;
    if proof.final_layer.len() != final_size {
        return false;
    }
    // Degree-0 (constant) final layer enforces the overall degree bound.
    if !proof.final_layer.iter().all(|v| *v == proof.final_layer[0]) {
        return false;
    }
    if proof.query_proofs.len() != num_queries {
        return false;
    }
    proof
        .query_proofs
        .iter()
        .all(|qp| qp.layers.len() == (log_bound as usize).saturating_sub(1))
}

/// Verify one query against the fold chain.
///
/// `index` is the canonical query position (`< |D|/2`); `value` and
/// `sibling_value` are the input-layer values at `index` and `index + |D|/2`
/// **after** the caller has removed the λ·v_b component (the caller computes
/// these from its own openings — this binds FRI to the outer commitments).
pub fn fri_verify_query(
    proof: &FriProof,
    domain: &CircleDomain,
    challenges: &[QM31],
    index: usize,
    value: QM31,
    sibling_value: QM31,
    query_proof: &FriQueryProof,
) -> bool {
    let half_scalar = M31::half();

    // Fold 0 (J-fold): twiddle is the y-coordinate of the domain point.
    let p = domain.at(index);
    let avg = (value + sibling_value).mul_base(half_scalar);
    let diff = (value - sibling_value).mul_base(half_scalar * p.y.inverse());
    let mut expected = avg + challenges[0] * diff;

    // x-coordinate for subsequent folds: repeated application of 2x^2 - 1 to
    // the x-coordinate of the *canonical pair index* of each layer.
    let mut pos = index;
    let mut layer_size = domain.size() / 2;

    for (fold, opening) in query_proof.layers.iter().enumerate() {
        let half = layer_size / 2;
        pos %= half;

        // Merkle checks against the committed layer.
        let root = &proof.layer_roots[fold];
        if !verify_path(root, pos, &opening.value.to_bytes(), &opening.path) {
            return false;
        }
        if !verify_path(root, pos + half, &opening.sibling.to_bytes(), &opening.sibling_path) {
            return false;
        }

        // Fold-consistency: the previous fold must land on this layer. The
        // expected value sits at `pos` or `pos + half` depending on which
        // side the query path came down.
        let came_low = (index % layer_size) < half;
        let at_pos = if came_low { opening.value } else { opening.sibling };
        if at_pos != expected {
            return false;
        }

        // Twiddle: x-coordinate of the canonical index in this layer, squared
        // `fold` times under π_x.
        let mut x = domain.at(pos).x;
        for _ in 0..fold {
            x = M31::TWO * x * x - M31::ONE;
        }

        let a = opening.value;
        let b = opening.sibling;
        let avg = (a + b).mul_base(half_scalar);
        let diff = (a - b).mul_base(half_scalar * x.inverse());
        expected = avg + challenges[fold + 1] * diff;

        layer_size = half;
    }

    // Land on the final layer.
    let final_half = proof.final_layer.len();
    let final_pos = pos % final_half;
    // After the last fold the expected value must equal the final layer entry
    // at the folded position.
    proof.final_layer[final_pos] == expected
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::circle::cfft::{embed_coeffs, evaluate, Twiddles};
    use rand::SeedableRng;

    /// Build the evaluations over `domain` of a random element of L_{2^b}:
    /// a random FFT-space element (order b) plus a random multiple of v_b.
    fn random_codeword(
        domain: &CircleDomain,
        tw: &Twiddles,
        log_bound: u32,
        seed: u64,
    ) -> (Vec<QM31>, QM31) {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let coeffs: Vec<QM31> = (0..1usize << log_bound).map(|_| QM31::random(&mut rng)).collect();
        let mut lde = embed_coeffs(&coeffs, log_bound, domain.log_size());
        evaluate(&mut lde, tw);
        let mu = QM31::random(&mut rng);
        let vb = vanishing_on_domain(domain, log_bound);
        let evals: Vec<QM31> = lde
            .iter()
            .zip(vb.iter())
            .map(|(g, v)| *g + mu.mul_base(*v))
            .collect();
        (evals, mu)
    }

    fn run_fri(evals: Vec<QM31>, domain: &CircleDomain, tw: &Twiddles, log_bound: u32) -> (FriProof, bool) {
        let mut t = Transcript::new(b"fri-test");
        let prover = fri_commit(&mut t, evals.clone(), domain, &tw.inverse, log_bound);
        let queries = t.draw_indices(8, domain.size() / 2);
        let proof = prover.into_proof(&queries);

        // Verify.
        let mut t = Transcript::new(b"fri-test");
        let challenges = fri_replay_transcript(&mut t, &proof, log_bound);
        let queries_v = t.draw_indices(8, domain.size() / 2);
        assert_eq!(queries, queries_v);
        if !fri_check_shape(&proof, domain, log_bound, 8) {
            return (proof, false);
        }
        let vb = vanishing_on_domain(domain, log_bound);
        let ok = queries_v.iter().enumerate().all(|(qi, &q)| {
            let sib = domain.conjugate_index(q);
            let val = evals[q] - proof.lambda.mul_base(vb[q]);
            let sval = evals[sib] - proof.lambda.mul_base(vb[sib]);
            fri_verify_query(
                &proof, domain, &challenges, q, val, sval,
                &proof.query_proofs[qi],
            )
        });
        (proof, ok)
    }

    #[test]
    fn honest_codeword_accepted_and_lambda_recovered() {
        let log_domain = 9u32;
        let log_bound = 5u32;
        let domain = CircleDomain::standard(log_domain);
        let tw = Twiddles::new(&domain);
        let (evals, mu) = random_codeword(&domain, &tw, log_bound, 1);
        let (proof, ok) = run_fri(evals, &domain, &tw, log_bound);
        assert!(ok, "honest proof rejected");
        assert_eq!(proof.lambda, mu, "decomposition lambda mismatch");
        // Final layer is constant.
        assert!(proof.final_layer.iter().all(|v| *v == proof.final_layer[0]));
    }

    #[test]
    fn tampered_final_layer_rejected() {
        let domain = CircleDomain::standard(8);
        let tw = Twiddles::new(&domain);
        let (evals, _) = random_codeword(&domain, &tw, 4, 2);
        let mut t = Transcript::new(b"fri-test");
        let prover = fri_commit(&mut t, evals, &domain, &tw.inverse, 4);
        let queries = t.draw_indices(8, domain.size() / 2);
        let mut proof = prover.into_proof(&queries);
        proof.final_layer[0] += QM31::ONE;
        assert!(!fri_check_shape(&proof, &domain, 4, 8));
    }

    #[test]
    fn high_degree_word_rejected() {
        // A random word (degree ~ |D|) should fail: the final layer of an
        // honest fold of it will not be constant.
        let mut rng = rand::rngs::StdRng::seed_from_u64(3);
        let domain = CircleDomain::standard(8);
        let tw = Twiddles::new(&domain);
        let evals: Vec<QM31> = (0..domain.size()).map(|_| QM31::random(&mut rng)).collect();
        let mut t = Transcript::new(b"fri-test");
        let prover = fri_commit(&mut t, evals, &domain, &tw.inverse, 4);
        assert!(
            !prover.final_layer.iter().all(|v| *v == prover.final_layer[0]),
            "random word folded to a constant — vanishingly unlikely"
        );
    }

    #[test]
    fn wrong_input_value_rejected() {
        let domain = CircleDomain::standard(8);
        let tw = Twiddles::new(&domain);
        let (evals, _) = random_codeword(&domain, &tw, 4, 4);

        let mut t = Transcript::new(b"fri-test");
        let prover = fri_commit(&mut t, evals.clone(), &domain, &tw.inverse, 4);
        let queries = t.draw_indices(4, domain.size() / 2);
        let proof = prover.into_proof(&queries);

        let mut t = Transcript::new(b"fri-test");
        let challenges = fri_replay_transcript(&mut t, &proof, 4);
        let queries_v = t.draw_indices(4, domain.size() / 2);
        let vb = vanishing_on_domain(&domain, 4);

        // Feed a corrupted first-layer value for query 0: must fail.
        let q = queries_v[0];
        let sib = domain.conjugate_index(q);
        let val = evals[q] - proof.lambda.mul_base(vb[q]) + QM31::ONE;
        let sval = evals[sib] - proof.lambda.mul_base(vb[sib]);
        assert!(!fri_verify_query(
            &proof, &domain, &challenges, q, val, sval, &proof.query_proofs[0]
        ));
    }
}
