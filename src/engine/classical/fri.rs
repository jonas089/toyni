//! Classical (2-adic) FRI over a multiplicative coset, no dimension-gap step.
//!
//! Folds a Reed-Solomon codeword `f` over the evaluation coset down to a
//! constant, one octave per round: `g(x^2) = (f(x)+f(-x))/2 +
//! beta*(f(x)-f(-x))/(2x)`, with `-x` the sibling at index `i + |D|/2`. The
//! final layer is sent in the clear and must be constant, which enforces the
//! degree bound.

use crate::engine::classical::domain::coset_elements;
use crate::field::{batch_inverse, BabyBear, BabyBearExt, Field};
use crate::merkle::{verify_path, Hash, MerkleTree};
use crate::transcript::Transcript;

type Ext = BabyBearExt;

#[derive(Debug, Clone)]
pub struct FriLayerOpening {
    pub value: Ext,
    pub sibling: Ext,
    pub path: Vec<Hash>,
    pub sibling_path: Vec<Hash>,
}

#[derive(Debug, Clone)]
pub struct FriQueryProof {
    pub layers: Vec<FriLayerOpening>,
}

#[derive(Debug, Clone)]
pub struct FriProof {
    pub layer_roots: Vec<Hash>,
    pub final_layer: Vec<Ext>,
    pub query_proofs: Vec<FriQueryProof>,
}

pub struct FriProver {
    layers: Vec<Vec<Ext>>,
    trees: Vec<MerkleTree>,
    pub layer_roots: Vec<Hash>,
    pub final_layer: Vec<Ext>,
}

fn ext_leaves(values: &[Ext]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len() * BabyBearExt::NUM_BYTES);
    for v in values {
        out.extend_from_slice(&v.to_bytes());
    }
    out
}

fn commit_layer(values: &[Ext]) -> MerkleTree {
    MerkleTree::from_packed(&ext_leaves(values), BabyBearExt::NUM_BYTES)
}

/// One fold: `g[i] = (a+b)/2 + beta*(a-b)/(2*x_i)`, with `(a,b)` at `i` and
/// `i+half`, `x_i` the layer's coset element.
fn fold_layer(values: &[Ext], xs_inv: &[BabyBear], beta: Ext) -> Vec<Ext> {
    let half = values.len() / 2;
    let half_inv = BabyBear::from_u64(2).inverse();
    (0..half)
        .map(|i| {
            let a = values[i];
            let b = values[i + half];
            let avg = (a + b).mul_base(half_inv);
            let diff = (a - b).mul_base(half_inv * xs_inv[i]);
            avg + beta * diff
        })
        .collect()
}

/// Commit phase: fold `word` `log_bound` times to a constant final layer.
pub fn commit(
    transcript: &mut Transcript,
    word: Vec<Ext>,
    log_eval: u32,
    log_bound: u32,
) -> FriProver {
    let size = 1usize << log_eval;
    assert_eq!(word.len(), size);
    let mut xs = coset_elements(log_eval);
    let mut layers = vec![word];
    let mut trees = Vec::new();
    let mut layer_roots = Vec::new();

    for fold in 0..log_bound {
        let cur = layers.last().unwrap();
        let half = cur.len() / 2;
        let xs_inv = batch_inverse(&xs[..half]);
        let beta = super::draw_ext(transcript);
        let folded = fold_layer(cur, &xs_inv, beta);
        // Square the coset for the next layer.
        xs.truncate(half);
        for x in &mut xs {
            *x = *x * *x;
        }
        if fold + 1 < log_bound {
            let tree = commit_layer(&folded);
            let root = tree.root();
            transcript.absorb_commitment(&root);
            layer_roots.push(root);
            trees.push(tree);
            layers.push(folded);
        } else {
            for v in &folded {
                transcript.absorb_field(*v);
            }
            return FriProver {
                layers,
                trees,
                layer_roots,
                final_layer: folded,
            };
        }
    }
    unreachable!("log_bound >= 1")
}

impl FriProver {
    pub fn into_proof(self, query_indices: &[usize]) -> FriProof {
        let query_proofs = query_indices
            .iter()
            .map(|&index| {
                let mut openings = Vec::new();
                let mut pos = index;
                for (layer, tree) in self.layers[1..].iter().zip(self.trees.iter()) {
                    let half = layer.len() / 2;
                    pos %= half;
                    openings.push(FriLayerOpening {
                        value: layer[pos],
                        sibling: layer[pos + half],
                        path: tree.prove(pos),
                        sibling_path: tree.prove(pos + half),
                    });
                }
                FriQueryProof { layers: openings }
            })
            .collect();
        FriProof {
            layer_roots: self.layer_roots,
            final_layer: self.final_layer,
            query_proofs,
        }
    }
}

pub fn replay(transcript: &mut Transcript, proof: &FriProof, log_bound: u32) -> Vec<Ext> {
    let mut challenges = Vec::with_capacity(log_bound as usize);
    for fold in 0..log_bound {
        challenges.push(super::draw_ext(transcript));
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

pub fn check_shape(proof: &FriProof, log_eval: u32, log_bound: u32, num_queries: usize) -> bool {
    if proof.layer_roots.len() != (log_bound as usize).saturating_sub(1) {
        return false;
    }
    let final_size = (1usize << log_eval) >> log_bound;
    if proof.final_layer.len() != final_size {
        return false;
    }
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

/// Verify one query against the fold chain. `word_q`, `word_sib` are the
/// DEEP-batched values at `index` and `index + |D|/2` (layer-0 values),
/// recomputed by the verifier from its own openings.
#[allow(clippy::too_many_arguments)]
pub fn verify_query(
    proof: &FriProof,
    log_eval: u32,
    challenges: &[Ext],
    index: usize,
    word_q: Ext,
    word_sib: Ext,
    qp: &FriQueryProof,
) -> bool {
    // Coset element `7 * omega^i`, computed lazily (materializing the whole
    // domain per query would make verification O(queries * |D|)).
    let shift = BabyBear::from_u64(super::domain::COSET_SHIFT);
    let omega = BabyBear::get_root_of_unity(log_eval);
    let coset_at = |i: usize| shift * omega.pow(i as u64);
    let half_inv = BabyBear::from_u64(2).inverse();

    // Fold 0: layer 0 → layer 1, x-twiddle at the query's coset element.
    let x0_inv = coset_at(index).inverse();
    let avg = (word_q + word_sib).mul_base(half_inv);
    let diff = (word_q - word_sib).mul_base(half_inv * x0_inv);
    let mut expected = avg + challenges[0] * diff;

    let mut pos = index;
    let mut layer_size = (1usize << log_eval) / 2;

    for (fold, opening) in qp.layers.iter().enumerate() {
        let half = layer_size / 2;
        pos %= half;
        let root = &proof.layer_roots[fold];
        if !verify_path(root, pos, &opening.value.to_bytes(), &opening.path) {
            return false;
        }
        if !verify_path(root, pos + half, &opening.sibling.to_bytes(), &opening.sibling_path) {
            return false;
        }
        let came_low = (index % layer_size) < half;
        let at_pos = if came_low { opening.value } else { opening.sibling };
        if at_pos != expected {
            return false;
        }
        let x = coset_at(pos).pow(1u64 << (fold + 1));
        let a = opening.value;
        let b = opening.sibling;
        let avg = (a + b).mul_base(half_inv);
        let diff = (a - b).mul_base(half_inv * x.inverse());
        expected = avg + challenges[fold + 1] * diff;
        layer_size = half;
    }

    let final_half = proof.final_layer.len();
    proof.final_layer[pos % final_half] == expected
}
