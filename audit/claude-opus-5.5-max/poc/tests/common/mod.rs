//! Tree builders and openings mirroring the private helpers in `fibonacci.rs`.

#![allow(dead_code)]

use rand::Rng;
use toyni::babybear::BabyBear;
use toyni::ext::Ext;
use toyni::fibonacci::{MerkleOpening, FriPairOpening, QueryProof, StarkProof, BLOWUP, MASK_DEGREE, NUM_QUERIES};
use toyni::math::fri::fri_fold_ext;
use toyni::merkle::MerkleTree;
use toyni::transcript::FiatShamirTranscript;

pub fn salted(evals: &[BabyBear], rng: &mut impl Rng) -> (MerkleTree, Vec<Vec<u8>>) {
    let salts: Vec<Vec<u8>> = (0..evals.len()).map(|_| rng.r#gen::<[u8; 16]>().to_vec()).collect();
    let leaves = evals.iter().zip(&salts).map(|(v, s)| [s.as_slice(), &v.to_bytes()[..]].concat()).collect();
    (MerkleTree::new(leaves), salts)
}

pub fn ext_tree(evals: &[Ext]) -> MerkleTree {
    MerkleTree::new(evals.iter().map(|v| v.to_bytes().to_vec()).collect())
}

pub fn pair_tree(layer: &[Ext]) -> MerkleTree {
    let half = layer.len() / 2;
    MerkleTree::new((0..half).map(|j| [layer[j].to_bytes(), layer[j + half].to_bytes()].concat()).collect())
}

pub fn open(tree: &MerkleTree, salts: &[Vec<u8>], evals: &[BabyBear], i: usize) -> MerkleOpening {
    MerkleOpening { index: i, value: evals[i], proof: tree.get_proof(i).unwrap(), salt: salts[i].clone() }
}

pub fn open_pair(tree: &MerkleTree, layer: &[Ext], lo: usize) -> FriPairOpening {
    FriPairOpening { index: lo, values: [layer[lo], layer[lo + layer.len() / 2]], proof: tree.get_proof(lo).unwrap() }
}

pub fn derive_z(ts: &mut FiatShamirTranscript) -> Ext {
    loop {
        let z = ts.squeeze_ext_challenge();
        if !z.is_base() {
            return z;
        }
    }
}

/// Commits `layer0` as the DEEP layer, runs FRI exactly as the honest prover
/// and assembles the proof around the given commitments and OOD claims.
#[allow(clippy::too_many_arguments)]
pub fn finish(
    mut ts: FiatShamirTranscript,
    trace_len: usize,
    xs: &[BabyBear],
    trace: (&MerkleTree, &[Vec<u8>], &[BabyBear]),
    quotient: (&MerkleTree, &[Vec<u8>], &[BabyBear]),
    ood: [Ext; 4],
    layer0: Vec<Ext>,
) -> StarkProof {
    let lde_size = trace_len * BLOWUP;
    let mut layers = vec![layer0];
    let mut trees = Vec::new();
    let mut fri_commitments: Vec<Vec<u8>> = Vec::new();
    let mut fold_xs = xs.to_vec();
    let final_size = lde_size / (trace_len + MASK_DEGREE).next_power_of_two();
    loop {
        let layer = layers.last().unwrap();
        let tree = if layer.len() > final_size { pair_tree(layer) } else { ext_tree(layer) };
        let root = tree.root().unwrap();
        ts.absorb_commitment(&root);
        fri_commitments.push(root);
        trees.push(tree);
        if layers.last().unwrap().len() <= final_size {
            break;
        }
        let beta = ts.squeeze_ext_challenge();
        let folded = fri_fold_ext(layers.last().unwrap(), &fold_xs, beta);
        fold_xs.truncate(folded.len());
        for x in &mut fold_xs {
            *x = *x * *x;
        }
        layers.push(folded);
    }
    let (t_tree, t_salts, t_vals) = trace;
    let (q_tree, q_salts, q_vals) = quotient;
    let query_proofs = ts
        .squeeze_indices(NUM_QUERIES, lde_size)
        .into_iter()
        .map(|qi| {
            let mut lo = qi;
            let fri_pairs = (0..layers.len() - 1)
                .map(|k| {
                    lo %= layers[k].len() / 2;
                    open_pair(&trees[k], &layers[k], lo)
                })
                .collect();
            QueryProof {
                index: qi,
                trace_opening: open(t_tree, t_salts, t_vals, qi),
                trace_opening_g: open(t_tree, t_salts, t_vals, (qi + BLOWUP) % lde_size),
                trace_opening_gg: open(t_tree, t_salts, t_vals, (qi + 2 * BLOWUP) % lde_size),
                quotient_opening: open(q_tree, q_salts, q_vals, qi),
                fri_pairs,
            }
        })
        .collect();
    StarkProof {
        trace_len,
        lde_size,
        trace_commitment: t_tree.root().unwrap(),
        quotient_commitment: q_tree.root().unwrap(),
        t_z: ood[0],
        t_gz: ood[1],
        t_ggz: ood[2],
        q_z: ood[3],
        fri_commitments,
        fri_final_layer: layers.last().unwrap().clone(),
        query_proofs,
    }
}
