use crate::babybear::BabyBear;
use crate::ext::Ext;
use crate::math::domain::BabyBearDomain;
use crate::math::fri::fri_fold_ext;
use crate::merkle::{MerkleProof, MerkleTree};
use crate::program::trace::ExecutionTrace;
use crate::transcript::FiatShamirTranscript;
use rand::Rng;

/// Spot-check queries. At rate 1/8 each gives ~3 bits, so 44 ~= 132 bits.
pub const NUM_QUERIES: usize = 44;
/// LDE blowup. Masking lifts deg(D) to ~4*trace_len, so blowup 32 keeps the
/// tested Reed-Solomon rate at 1/8.
pub const BLOWUP: usize = 32;
/// Coset shift used for the LDE domain.
pub const COSET_SHIFT: u64 = 7;
/// Random blinding coefficients per trace polynomial (T_hat = T + Z_H * R).
/// Covers every revealed trace evaluation: 3 openings per query + 3 OOD points.
pub const MASK_DEGREE: usize = 3 * NUM_QUERIES + 8;

/// Blind a base-field column in place: P += Z_H * R with fresh random R.
/// Z_H = x^n - 1, so Z_H*R = x^n*R - R (subtract R low, add it back shifted n).
/// Z_H vanishes on the trace domain, so the masked polynomial equals the
/// original there (constraints unchanged) but its off-domain openings are
/// uniformly random.
fn mask_poly_base(poly: &mut Vec<BabyBear>, n: usize, rng: &mut impl rand::Rng) {
    let r: Vec<BabyBear> = (0..MASK_DEGREE).map(|_| BabyBear::random(rng)).collect();
    if poly.len() < n + MASK_DEGREE {
        poly.resize(n + MASK_DEGREE, BabyBear::zero());
    }
    for i in 0..MASK_DEGREE {
        poly[i] = poly[i] - r[i];
        poly[n + i] = poly[n + i] + r[i];
    }
}

// ── proof data structures ──────────────────────────────────────────────

/// Opening of a base-field committed position (trace / quotient). The leaf is
/// committed as H(salt || value); the salt hides unopened sibling values, which
/// are small base-field elements and otherwise recoverable by brute force.
#[derive(Debug, Clone)]
pub struct MerkleOpening {
    pub index: usize,
    pub value: BabyBear,
    pub proof: MerkleProof,
    pub salt: Vec<u8>,
}

/// Opening of an extension-field committed position (DEEP / FRI layers). No
/// salt: Ext values are 32 bytes of high entropy, so they are not brute-forceable.
#[derive(Debug, Clone)]
pub struct MerkleOpeningExt {
    pub index: usize,
    pub value: Ext,
    pub proof: MerkleProof,
}

/// All the data the verifier needs to check one query position across
/// every FRI layer.
#[derive(Debug, Clone)]
pub struct QueryProof {
    /// The original query index into the first FRI layer (= DEEP evals).
    pub index: usize,

    /// Openings into the DEEP layer (position and its pair).
    pub deep_opening: MerkleOpeningExt,
    pub deep_opening_pair: MerkleOpeningExt,

    /// Trace polynomial openings on the LDE domain:
    ///   trace_opening    → T(x_i)     at position qi
    ///   trace_opening_g  → T(g·x_i)   at position (qi + BLOWUP) % lde_size
    ///   trace_opening_gg → T(g²·x_i)  at position (qi + 2·BLOWUP) % lde_size
    pub trace_opening: MerkleOpening,
    pub trace_opening_g: MerkleOpening,
    pub trace_opening_gg: MerkleOpening,

    /// Quotient polynomial opening at the query position on the LDE domain.
    pub quotient_opening: MerkleOpening,

    /// For each intermediate FRI layer: openings of a position and its pair.
    pub fri_openings: Vec<(MerkleOpeningExt, MerkleOpeningExt)>,
}

/// A complete STARK proof that the verifier can check.
///
/// The trace and quotient stay in the base field; the out-of-domain point `z`,
/// FRI folding challenges, and everything derived from them (OOD evaluations,
/// DEEP composition, FRI layers) live in the quartic extension for ~124-bit
/// soundness.
#[derive(Debug)]
pub struct StarkProof {
    // protocol parameters baked into the proof
    pub trace_len: usize,
    pub lde_size: usize,

    // commitments
    pub trace_commitment: Vec<u8>,
    pub quotient_commitment: Vec<u8>,

    // out-of-domain evaluations (extension field)
    pub t_z: Ext,
    pub t_gz: Ext,
    pub t_ggz: Ext,
    pub q_z: Ext,

    // FRI commitments (layer 0 = DEEP evals, then each folded layer)
    pub fri_commitments: Vec<Vec<u8>>,
    /// Full final FRI layer, sent in the clear; the verifier checks it is a
    /// constant codeword. This is what enforces the low-degree bound.
    pub fri_final_layer: Vec<Ext>,

    // per-query openings
    pub query_proofs: Vec<QueryProof>,
}

// ── prover ──────────────────────────────────────────────────────────────

pub struct StarkProver {
    trace: ExecutionTrace,
}

impl StarkProver {
    pub fn new(trace: ExecutionTrace) -> Self {
        Self { trace }
    }

    pub fn generate_proof(&self, use_gpu: bool) -> StarkProof {
        let trace_len = self.trace.trace.len();
        let domain = BabyBearDomain::new(trace_len).with_gpu(use_gpu);
        let lde_size = trace_len * BLOWUP;
        let extended_domain = BabyBearDomain::new(lde_size).with_gpu(use_gpu);
        let shift = BabyBear::new(COSET_SHIFT);
        let shifted_domain = extended_domain.get_coset(shift);

        let g = domain.group_gen();
        let shifted_elements = shifted_domain.elements();
        let mut rng = rand::thread_rng();

        // ── 1. trace polynomial (+ zero-knowledge masking) ─────────────
        // Interpolate the column off the trace domain (ifft), blind it, then
        // evaluate over the shifted LDE coset (fft).
        let column = self.trace.get_column(0);
        let mut trace_poly = domain.ifft(&column);
        mask_poly_base(&mut trace_poly, trace_len, &mut rng);
        let trace_lde = shifted_domain.fft(&trace_poly);

        let trace_tree = build_merkle_tree(&trace_lde, &mut rng);
        let trace_commitment = trace_tree.root().unwrap();

        // ── 2. constraint & quotient ───────────────────────────────────
        // C(x) evaluated over the LDE coset by indexing the trace LDE: g·x_i is
        // BLOWUP positions further along the domain (g = ω_n = ω_N^BLOWUP).
        let c_evals: Vec<BabyBear> = (0..lde_size)
            .map(|i| {
                let t_x = trace_lde[i];
                let t_gx = trace_lde[(i + BLOWUP) % lde_size];
                let t_ggx = trace_lde[(i + 2 * BLOWUP) % lde_size];
                let x = shifted_elements[i];
                fibonacci_constraint(t_ggx, t_gx, t_x)
                    * boundary_constraint_1(x, g, trace_len)
                    * boundary_constraint_2(x, g, trace_len)
            })
            .collect();

        // Q = C / Z_H. Z_H(x) = x^n - 1 on the trace domain; divide the
        // constraint evals by it directly on the coset (no re-interpolation).
        let n = trace_len as u64;
        let q_evals: Vec<BabyBear> = (0..lde_size)
            .map(|i| c_evals[i] / (shifted_elements[i].pow(n) - BabyBear::one()))
            .collect();
        // Coefficients of Q, for the out-of-domain evaluation q(z).
        let q_coeffs = shifted_domain.ifft(&q_evals);

        let quotient_tree = build_merkle_tree(&q_evals, &mut rng);
        let quotient_commitment = quotient_tree.root().unwrap();

        // ── 3. Fiat-Shamir: derive OOD point z (extension field) ───────
        let mut transcript = FiatShamirTranscript::new();
        transcript.absorb_commitment(&trace_commitment);
        transcript.absorb_commitment(&quotient_commitment);

        let z = derive_z(&mut transcript);
        let gz = z.mul_base(g);
        let ggz = z.mul_base(g).mul_base(g);

        // ── 4. OOD evaluations (base-coefficient polys at the Ext point) ─
        let t_z = z.eval_base_at_ext(&trace_poly);
        let t_gz = gz.eval_base_at_ext(&trace_poly);
        let t_ggz = ggz.eval_base_at_ext(&trace_poly);
        let q_z = z.eval_base_at_ext(&q_coeffs);

        // Sanity: the constraint relation C(z) = Q(z)·Z(z) holds at z.
        let vanishing = domain.vanishing_poly_coeffs();
        let c_z = eval_fibonacci_constraint_ext(t_ggz, t_gz, t_z)
            * eval_boundary_1_ext(z, g, trace_len)
            * eval_boundary_2_ext(z, g, trace_len);
        assert_eq!(
            c_z,
            q_z * z.eval_base_at_ext(&vanishing),
            "Constraint check at z failed"
        );

        // Feed OOD values into transcript.
        transcript.absorb_ext(t_z);
        transcript.absorb_ext(t_gz);
        transcript.absorb_ext(t_ggz);
        transcript.absorb_ext(q_z);

        // ── 5. DEEP polynomial (extension field) ───────────────────────
        // D(x) = Σ (P(x) - P(z)) / (x - z) over {Q, T, T∘g, T∘g²}. x is base,
        // z is Ext, so the terms are Ext-valued; base openings are lifted.
        let d_evals: Vec<Ext> = (0..lde_size)
            .map(|i| {
                let x = shifted_elements[i];
                let inv_x_z = (Ext::from(x) - z).inverse();
                let t_x = Ext::from(trace_lde[i]);
                let t_gx = Ext::from(trace_lde[(i + BLOWUP) % lde_size]);
                let t_ggx = Ext::from(trace_lde[(i + 2 * BLOWUP) % lde_size]);
                let q_x = Ext::from(q_evals[i]);
                (q_x - q_z) * inv_x_z
                    + (t_ggx - t_ggz) * inv_x_z
                    + (t_gx - t_gz) * inv_x_z
                    + (t_x - t_z) * inv_x_z
            })
            .collect();

        // ── 6. FRI folding with Merkle commits (extension field) ───────
        let mut fri_layers: Vec<Vec<Ext>> = Vec::new();
        let mut fri_trees: Vec<MerkleTree> = Vec::new();
        let mut fri_commitments: Vec<Vec<u8>> = Vec::new();

        // Layer 0 = DEEP evaluations.
        fri_layers.push(d_evals.clone());
        let tree0 = build_merkle_tree_ext(&d_evals);
        let root0 = tree0.root().unwrap();
        transcript.absorb_commitment(&root0);
        fri_commitments.push(root0);
        fri_trees.push(tree0);

        let mut current = d_evals;
        let mut xs: Vec<BabyBear> = shifted_elements.clone();

        // Fold down to the degree-bound layer (size lde/D_BOUND). For an honest
        // codeword that layer is constant; the verifier checks that, which is
        // what enforces the degree bound. The round count is fixed, not
        // data-dependent.
        let fri_degree_bound = (trace_len + MASK_DEGREE).next_power_of_two();
        let final_layer_size = lde_size / fri_degree_bound;
        while current.len() > final_layer_size {
            let beta = transcript.squeeze_ext_challenge();

            let folded = fri_fold_ext(&current, &xs, beta);

            // Square the x-coordinates for the next domain.
            xs.truncate(folded.len());
            for x in &mut xs {
                *x = *x * *x;
            }

            fri_layers.push(folded.clone());
            let tree = build_merkle_tree_ext(&folded);
            let root = tree.root().unwrap();
            transcript.absorb_commitment(&root);
            fri_commitments.push(root);
            fri_trees.push(tree);

            current = folded;
        }

        let fri_final_layer = current;

        // ── 7. Query phase ─────────────────────────────────────────────
        let first_layer_half = fri_layers[0].len() / 2;
        let query_indices = transcript.squeeze_indices(NUM_QUERIES, first_layer_half);

        let mut query_proofs = Vec::with_capacity(NUM_QUERIES);

        for &qi in &query_indices {
            // Trace openings: T(x_i), T(g·x_i), T(g²·x_i).
            let idx_g = (qi + BLOWUP) % lde_size;
            let idx_gg = (qi + 2 * BLOWUP) % lde_size;
            let trace_opening = open_merkle(&trace_tree, &trace_lde, qi);
            let trace_opening_g = open_merkle(&trace_tree, &trace_lde, idx_g);
            let trace_opening_gg = open_merkle(&trace_tree, &trace_lde, idx_gg);

            // Quotient opening: Q(x_i).
            let quotient_opening = open_merkle(&quotient_tree, &q_evals, qi);

            // DEEP layer (layer 0) openings: position qi and its pair qi + half.
            let half0 = fri_layers[0].len() / 2;
            let deep_opening = open_merkle_ext(&fri_trees[0], &fri_layers[0], qi);
            let deep_opening_pair = open_merkle_ext(&fri_trees[0], &fri_layers[0], qi + half0);

            // Intermediate FRI layers.
            let mut fri_openings = Vec::new();
            let mut idx = qi;
            for layer_idx in 1..fri_layers.len() - 1 {
                let half = fri_layers[layer_idx].len() / 2;
                idx %= half;
                let op = open_merkle_ext(&fri_trees[layer_idx], &fri_layers[layer_idx], idx);
                let op_pair =
                    open_merkle_ext(&fri_trees[layer_idx], &fri_layers[layer_idx], idx + half);
                fri_openings.push((op, op_pair));
            }

            query_proofs.push(QueryProof {
                index: qi,
                deep_opening,
                deep_opening_pair,
                trace_opening,
                trace_opening_g,
                trace_opening_gg,
                quotient_opening,
                fri_openings,
            });
        }

        StarkProof {
            trace_len,
            lde_size,
            trace_commitment,
            quotient_commitment,
            t_z,
            t_gz,
            t_ggz,
            q_z,
            fri_commitments,
            fri_final_layer,
            query_proofs,
        }
    }
}

// ── helpers ────────────────────────────────────────────────────────────

fn fibonacci_constraint(t2: BabyBear, t1: BabyBear, t0: BabyBear) -> BabyBear {
    t2 - (t1 + t0)
}

fn boundary_constraint_1(x: BabyBear, g: BabyBear, n: usize) -> BabyBear {
    x - g.pow((n - 1) as u64)
}

fn boundary_constraint_2(x: BabyBear, g: BabyBear, n: usize) -> BabyBear {
    x - g.pow((n - 2) as u64)
}

/// A Merkle tree plus the per-leaf salts used to build it.
struct SaltedTree {
    tree: MerkleTree,
    salts: Vec<Vec<u8>>,
}

impl SaltedTree {
    fn root(&self) -> Option<Vec<u8>> {
        self.tree.root()
    }
}

/// Build a hiding base-field Merkle tree: each leaf is H(salt || value) with a
/// fresh salt.
fn build_merkle_tree(evals: &[BabyBear], rng: &mut impl Rng) -> SaltedTree {
    let salts: Vec<Vec<u8>> = (0..evals.len())
        .map(|_| rng.r#gen::<[u8; 16]>().to_vec())
        .collect();
    let leaves: Vec<Vec<u8>> = evals
        .iter()
        .zip(&salts)
        .map(|(v, s)| [s.as_slice(), &v.to_bytes()].concat())
        .collect();
    SaltedTree {
        tree: MerkleTree::new(leaves),
        salts,
    }
}

/// Build an extension-field Merkle tree (leaf = 32-byte Ext encoding). No salt:
/// Ext values carry full entropy, and the final layer's root is recomputed
/// directly by the verifier.
fn build_merkle_tree_ext(evals: &[Ext]) -> MerkleTree {
    MerkleTree::new(evals.iter().map(|v| v.to_bytes().to_vec()).collect())
}

/// Open a salted base-field Merkle tree at a given index.
fn open_merkle(tree: &SaltedTree, evals: &[BabyBear], index: usize) -> MerkleOpening {
    let proof = tree.tree.get_proof(index).expect("Index out of bounds");
    MerkleOpening {
        index,
        value: evals[index],
        proof,
        salt: tree.salts[index].clone(),
    }
}

/// Open an extension-field Merkle tree at a given index.
fn open_merkle_ext(tree: &MerkleTree, evals: &[Ext], index: usize) -> MerkleOpeningExt {
    let proof = tree.get_proof(index).expect("Index out of bounds");
    MerkleOpeningExt {
        index,
        value: evals[index],
        proof,
    }
}

/// Derive the out-of-domain point in the extension field. A random extension
/// element is outside the base evaluation domain unless it happens to be a base
/// element (negligible probability), so we only reject the base case to
/// guarantee z, g·z, g²·z ∉ domain (keeping the DEEP denominators invertible).
fn derive_z(transcript: &mut FiatShamirTranscript) -> Ext {
    loop {
        let z = transcript.squeeze_ext_challenge();
        if !z.is_base() {
            return z;
        }
    }
}

// ── public constraint helpers (for the verifier) ───────────────────────
// Evaluated at the extension-field OOD point.

pub fn eval_fibonacci_constraint_ext(t2: Ext, t1: Ext, t0: Ext) -> Ext {
    t2 - (t1 + t0)
}

pub fn eval_boundary_1_ext(x: Ext, g: BabyBear, n: usize) -> Ext {
    x - Ext::from(g.pow((n - 1) as u64))
}

pub fn eval_boundary_2_ext(x: Ext, g: BabyBear, n: usize) -> Ext {
    x - Ext::from(g.pow((n - 2) as u64))
}

#[cfg(test)]
mod tests {
    use super::StarkProver;
    use crate::{babybear::BabyBear, program::trace::ExecutionTrace};

    #[test]
    fn test_fibonacci() {
        let mut execution_trace = ExecutionTrace::new();
        let trace: Vec<u64> = fibonacci_list(64);
        let trace_field: Vec<BabyBear> = trace.iter().map(|x| BabyBear::new(*x)).collect();
        execution_trace.insert_column(trace_field);
        let stark = StarkProver::new(execution_trace.clone());
        let _proof = stark.generate_proof(false);
    }

    #[test]
    #[should_panic]
    fn test_invalid_trace_should_fail() {
        let mut execution_trace = ExecutionTrace::new();
        let mut trace: Vec<u64> = fibonacci_list(64);
        for i in 1..50 {
            trace[i] = i as u64 * 3143;
        }
        let trace_field: Vec<BabyBear> = trace.iter().map(|x| BabyBear::new(*x)).collect();
        execution_trace.insert_column(trace_field);
        let stark = StarkProver::new(execution_trace.clone());
        let _proof = stark.generate_proof(false);
    }

    fn fibonacci_list(n: usize) -> Vec<u64> {
        let mut fibs: Vec<u64> = Vec::with_capacity(n);
        let mut a = 1u64;
        let mut b = 1u64;
        for _ in 0..n {
            fibs.push(a);
            let next = a.wrapping_add(b);
            a = b;
            b = next;
        }
        fibs
    }
}
