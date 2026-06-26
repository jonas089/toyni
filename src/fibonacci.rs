use crate::babybear::BabyBear;
use crate::math::domain::BabyBearDomain;
use crate::math::fri::fri_fold;
use crate::math::polynomial::Polynomial;
use crate::merkle::{MerkleProof, MerkleTree};
use crate::program::trace::ExecutionTrace;
use crate::transcript::FiatShamirTranscript;
use rand::Rng;

/// Number of spot-check queries the verifier performs.
/// The masked DEEP composition is tested at rate 1/8 (D_BOUND/lde), so each
/// query gives ~3 bits and 44 queries → ~132 bits of security.
pub const NUM_QUERIES: usize = 44;
/// LDE blowup factor. Zero-knowledge masking raises the trace degree to
/// roughly 3·trace_len (see MASK_DEGREE), so the DEEP composition has degree
/// < 4·trace_len. Blowup 32 keeps the tested Reed–Solomon rate at 1/8
/// (D_BOUND/lde = 4n/32n) — i.e. masking is paid for in domain size, not
/// soundness.
pub const BLOWUP: usize = 32;
/// Coset shift used for the LDE domain.
pub const COSET_SHIFT: u64 = 7;
/// Number of random blinding coefficients added to the trace polynomial for
/// zero-knowledge: T̂(x) = T(x) + Z_H(x)·R(x) with R uniformly random of this
/// many coefficients. It must cover every trace evaluation the verifier sees —
/// 3 openings per query (T̂ at x, g·x, g²·x) plus the 3 out-of-domain points
/// (z, g·z, g²·z) — so that those evaluations are jointly uniform and reveal
/// nothing about the witness. 3·NUM_QUERIES + 3 = 135 suffices; a little margin.
pub const MASK_DEGREE: usize = 3 * NUM_QUERIES + 8;

// ── proof data structures ──────────────────────────────────────────────

/// Opening of a single position inside a Merkle-committed layer.
#[derive(Debug, Clone)]
pub struct MerkleOpening {
    pub index: usize,
    pub value: BabyBear,
    pub proof: MerkleProof,
    /// Per-leaf zero-knowledge salt: the leaf is committed as H(salt ‖ value).
    /// Revealed only for opened leaves, so an unopened sibling hash over a small
    /// base-field value cannot be brute-forced to recover the masked witness.
    pub salt: Vec<u8>,
}

/// All the data the verifier needs to check one query position across
/// every FRI layer.
#[derive(Debug, Clone)]
pub struct QueryProof {
    /// The original query index into the first FRI layer (= DEEP evals).
    pub index: usize,

    /// Openings into the DEEP layer (position and its pair).
    pub deep_opening: MerkleOpening,
    pub deep_opening_pair: MerkleOpening,

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
    pub fri_openings: Vec<(MerkleOpening, MerkleOpening)>,
}

/// A complete STARK proof that the verifier can check.
#[derive(Debug)]
pub struct StarkProof {
    // protocol parameters baked into the proof
    pub trace_len: usize,
    pub lde_size: usize,

    // commitments
    pub trace_commitment: Vec<u8>,
    pub quotient_commitment: Vec<u8>,

    // out-of-domain evaluations
    pub t_z: BabyBear,
    pub t_gz: BabyBear,
    pub t_ggz: BabyBear,
    pub q_z: BabyBear,

    // FRI commitments (layer 0 = DEEP evals, then each folded layer)
    pub fri_commitments: Vec<Vec<u8>>,
    /// The full final FRI layer, sent in the clear. The verifier checks it is a
    /// genuine low-degree (constant) codeword — the low-degree enforcement that
    /// per-query fold-consistency does not provide on its own.
    pub fri_final_layer: Vec<BabyBear>,

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

        let z_poly = Polynomial::new(domain.vanishing_poly_coeffs());
        let g = domain.group_gen();

        // ── 1. trace polynomial (+ zero-knowledge masking) ─────────────
        let domain_elements = domain.elements();
        let trace_poly = self.trace.interpolate_column(&domain_elements, 0);

        // Mask: T̂ = T + Z_H·R with R uniformly random. Z_H vanishes on the
        // trace domain H, so T̂ = T on H — the constraints (and thus soundness
        // and completeness) are untouched — but off H the LDE, query, and OOD
        // openings of T̂ are uniformly random and reveal nothing about the trace.
        let mut rng = rand::thread_rng();
        let r_poly = Polynomial::new(
            (0..MASK_DEGREE).map(|_| BabyBear::random(&mut rng)).collect(),
        );
        let trace_poly = trace_poly.add(&z_poly.multiply(&r_poly));

        // Evaluate trace poly on the shifted LDE domain and commit
        let shifted_elements = shifted_domain.elements();
        let trace_lde: Vec<BabyBear> = shifted_elements
            .iter()
            .map(|&x| trace_poly.evaluate(x))
            .collect();
        let trace_tree = build_merkle_tree(&trace_lde, &mut rng);
        let trace_commitment = trace_tree.root().unwrap();

        // ── 2. constraint & quotient ───────────────────────────────────
        let c_evals: Vec<BabyBear> = shifted_elements
            .iter()
            .map(|&x| {
                fibonacci_constraint(
                    trace_poly.evaluate(g * g * x),
                    trace_poly.evaluate(g * x),
                    trace_poly.evaluate(x),
                ) * boundary_constraint_1(x, g, trace_len)
                    * boundary_constraint_2(x, g, trace_len)
            })
            .collect();

        let c_poly = Polynomial::new(shifted_domain.ifft(&c_evals));

        let q_evals: Vec<BabyBear> = shifted_elements
            .iter()
            .map(|&x| c_poly.evaluate(x) / z_poly.evaluate(x))
            .collect();
        let q_poly = Polynomial::new(shifted_domain.ifft(&q_evals));

        let quotient_tree = build_merkle_tree(&q_evals, &mut rng);
        let quotient_commitment = quotient_tree.root().unwrap();

        // ── 3. Fiat-Shamir: derive OOD point z ────────────────────────
        let mut transcript = FiatShamirTranscript::new();
        transcript.absorb_commitment(&trace_commitment);
        transcript.absorb_commitment(&quotient_commitment);

        let z = derive_z_from_transcript(&mut transcript, &extended_domain, &shifted_domain);

        // ── 4. OOD evaluations ─────────────────────────────────────────
        let t_z = trace_poly.evaluate(z);
        let t_gz = trace_poly.evaluate(g * z);
        let t_ggz = trace_poly.evaluate(g * g * z);
        let q_z = q_poly.evaluate(z);

        // Sanity: constraint relation holds at z
        let c_z = fibonacci_constraint(t_ggz, t_gz, t_z)
            * boundary_constraint_1(z, g, trace_len)
            * boundary_constraint_2(z, g, trace_len);
        assert_eq!(
            c_z,
            q_z * z_poly.evaluate(z),
            "Constraint check at z failed"
        );

        // Feed OOD values into transcript
        transcript.absorb_field(t_z);
        transcript.absorb_field(t_gz);
        transcript.absorb_field(t_ggz);
        transcript.absorb_field(q_z);

        // ── 5. DEEP polynomial ─────────────────────────────────────────
        let d_evals: Vec<BabyBear> = shifted_elements
            .iter()
            .map(|&x| {
                let q_x = q_poly.evaluate(x);
                let t_x = trace_poly.evaluate(x);
                let t_gx = trace_poly.evaluate(g * x);
                let t_ggx = trace_poly.evaluate(g * g * x);
                (q_x - q_z) / (x - z)
                    + (t_ggx - t_ggz) / (x - z)
                    + (t_gx - t_gz) / (x - z)
                    + (t_x - t_z) / (x - z)
            })
            .collect();

        // ── 6. FRI folding with Merkle commits ────────────────────────
        let mut fri_layers: Vec<Vec<BabyBear>> = Vec::new();
        let mut fri_trees: Vec<SaltedTree> = Vec::new();
        let mut fri_commitments: Vec<Vec<u8>> = Vec::new();

        // Layer 0 = DEEP evaluations
        fri_layers.push(d_evals.clone());
        let tree0 = build_merkle_tree(&d_evals, &mut rng);
        let root0 = tree0.root().unwrap();
        transcript.absorb_commitment(&root0);
        fri_commitments.push(root0);
        fri_trees.push(tree0);

        let mut current = d_evals;
        let mut xs: Vec<BabyBear> = shifted_elements.clone();

        // Fold a FIXED number of rounds down to the degree-bound layer of size
        // lde_size / D_BOUND, where D_BOUND bounds deg(D). Masking raises the
        // trace degree, so deg(D) < D_BOUND = next_pow2(trace_len + MASK_DEGREE)
        // (≈ 4·trace_len). For an honest codeword that final layer is constant;
        // the verifier reads it whole and checks that constancy — THAT is what
        // enforces the low-degree bound. The round count is data-independent (no
        // early "constant yet?" break), so a malicious prover cannot stop early.
        let fri_degree_bound = (trace_len + MASK_DEGREE).next_power_of_two();
        let final_layer_size = lde_size / fri_degree_bound;
        while current.len() > final_layer_size {
            let beta = transcript.squeeze_challenge();

            let folded = fri_fold(&current, &xs, beta);

            // Square the x-coordinates for the next domain
            xs.truncate(folded.len());
            for x in &mut xs {
                *x = *x * *x;
            }

            fri_layers.push(folded.clone());
            let tree = if folded.len() == final_layer_size {
                build_unsalted_tree(&folded)
            } else {
                build_merkle_tree(&folded, &mut rng)
            };
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
            // Trace openings: T(x_i), T(g·x_i), T(g²·x_i)
            // g = ω_n = ω_N^BLOWUP, so g·x_qi corresponds to position
            // (qi + BLOWUP) % lde_size on the shifted domain.
            let idx_g = (qi + BLOWUP) % lde_size;
            let idx_gg = (qi + 2 * BLOWUP) % lde_size;
            let trace_opening = open_merkle(&trace_tree, &trace_lde, qi);
            let trace_opening_g = open_merkle(&trace_tree, &trace_lde, idx_g);
            let trace_opening_gg = open_merkle(&trace_tree, &trace_lde, idx_gg);

            // Quotient opening: Q(x_i)
            let quotient_opening = open_merkle(&quotient_tree, &q_evals, qi);

            // DEEP layer (layer 0) openings: position qi and its pair qi + half
            let half0 = fri_layers[0].len() / 2;
            let deep_opening = open_merkle(&fri_trees[0], &fri_layers[0], qi);
            let deep_opening_pair = open_merkle(&fri_trees[0], &fri_layers[0], qi + half0);

            // Intermediate FRI layers
            let mut fri_openings = Vec::new();
            let mut idx = qi;
            for layer_idx in 1..fri_layers.len() - 1 {
                let half = fri_layers[layer_idx].len() / 2;
                idx = idx % half;
                let op = open_merkle(&fri_trees[layer_idx], &fri_layers[layer_idx], idx);
                let op_pair =
                    open_merkle(&fri_trees[layer_idx], &fri_layers[layer_idx], idx + half);
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

/// Build a salted Merkle tree (leaf = salt ‖ value) — a hiding commitment.
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

/// Build an UNsalted tree (leaf = value), for the in-clear final FRI layer whose
/// root the verifier recomputes directly. Its values are public, so salting
/// would only block that recomputation.
fn build_unsalted_tree(evals: &[BabyBear]) -> SaltedTree {
    let leaves: Vec<Vec<u8>> = evals.iter().map(|v| v.to_bytes().to_vec()).collect();
    SaltedTree {
        tree: MerkleTree::new(leaves),
        salts: vec![Vec::new(); evals.len()],
    }
}

/// Open a salted Merkle tree at a given index.
fn open_merkle(tree: &SaltedTree, evals: &[BabyBear], index: usize) -> MerkleOpening {
    let proof = tree.tree.get_proof(index).expect("Index out of bounds");
    MerkleOpening {
        index,
        value: evals[index],
        proof,
        salt: tree.salts[index].clone(),
    }
}

/// Derive z via Fiat-Shamir, ensuring it avoids both the standard and shifted
/// LDE domains (and their g / g² shifts).
fn derive_z_from_transcript(
    transcript: &mut FiatShamirTranscript,
    extended_domain: &BabyBearDomain,
    shifted_domain: &BabyBearDomain,
) -> BabyBear {
    let ext_set: std::collections::HashSet<BabyBear> =
        extended_domain.elements().into_iter().collect();
    let shift_set: std::collections::HashSet<BabyBear> =
        shifted_domain.elements().into_iter().collect();
    let g = extended_domain.group_gen();

    loop {
        let z = transcript.squeeze_challenge();
        if !ext_set.contains(&z)
            && !shift_set.contains(&z)
            && !shift_set.contains(&(g * z))
            && !shift_set.contains(&(g * g * z))
        {
            return z;
        }
    }
}

// ── public constraint helpers (for the verifier) ───────────────────────

pub fn eval_fibonacci_constraint(t2: BabyBear, t1: BabyBear, t0: BabyBear) -> BabyBear {
    fibonacci_constraint(t2, t1, t0)
}

pub fn eval_boundary_1(x: BabyBear, g: BabyBear, n: usize) -> BabyBear {
    boundary_constraint_1(x, g, n)
}

pub fn eval_boundary_2(x: BabyBear, g: BabyBear, n: usize) -> BabyBear {
    boundary_constraint_2(x, g, n)
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
