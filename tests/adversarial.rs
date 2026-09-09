//! A malicious prover, not a tamper of an honest proof.
//!
//! This attack verified successfully before the DEEP composition was given
//! per-term coefficients. Every rejection test in the unit suite mutates a
//! *valid* proof, which is exactly why it was not caught: the forged proof is
//! internally consistent by construction.

use rand::Rng;
use toyni::babybear::BabyBear;
use toyni::ext::Ext;
use toyni::fibonacci::{
    deep_value, MerkleOpening, MerkleOpeningExt, QueryProof, StarkProof, BLOWUP, COSET_SHIFT,
    MASK_DEGREE, NUM_DEEP_TERMS, NUM_QUERIES,
};
use toyni::math::domain::BabyBearDomain;
use toyni::math::fri::fri_fold_ext;
use toyni::math::mask::mask_poly_base;
use toyni::merkle::MerkleTree;
use toyni::transcript::FiatShamirTranscript;
use toyni::verifier::StarkVerifier;

/// **DEEP composition forgery.**
///
/// With every DEEP term weighted 1, the composition is
/// `(Q + T + T.g + T.g^2 - sigma)/(x - z)`, so the commitments constrain only
/// the *sum* `sigma` of the four claimed out-of-domain values. Together with
/// the constraint identity at `z` that is two equations in four unknowns, and
/// the remaining freedom is enough to make a uniformly random trace verify.
///
/// This runs the attack in full — garbage trace, garbage quotient, an honest
/// DEEP layer and honest FRI over the forged claims. With independent
/// coefficients drawn *after* the claims are absorbed, the DEEP numerator no
/// longer vanishes at `z`, so the layer is not a low-degree codeword and the
/// final-layer constancy check fires.
#[test]
fn deep_forgery_is_rejected() {
    let trace_len = 64usize;
    let lde_size = trace_len * BLOWUP;
    let domain = BabyBearDomain::new(trace_len);
    let shifted = BabyBearDomain::new(lde_size).get_coset(BabyBear::new(COSET_SHIFT));
    let g = domain.group_gen();
    let xs = shifted.elements();
    let mut rng = rand::thread_rng();

    // A trace that is not a Fibonacci sequence at all.
    let column: Vec<BabyBear> = (0..trace_len).map(|_| BabyBear::random(&mut rng)).collect();
    let mut trace_poly = domain.ifft(&column);
    mask_poly_base(&mut trace_poly, trace_len, MASK_DEGREE, &mut rng);
    let trace_lde = shifted.fft(&trace_poly);
    let (trace_tree, trace_salts) = salted(&trace_lde, &mut rng);

    // A "quotient" that is not C/Z_H: a random polynomial of the same degree.
    let q_coeffs: Vec<BabyBear> = (0..MASK_DEGREE + 2)
        .map(|_| BabyBear::random(&mut rng))
        .collect();
    let q_evals = shifted.fft(&q_coeffs);
    let (q_tree, q_salts) = salted(&q_evals, &mut rng);

    let trace_commitment = trace_tree.root().unwrap();
    let quotient_commitment = q_tree.root().unwrap();
    let mut ts = FiatShamirTranscript::new();
    ts.absorb_commitment(&trace_commitment);
    ts.absorb_commitment(&quotient_commitment);
    let z = derive_z(&mut ts);
    let (gz, ggz) = (z.mul_base(g), z.mul_base(g).mul_base(g));

    // True evaluations of what was actually committed.
    let t_z_true = z.eval_base_at_ext(&trace_poly);
    let t_gz_true = gz.eval_base_at_ext(&trace_poly);
    let t_ggz_true = ggz.eval_base_at_ext(&trace_poly);
    let sigma = t_z_true + t_gz_true + t_ggz_true + z.eval_base_at_ext(&q_coeffs);

    // Forge the claims. Keep t_z and t_gz true and solve the two equations
    //   (i)  (t_ggz - t_gz - t_z) * B(z) = q_z * Z_H(z)   [the AIR at z]
    //   (ii) t_z + t_gz + t_ggz + q_z = sigma             [the old DEEP's only tie]
    // for (t_ggz, q_z). Both are linear in t_ggz, so this is a 2x2 solve.
    let b_z = (z - Ext::from(g.pow((trace_len - 1) as u64)))
        * (z - Ext::from(g.pow((trace_len - 2) as u64)));
    let k = b_z * z.eval_base_at_ext(&domain.vanishing_poly_coeffs()).inverse();
    let c = t_z_true + t_gz_true;
    let t_ggz = (sigma - c + c * k) * (Ext::one() + k).inverse();
    let q_z = (t_ggz - c) * k;
    let (t_z, t_gz) = (t_z_true, t_gz_true);

    assert_eq!(t_z + t_gz + t_ggz + q_z, sigma, "the sum condition holds");
    assert_ne!(t_ggz, t_ggz_true, "the forged claim differs from the truth");

    for v in [t_z, t_gz, t_ggz, q_z] {
        ts.absorb_ext(v);
    }
    // The coefficients are only available *now*, after the claims are
    // committed. That is what defeats the attack.
    let mut coeffs = [Ext::zero(); NUM_DEEP_TERMS];
    for c in coeffs.iter_mut() {
        *c = ts.squeeze_ext_challenge();
    }

    // Built through the same `deep_value` the honest prover and the verifier
    // use, so this really is the strongest DEEP layer the attacker can commit
    // to for these claims.
    let d_evals: Vec<Ext> = (0..lde_size)
        .map(|i| {
            deep_value(
                Ext::from(xs[i]),
                z,
                Ext::from(trace_lde[i]),
                Ext::from(trace_lde[(i + BLOWUP) % lde_size]),
                Ext::from(trace_lde[(i + 2 * BLOWUP) % lde_size]),
                Ext::from(q_evals[i]),
                t_z,
                t_gz,
                t_ggz,
                q_z,
                &coeffs,
            )
        })
        .collect();

    // FRI, exactly as the honest prover does it.
    let mut layers = vec![d_evals];
    let mut trees = Vec::new();
    let mut fri_commitments: Vec<Vec<u8>> = Vec::new();
    let mut fold_xs = xs.clone();
    let final_size = lde_size / (trace_len + MASK_DEGREE).next_power_of_two();
    loop {
        let tree = ext_tree(layers.last().unwrap());
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

    let query_proofs = ts
        .squeeze_indices(NUM_QUERIES, lde_size / 2)
        .into_iter()
        .map(|qi| {
            let mut idx = qi;
            let fri_openings = (1..layers.len() - 1)
                .map(|l| {
                    let half = layers[l].len() / 2;
                    idx %= half;
                    (
                        open_ext(&trees[l], &layers[l], idx),
                        open_ext(&trees[l], &layers[l], idx + half),
                    )
                })
                .collect();
            QueryProof {
                index: qi,
                deep_opening: open_ext(&trees[0], &layers[0], qi),
                deep_opening_pair: open_ext(&trees[0], &layers[0], qi + lde_size / 2),
                trace_opening: open(&trace_tree, &trace_salts, &trace_lde, qi),
                trace_opening_g: open(
                    &trace_tree,
                    &trace_salts,
                    &trace_lde,
                    (qi + BLOWUP) % lde_size,
                ),
                trace_opening_gg: open(
                    &trace_tree,
                    &trace_salts,
                    &trace_lde,
                    (qi + 2 * BLOWUP) % lde_size,
                ),
                quotient_opening: open(&q_tree, &q_salts, &q_evals, qi),
                fri_openings,
            }
        })
        .collect();

    let forged = StarkProof {
        trace_len,
        lde_size,
        trace_commitment,
        quotient_commitment,
        t_z,
        t_gz,
        t_ggz,
        q_z,
        fri_commitments,
        fri_final_layer: layers.last().unwrap().clone(),
        query_proofs,
    };

    assert!(
        !StarkVerifier.verify(&forged),
        "a uniformly random trace was accepted"
    );
}

// ── helpers mirroring the prover's private tree builders ───────────────

fn salted(evals: &[BabyBear], rng: &mut impl Rng) -> (MerkleTree, Vec<Vec<u8>>) {
    let salts: Vec<Vec<u8>> = (0..evals.len())
        .map(|_| rng.r#gen::<[u8; 16]>().to_vec())
        .collect();
    let leaves = evals
        .iter()
        .zip(&salts)
        .map(|(v, s)| [s.as_slice(), &v.to_bytes()[..]].concat())
        .collect();
    (MerkleTree::new(leaves), salts)
}

fn ext_tree(evals: &[Ext]) -> MerkleTree {
    MerkleTree::new(evals.iter().map(|v| v.to_bytes().to_vec()).collect())
}

fn open(
    tree: &MerkleTree,
    salts: &[Vec<u8>],
    evals: &[BabyBear],
    i: usize,
) -> MerkleOpening {
    MerkleOpening {
        index: i,
        value: evals[i],
        proof: tree.get_proof(i).unwrap(),
        salt: salts[i].clone(),
    }
}

fn open_ext(tree: &MerkleTree, evals: &[Ext], i: usize) -> MerkleOpeningExt {
    MerkleOpeningExt {
        index: i,
        value: evals[i],
        proof: tree.get_proof(i).unwrap(),
    }
}

fn derive_z(ts: &mut FiatShamirTranscript) -> Ext {
    loop {
        let z = ts.squeeze_ext_challenge();
        if !z.is_base() {
            return z;
        }
    }
}
