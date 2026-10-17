//! TOYNI-05: recover the whole witness from one proof by linear algebra (see REPORT.md).

use toyni::babybear::BabyBear;
use toyni::ext::Ext;
use toyni::fibonacci::{StarkProver, BLOWUP, COSET_SHIFT, MASK_DEGREE, NUM_DEEP_TERMS};
use toyni::math::domain::BabyBearDomain;
use toyni::program::trace::ExecutionTrace;
use toyni::transcript::FiatShamirTranscript;

/// Gaussian elimination over BabyBear; returns the unique solution if the system has full column rank.
fn solve(mut rows: Vec<Vec<BabyBear>>, unknowns: usize) -> Option<Vec<BabyBear>> {
    let mut r = 0;
    for col in 0..unknowns {
        let pivot = (r..rows.len()).find(|&i| !rows[i][col].is_zero())?;
        rows.swap(r, pivot);
        let inv = rows[r][col].inverse();
        for v in rows[r].iter_mut() {
            *v = *v * inv;
        }
        for i in 0..rows.len() {
            if i != r && !rows[i][col].is_zero() {
                let f = rows[i][col];
                for j in col..=unknowns {
                    let t = rows[r][j] * f;
                    rows[i][j] = rows[i][j] - t;
                }
            }
        }
        r += 1;
    }
    // Remaining rows must be consistent (0 = 0).
    if rows[r..].iter().any(|row| !row[unknowns].is_zero()) {
        return None;
    }
    Some((0..unknowns).map(|c| rows[c][unknowns]).collect())
}

fn powers(x: BabyBear, k: usize) -> Vec<BabyBear> {
    let mut out = Vec::with_capacity(k);
    let mut cur = BabyBear::one();
    for _ in 0..k {
        out.push(cur);
        cur = cur * x;
    }
    out
}

#[test]
fn witness_is_recovered_from_one_proof() {
    // Secret: a recurrence with random seeds.
    let n = 64;
    let mut rng = rand::thread_rng();
    let mut secret = vec![BabyBear::random(&mut rng), BabyBear::random(&mut rng)];
    while secret.len() < n {
        let k = secret.len();
        secret.push(secret[k - 1] + secret[k - 2]);
    }
    let mut trace = ExecutionTrace::new();
    trace.insert_column(secret.clone());
    let proof = StarkProver::new(trace).generate_proof(false);

    // Everything below uses only public parameters and the proof.
    let lde = proof.lde_size;
    let unknowns = n + MASK_DEGREE;
    let domain = BabyBearDomain::new(n);
    let xs = BabyBearDomain::new(lde).get_coset(BabyBear::new(COSET_SHIFT)).elements();
    let g = domain.group_gen();

    let mut ts = FiatShamirTranscript::new();
    ts.absorb_commitment(&proof.trace_commitment);
    ts.absorb_commitment(&proof.quotient_commitment);
    let z = loop {
        let z = ts.squeeze_ext_challenge();
        if !z.is_base() {
            break z;
        }
    };
    for v in [proof.t_z, proof.t_gz, proof.t_ggz, proof.q_z] {
        ts.absorb_ext(v);
    }
    let mut c = [Ext::zero(); NUM_DEEP_TERMS];
    for ci in c.iter_mut() {
        *ci = ts.squeeze_ext_challenge();
    }

    let mut rows: Vec<Vec<BabyBear>> = Vec::new();
    // Ext-valued linear equation sum_j coef[j] c_j = rhs -> 4 base equations.
    let mut push_ext = |coef: Vec<Ext>, rhs: Ext| {
        for k in 0..4 {
            let mut row: Vec<BabyBear> = coef.iter().map(|e| e.c[k]).collect();
            row.push(rhs.c[k]);
            rows.push(row);
        }
    };

    // OOD values: T^(z g^s) = t_s.
    for (s, t) in [(0u64, proof.t_z), (1, proof.t_gz), (2, proof.t_ggz)] {
        let pt = z.mul_base(g.pow(s));
        let mut coef = Vec::with_capacity(unknowns);
        let mut cur = Ext::one();
        for _ in 0..unknowns {
            coef.push(cur);
            cur = cur * pt;
        }
        push_ext(coef, t);
    }
    let mut base_rows: Vec<Vec<BabyBear>> = Vec::new();
    for qp in &proof.query_proofs {
        // Opened trace values.
        for op in [&qp.trace_opening, &qp.trace_opening_g, &qp.trace_opening_gg] {
            let mut row = powers(xs[op.index], unknowns);
            row.push(op.value);
            base_rows.push(row);
        }
        // DEEP pair value at x' = x_{qi ^ lde/2}:
        // d (x'-z) + c0 t_z + c1 t_gz + c2 t_ggz + c3 q_z
        //   = T^(x')(c0 - c3 r) + T^(gx')(c1 - c3 r) + T^(g^2 x')(c2 + c3 r),  r = B(x')/Z_H(x').
        let i = qp.index ^ (lde / 2);
        let (x0, x1, x2) = (xs[i], xs[(i + BLOWUP) % lde], xs[(i + 2 * BLOWUP) % lde]);
        let b = (x0 - g.pow((n - 1) as u64)) * (x0 - g.pow((n - 2) as u64));
        let r = b * (x0.pow(n as u64) - BabyBear::one()).inverse();
        let k0 = c[0] - c[3].mul_base(r);
        let k1 = c[1] - c[3].mul_base(r);
        let k2 = c[2] + c[3].mul_base(r);
        let (p0, p1, p2) = (powers(x0, unknowns), powers(x1, unknowns), powers(x2, unknowns));
        let coef: Vec<Ext> =
            (0..unknowns).map(|j| k0.mul_base(p0[j]) + k1.mul_base(p1[j]) + k2.mul_base(p2[j])).collect();
        let rhs = qp.deep_opening_pair.value * (Ext::from(x0) - z)
            + c[0] * proof.t_z
            + c[1] * proof.t_gz
            + c[2] * proof.t_ggz
            + c[3] * proof.q_z;
        push_ext(coef, rhs);
    }
    rows.extend(base_rows);
    assert!(rows.len() > unknowns, "{} equations for {} unknowns", rows.len(), unknowns);

    let coeffs = solve(rows, unknowns).expect("system not full rank");
    let recovered: Vec<BabyBear> = domain
        .elements()
        .iter()
        .map(|&w| coeffs.iter().rev().fold(BabyBear::zero(), |acc, &cj| acc * w + cj))
        .collect();
    assert_eq!(recovered, secret, "witness not recovered");
}
