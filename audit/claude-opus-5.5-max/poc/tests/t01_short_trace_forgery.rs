//! TOYNI-01 regression: random traces are rejected at trace_len 8 and 16 (rate check).

mod common;

use common::{derive_z, finish, salted};
use toyni::babybear::BabyBear;
use toyni::ext::Ext;
use toyni::fibonacci::{deep_value, BLOWUP, COSET_SHIFT, NUM_DEEP_TERMS};
use toyni::math::domain::BabyBearDomain;
use toyni::transcript::FiatShamirTranscript;
use toyni::verifier::StarkVerifier;

fn forge(trace_len: usize) -> bool {
    let lde_size = trace_len * BLOWUP;
    let domain = BabyBearDomain::new(trace_len);
    let xs = BabyBearDomain::new(lde_size).get_coset(BabyBear::new(COSET_SHIFT)).elements();
    let g = domain.group_gen();
    let mut rng = rand::thread_rng();

    // Garbage trace and quotient values, committed as the prover would.
    let trace: Vec<BabyBear> = (0..lde_size).map(|_| BabyBear::random(&mut rng)).collect();
    let quot: Vec<BabyBear> = (0..lde_size).map(|_| BabyBear::random(&mut rng)).collect();
    let (t_tree, t_salts) = salted(&trace, &mut rng);
    let (q_tree, q_salts) = salted(&quot, &mut rng);

    let mut ts = FiatShamirTranscript::new();
    ts.absorb_commitment(&t_tree.root().unwrap());
    ts.absorb_commitment(&q_tree.root().unwrap());
    let z = derive_z(&mut ts);

    // Arbitrary trace claims; q(z) solved from the constraint identity at z.
    let (t_z, t_gz, t_ggz) = (Ext::random(&mut rng), Ext::random(&mut rng), Ext::random(&mut rng));
    let b = (z - Ext::from(g.pow((trace_len - 1) as u64))) * (z - Ext::from(g.pow((trace_len - 2) as u64)));
    let z_h = z.eval_base_at_ext(&domain.vanishing_poly_coeffs());
    let q_z = (t_ggz - (t_gz + t_z)) * b * z_h.inverse();
    for v in [t_z, t_gz, t_ggz, q_z] {
        ts.absorb_ext(v);
    }
    let mut coeffs = [Ext::zero(); NUM_DEEP_TERMS];
    for c in coeffs.iter_mut() {
        *c = ts.squeeze_ext_challenge();
    }

    // DEEP values the verifier will recompute on the first half...
    let half = lde_size / 2;
    let d = |i: usize| {
        deep_value(
            Ext::from(xs[i]),
            z,
            Ext::from(trace[i]),
            Ext::from(trace[(i + BLOWUP) % lde_size]),
            Ext::from(trace[(i + 2 * BLOWUP) % lde_size]),
            Ext::from(quot[i]),
            t_z,
            t_gz,
            t_ggz,
            q_z,
            &coeffs,
        )
    };
    let h1: Vec<Ext> = (0..half).map(d).collect();
    // ...and their unique degree < half interpolant on the second half.
    let w: Vec<BabyBear> = (0..half)
        .map(|j| (0..half).filter(|&k| k != j).fold(BabyBear::one(), |a, k| a * (xs[j] - xs[k])).inverse())
        .collect();
    let mut layer0 = h1.clone();
    for i in half..lde_size {
        let v = (0..half).fold(BabyBear::one(), |a, k| a * (xs[i] - xs[k]));
        let s = (0..half).fold(Ext::zero(), |a, j| a + h1[j].mul_base(w[j] * (xs[i] - xs[j]).inverse()));
        layer0.push(s.mul_base(v));
    }

    let proof = finish(
        ts,
        trace_len,
        &xs,
        (&t_tree, &t_salts, &trace),
        (&q_tree, &q_salts, &quot),
        [t_z, t_gz, t_ggz, q_z],
        layer0,
    );
    StarkVerifier.verify(&proof)
}

#[test]
fn random_trace_rejected_at_trace_len_16() {
    assert!(!forge(16), "forged proof accepted");
}

#[test]
fn random_trace_rejected_at_trace_len_8() {
    assert!(!forge(8), "forged proof accepted");
}

#[test]
fn same_attack_fails_at_trace_len_64() {
    assert!(!forge(64));
}
