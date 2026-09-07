pub fn ntt(values: Vec<u32>, roots: Vec<u32>, modulus: u32) -> Vec<u32>{
    assert_eq!(values.len(), roots.len());
    // values: polynomial coefficients
    // roots: n roots of unity 
    let mut coeffs = reverse(values);
    // len starts at 2 for the first butterfly, then *2 each step; blocks grow by 2x each round and start at 1
    let mut len = 2;
    while len<=coeffs.len(){
        let half = len/2;
        // the stride is the steps which iterate over the roots 1x for each of the current stage's blocks
        // first we go over len = 2; stride = n/len = 4/2 with j=0 because half=1 and 0..1; => 0 idx of roots for each block
        // then we go over len = 4; stride = n/len = 4/4 with j=0,1 because half=2 and 0..2; => 0,1 idx of roots for each block
        let stride = coeffs.len()/len;
        // chunks_mut hands mutable slice, not a copy
        for block in coeffs.chunks_mut(len){
            let (lo, hi) = block.split_at_mut(half);
            for j in 0..half{
                // the current root
                let w = roots[j*stride] as u64;
                let a = lo[j]; // E[j]
                let b = hi[j] as u64 * w % modulus as u64; // w * O[j]
                // for each stage we evaluate the E and O side and update the slots in place, then re-use them for E and O in the 
                // next rounds
                lo[j] = (a as u32 + b as u32) % modulus; 
                // +modulus vanishes in reduction and handles case where b > a
                hi[j] = (a as u32 + modulus - b as u32) % modulus;
            }
        }
        len *=2;
    }
    coeffs
}

fn reverse(values: Vec<u32>) -> Vec<u32>{
    if values.len() <= 1{
        return values
    };

    let evens: Vec<u32> = values.iter().step_by(2).copied().collect();
    let odds: Vec<u32> = values.iter().skip(1).step_by(2).copied().collect();

    // build new vector from recursive results
    let mut out = reverse(evens);
    out.append(&mut reverse(odds));
    out
}

pub fn intt(values: Vec<u32>, roots: Vec<u32>, modulus: u32) -> Vec<u32> {
    let n = values.len();
    let mut inv_roots = Vec::with_capacity(n);
    inv_roots.push(roots[0]);
    inv_roots.extend(roots[1..].iter().rev());

    let coeffs = ntt(values, inv_roots, modulus);

    let n_inv = pow_mod(n as u32, modulus - 2, modulus) as u64;
    coeffs.into_iter()
        .map(|c| (c as u64 * n_inv % modulus as u64) as u32)
        .collect()
}

fn pow_mod(mut b: u32, mut e: u32, p: u32) -> u32 {
    let mut r = 1u32;
    b %= p;
    while e > 0 {
        if e & 1 == 1 { r = ((r as u64 * b as u64) % p as u64) as u32; }
        b = ((b as u64 * b as u64) % p as u64) as u32;
        e >>= 1;
    }
    r
}

/// Distinct prime factors of m, by trial division.
fn prime_factors(mut m: u32) -> Vec<u32> {
    let mut out = Vec::new();
    let mut q = 2u32;
    while q * q <= m {
        if m % q == 0 {
            out.push(q);
            while m % q == 0 { m /= q; }
        }
        q += 1;
    }
    if m > 1 { out.push(m); }
    out
}

/// Smallest multiplicative generator of the field F_p.
fn find_generator(p: u32) -> u32 {
    let factors = prime_factors(p - 1);
    (2..p)
        .find(|&g| factors.iter().all(|&q| pow_mod(g, (p - 1) / q, p) != 1))
        .expect("p must be prime")
}

/// The n roots of unity in F_p: omega^0, omega^1, ..., omega^(n-1).
pub fn roots_of_unity(n: u32, p: u32) -> Vec<u32> {
    assert!(n >= 1 && (p - 1) % n == 0, "n must divide p - 1");
    let g = find_generator(p);
    // power (p-1)/n since we want a sub-group of size n
    let omega = pow_mod(g, (p - 1) / n, p);

    let mut roots = Vec::with_capacity(n as usize);
    let mut cur = 1u32;
    for _ in 0..n {
        roots.push(cur);
        cur = ((cur as u64 * omega as u64) % p as u64) as u32;
    }
    roots
}

#[test]
fn roundtrip(){
    let prime = 257;
    let poly = vec![1,3,4,5,6,7,9,10];
    let n = poly.len() as u32;
    let roots = roots_of_unity(n, prime);
    let ntt_out = ntt(poly.clone(), roots.clone(), prime);
    let intt_out = intt(ntt_out, roots, prime);
    assert_eq!(poly, intt_out);
}
// ---------------------------------------------------------------------------
// BabyBear glue. Everything above is the plain u32 reference NTT copied
// verbatim; the wrappers below adapt it to the BabyBear type the rest of this
// crate uses so callers do not have to hand-roll the conversions.
// ---------------------------------------------------------------------------

use crate::babybear::BabyBear;

const P: u32 = crate::babybear::BABYBEAR_PRIME;

/// Forward NTT in place over BabyBear. Length must be a power of two.
pub fn ntt_babybear(values: &mut [BabyBear]) {
    let n = values.len();
    assert!(n.is_power_of_two(), "NTT size must be power of 2");
    let raw: Vec<u32> = values.iter().map(|v| v.value).collect();
    let out = ntt(raw, roots_of_unity(n as u32, P), P);
    for (v, r) in values.iter_mut().zip(out) {
        *v = BabyBear { value: r };
    }
}

/// Inverse NTT in place over BabyBear. Length must be a power of two.
pub fn intt_babybear(values: &mut [BabyBear]) {
    let n = values.len();
    assert!(n.is_power_of_two(), "NTT size must be power of 2");
    let raw: Vec<u32> = values.iter().map(|v| v.value).collect();
    let out = intt(raw, roots_of_unity(n as u32, P), P);
    for (v, r) in values.iter_mut().zip(out) {
        *v = BabyBear { value: r };
    }
}

#[cfg(test)]
mod babybear_tests {
    use super::*;

    #[test]
    fn babybear_roundtrip() {
        for log_n in [0u32, 1, 4, 8, 12] {
            let n = 1usize << log_n;
            let original: Vec<BabyBear> = (0..n).map(|i| BabyBear::new((i * 7 + 3) as u64)).collect();
            let mut values = original.clone();
            ntt_babybear(&mut values);
            intt_babybear(&mut values);
            assert_eq!(original, values, "roundtrip failed for n={n}");
        }
    }

    #[test]
    fn babybear_matches_naive_evaluation() {
        let n = 8usize;
        let roots = roots_of_unity(n as u32, P);
        let coeffs: Vec<BabyBear> = (0..n).map(|i| BabyBear::new((i * 13 + 5) as u64)).collect();
        let mut evals = coeffs.clone();
        ntt_babybear(&mut evals);
        for (k, root) in roots.iter().enumerate() {
            let x = BabyBear::new(*root as u64);
            let mut expected = BabyBear::zero();
            let mut xp = BabyBear::one();
            for c in &coeffs {
                expected = expected + *c * xp;
                xp = xp * x;
            }
            assert_eq!(evals[k], expected, "mismatch at evaluation point {k}");
        }
    }
}
