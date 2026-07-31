//! Zero-knowledge blinding of committed column polynomials.
//!
//! Adds `Z_H * R` in place, with `Z_H = x^n - 1` and `R` a fresh uniform
//! polynomial of `mask_degree` coefficients. `Z_H` vanishes on the size-`n`
//! trace domain, so the masked polynomial equals the original there
//! (constraints unchanged) but its off-domain openings are uniformly random.
//! Because `Z_H * R = x^n * R - R`, this subtracts `R` into the low
//! coefficients and adds it back shifted by `n`.

use crate::field::babybear::BabyBear;
use crate::field::babybear_ext::Ext;
use rand::Rng;

/// Blind a base-field column polynomial in place.
pub fn mask_poly_base(poly: &mut Vec<BabyBear>, n: usize, mask_degree: usize, rng: &mut impl Rng) {
    let r: Vec<BabyBear> = (0..mask_degree).map(|_| BabyBear::random(rng)).collect();
    if poly.len() < n + mask_degree {
        poly.resize(n + mask_degree, BabyBear::zero());
    }
    for i in 0..mask_degree {
        poly[i] = poly[i] - r[i];
        poly[n + i] = poly[n + i] + r[i];
    }
}

/// Blind an extension-field column polynomial in place.
pub fn mask_poly_ext(poly: &mut Vec<Ext>, n: usize, mask_degree: usize, rng: &mut impl Rng) {
    let r: Vec<Ext> = (0..mask_degree).map(|_| Ext::random(rng)).collect();
    if poly.len() < n + mask_degree {
        poly.resize(n + mask_degree, Ext::zero());
    }
    for i in 0..mask_degree {
        poly[i] = poly[i] - r[i];
        poly[n + i] = poly[n + i] + r[i];
    }
}
