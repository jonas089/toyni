//! TOYNI-02: BabyBear arithmetic/equality are wrong for non-canonical values (see REPORT.md).

use toyni::babybear::{BabyBear, BABYBEAR_PRIME as P};
use toyni::ext::Ext;

const DELTA: u64 = (1u64 << 32) % P as u64;

fn canon(x: BabyBear) -> u32 {
    (x * BabyBear::one()).value
}

#[test]
fn equality_is_not_field_equality() {
    let p_as_zero = BabyBear { value: P };
    assert_eq!(canon(p_as_zero), 0);
    assert_ne!(p_as_zero, BabyBear::zero());
}

#[test]
fn byte_decoding_is_malleable() {
    let a = BabyBear::from_bytes(&5u32.to_le_bytes());
    let b = BabyBear::from_bytes(&(5 + P).to_le_bytes());
    assert_eq!(a, b, "two different encodings decode to the same element");
}

#[test]
#[cfg_attr(debug_assertions, ignore = "panics with overflow in debug builds")]
fn addition_wraps_and_is_wrong() {
    let a = BabyBear { value: u32::MAX }; // the element 2^32 - 1 mod p
    let sum = a + BabyBear::one();
    let expected = (((u32::MAX as u64) + 1) % P as u64) as u32; // = DELTA
    assert_eq!(expected as u64, DELTA);
    assert_ne!(canon(sum), expected, "u32 wrap changed the field result");
}

#[test]
#[cfg_attr(debug_assertions, ignore = "panics with overflow in debug builds")]
fn deep_and_constraint_see_different_values() {
    // OOD coordinate 1 of t_z sent as t + p (same element t).
    let t = 12345u32;
    let t_z = Ext::new([BabyBear::zero(), BabyBear { value: t + P }, BabyBear::zero(), BabyBear::zero()]);
    // DEEP numerator for an opened base-field value 0: (0 - t_z), coordinate 1.
    let deep = canon((Ext::zero() - t_z).c[1]);
    // Constraint side: (0 - (0 + t_z)) uses an addition first, which does not wrap here.
    let constraint = canon((Ext::zero() - (Ext::zero() + t_z)).c[1]);
    let field_neg_t = P - t;
    assert_eq!(constraint, field_neg_t);
    assert_eq!(deep as u64, (field_neg_t as u64 + DELTA) % P as u64, "DEEP is off by DELTA");
}
