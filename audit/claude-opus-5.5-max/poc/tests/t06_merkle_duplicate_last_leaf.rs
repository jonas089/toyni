//! TOYNI-08: odd-level duplication lets [a,b,c] and [a,b,c,c] share a root (see REPORT.md).

use toyni::merkle::{verify_merkle_proof, MerkleTree};

#[test]
fn three_and_four_leaf_trees_share_a_root() {
    let (a, b, c) = (vec![1u8], vec![2u8], vec![3u8]);
    let t3 = MerkleTree::new(vec![a.clone(), b.clone(), c.clone()]);
    let t4 = MerkleTree::new(vec![a, b, c.clone(), c.clone()]);
    assert_eq!(t3.root(), t4.root());

    let root3 = t3.root().unwrap();
    let proof = t4.get_proof(3).unwrap();
    assert!(verify_merkle_proof(&c, 3, 4, &proof, &root3), "index 3 opened in a 3-leaf commitment");
}
