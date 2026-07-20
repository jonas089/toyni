//! SHA-256 Merkle tree commitments with domain separation and optional
//! per-leaf salts (for zero-knowledge hiding of committed evaluations).
//!
//! Leaf hashes are tagged `0x00` and internal nodes `0x01`, so a node hash
//! can never be presented as a leaf or vice versa. Leaf counts are always
//! powers of two in this library, so the tree is perfect.

use sha2::{Digest, Sha256};

pub type Hash = [u8; 32];

const LEAF_TAG: u8 = 0x00;
const NODE_TAG: u8 = 0x01;

pub fn hash_leaf(data: &[u8]) -> Hash {
    let mut h = Sha256::new();
    h.update([LEAF_TAG]);
    h.update(data);
    h.finalize().into()
}

pub fn hash_node(left: &Hash, right: &Hash) -> Hash {
    let mut h = Sha256::new();
    h.update([NODE_TAG]);
    h.update(left);
    h.update(right);
    h.finalize().into()
}

/// A perfect binary Merkle tree over pre-serialized leaves.
#[derive(Debug, Clone)]
pub struct MerkleTree {
    /// `levels[0]` = leaf hashes, `levels.last()` = `[root]`.
    pub levels: Vec<Vec<Hash>>,
}

impl MerkleTree {
    pub fn new(leaves: &[Vec<u8>]) -> Self {
        assert!(leaves.len().is_power_of_two(), "leaf count must be a power of two");
        let leaf_hashes: Vec<Hash> = hash_leaves(leaves);
        Self::from_leaf_hashes(leaf_hashes)
    }

    /// Build from a packed buffer of `count` fixed-size leaves.
    pub fn from_packed(packed: &[u8], leaf_len: usize) -> Self {
        assert_eq!(packed.len() % leaf_len, 0);
        let count = packed.len() / leaf_len;
        assert!(count.is_power_of_two(), "leaf count must be a power of two");
        let leaf_hashes: Vec<Hash> = {
            #[cfg(feature = "parallel")]
            {
                use rayon::prelude::*;
                if count >= 1 << 12 {
                    packed.par_chunks(leaf_len).map(hash_leaf).collect()
                } else {
                    packed.chunks(leaf_len).map(hash_leaf).collect()
                }
            }
            #[cfg(not(feature = "parallel"))]
            packed.chunks(leaf_len).map(hash_leaf).collect()
        };
        Self::from_leaf_hashes(leaf_hashes)
    }

    pub fn from_leaf_hashes(leaf_hashes: Vec<Hash>) -> Self {
        assert!(leaf_hashes.len().is_power_of_two());
        let mut levels = vec![leaf_hashes];
        while levels.last().unwrap().len() > 1 {
            let prev = levels.last().unwrap();
            let next: Vec<Hash> = build_level(prev);
            levels.push(next);
        }
        Self { levels }
    }

    pub fn root(&self) -> Hash {
        self.levels.last().unwrap()[0]
    }

    pub fn leaf_count(&self) -> usize {
        self.levels[0].len()
    }

    /// Authentication path for a leaf: sibling hashes bottom-up. The side of
    /// each sibling is implied by the bits of `index`.
    pub fn prove(&self, index: usize) -> Vec<Hash> {
        assert!(index < self.leaf_count());
        let mut path = Vec::with_capacity(self.levels.len() - 1);
        let mut i = index;
        for level in &self.levels[..self.levels.len() - 1] {
            path.push(level[i ^ 1]);
            i >>= 1;
        }
        path
    }
}

fn hash_leaves(leaves: &[Vec<u8>]) -> Vec<Hash> {
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        if leaves.len() >= 1 << 12 {
            return leaves.par_iter().map(|l| hash_leaf(l)).collect();
        }
    }
    leaves.iter().map(|l| hash_leaf(l)).collect()
}

fn build_level(prev: &[Hash]) -> Vec<Hash> {
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        if prev.len() >= 1 << 12 {
            return prev
                .par_chunks(2)
                .map(|pair| hash_node(&pair[0], &pair[1]))
                .collect();
        }
    }
    prev.chunks(2).map(|pair| hash_node(&pair[0], &pair[1])).collect()
}

/// Verify an authentication path produced by [`MerkleTree::prove`].
pub fn verify_path(root: &Hash, index: usize, leaf_data: &[u8], path: &[Hash]) -> bool {
    if index >= (1usize << path.len()) {
        return false;
    }
    let mut hash = hash_leaf(leaf_data);
    let mut i = index;
    for sibling in path {
        hash = if i & 1 == 0 {
            hash_node(&hash, sibling)
        } else {
            hash_node(sibling, &hash)
        };
        i >>= 1;
    }
    hash == *root
}

#[cfg(test)]
mod tests {
    use super::*;

    fn leaves(n: usize) -> Vec<Vec<u8>> {
        (0..n as u64).map(|i| i.to_le_bytes().to_vec()).collect()
    }

    #[test]
    fn prove_verify_all_positions() {
        for n in [1usize, 2, 8, 64] {
            let ls = leaves(n);
            let tree = MerkleTree::new(&ls);
            let root = tree.root();
            for (i, leaf) in ls.iter().enumerate() {
                let path = tree.prove(i);
                assert!(verify_path(&root, i, leaf, &path));
                // Wrong leaf fails.
                assert!(!verify_path(&root, i, b"junk", &path));
                // Wrong index fails.
                if n > 1 {
                    assert!(!verify_path(&root, i ^ 1, leaf, &path));
                }
            }
        }
    }

    #[test]
    fn leaf_node_domain_separation() {
        let ls = leaves(2);
        let tree = MerkleTree::new(&ls);
        let root = tree.root();
        // Presenting the root as a single leaf yields a different commitment.
        let masquerade = MerkleTree::new(&[root.to_vec()]);
        assert_ne!(masquerade.root(), root);
    }

    #[test]
    #[should_panic]
    fn non_power_of_two_rejected() {
        MerkleTree::new(&leaves(3));
    }
}
