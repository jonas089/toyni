use sha2::{Digest, Sha256};

/// The sibling hashes from a leaf to the root, leaf level first.
///
/// The left/right decision at each level is deliberately *not* stored here: it
/// is derived from the leaf index the verifier asks about. An earlier version
/// kept prover-chosen direction bits and ignored the index, which let a prover
/// answer a query at position `i` with the leaf and path of any other position
/// `j` — making every query check vacuous. See
/// `tests::index_substitution_is_rejected`.
#[derive(Debug, Clone)]
pub struct MerkleProof {
    pub path: Vec<Vec<u8>>,
}

#[derive(Debug)]
pub struct MerkleTree {
    pub leaves: Vec<Vec<u8>>,
    pub levels: Vec<Vec<Vec<u8>>>,
}

impl MerkleTree {
    pub fn new(leaves: Vec<Vec<u8>>) -> Self {
        let mut tree = MerkleTree {
            leaves: leaves.clone(),
            levels: Vec::new(),
        };
        tree.build_tree();
        tree
    }

    pub fn build_tree(&mut self) {
        // Leaf level is the domain-separated hash of each supplied leaf. This
        // both hashes raw (e.g. unhashed field-element) leaves and tags them so
        // a leaf can never be reinterpreted as an internal node (see hash_leaf /
        // hash_node).
        let mut current_level: Vec<Vec<u8>> =
            self.leaves.iter().map(|leaf| hash_leaf(leaf)).collect();
        self.levels.push(current_level.clone());

        while current_level.len() > 1 {
            let mut next_level = Vec::new();
            for i in (0..current_level.len()).step_by(2) {
                let left = current_level.get(i).unwrap();
                let right = if i + 1 < current_level.len() {
                    current_level.get(i + 1).unwrap()
                } else {
                    current_level.get(i).unwrap() // Duplicate last node if odd number
                };
                next_level.push(hash_node(left, right));
            }
            current_level = next_level;
            self.levels.push(current_level.clone());
        }
    }

    pub fn get_proof(&self, index: usize) -> Option<MerkleProof> {
        if index >= self.leaves.len() {
            return None;
        }

        let mut path = Vec::new();
        let mut current_index = index;

        // Start from the leaf level
        for level in &self.levels[..self.levels.len() - 1] {
            // If we're at the last node in an odd-sized level, it is its own sibling
            let sibling = level
                .get(current_index ^ 1)
                .unwrap_or_else(|| level.get(current_index).unwrap());
            path.push(sibling.clone());
            current_index /= 2;
        }

        Some(MerkleProof { path })
    }

    pub fn root(&self) -> Option<Vec<u8>> {
        self.levels.last().unwrap().first().cloned().to_owned()
    }
}

/// Verify that `leaf` sits at `index` of a tree of `num_leaves` leaves.
///
/// The index is what makes this binding: the direction taken at each level is
/// computed from it, and the path must have exactly the depth `num_leaves`
/// implies. A path for a different position therefore cannot verify.
pub fn verify_merkle_proof(
    leaf: &[u8],
    index: usize,
    num_leaves: usize,
    proof: &MerkleProof,
    root: &[u8],
) -> bool {
    if num_leaves == 0 || index >= num_leaves || proof.path.len() != depth_for(num_leaves) {
        return false;
    }

    // Mirror build_tree: the supplied leaf is first domain-separated as a leaf,
    // then combined upward as internal nodes.
    let mut current_hash = hash_leaf(leaf);
    let mut idx = index;
    let mut level_size = num_leaves;

    for sibling in &proof.path {
        // The last node of an odd-sized level is its own sibling, and must be
        // presented as such; anything else is a malformed path.
        current_hash = if idx == level_size - 1 && level_size % 2 == 1 {
            if *sibling != current_hash {
                return false;
            }
            hash_node(&current_hash, &current_hash)
        } else if idx % 2 == 0 {
            hash_node(&current_hash, sibling)
        } else {
            hash_node(sibling, &current_hash)
        };
        idx /= 2;
        level_size = (level_size + 1) / 2;
    }

    current_hash == *root
}

/// Number of sibling hashes on a root path for a tree with `num_leaves` leaves.
pub fn depth_for(num_leaves: usize) -> usize {
    let mut depth = 0;
    let mut size = num_leaves;
    while size > 1 {
        size = (size + 1) / 2;
        depth += 1;
    }
    depth
}

/// Domain-separation tags keep the leaf and internal-node hash spaces disjoint,
/// so an internal node hash can never be presented as a leaf (or vice versa).
const LEAF_TAG: u8 = 0x00;
const NODE_TAG: u8 = 0x01;

/// Hash a leaf: `SHA256(0x00 || leaf)`.
fn hash_leaf(data: &[u8]) -> Vec<u8> {
    let mut hasher = Sha256::new();
    hasher.update([LEAF_TAG]);
    hasher.update(data);
    hasher.finalize().to_vec()
}

/// Hash an internal node: `SHA256(0x01 || left || right)`.
fn hash_node(left: &[u8], right: &[u8]) -> Vec<u8> {
    let mut hasher = Sha256::new();
    hasher.update([NODE_TAG]);
    hasher.update(left);
    hasher.update(right);
    hasher.finalize().to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn leaf(n: u64) -> Vec<u8> {
        n.to_le_bytes().to_vec()
    }

    #[test]
    fn test_merkle_proof_verification() {
        for count in [1usize, 2, 3, 4, 5, 8, 9] {
            let leaves: Vec<Vec<u8>> = (1..=count as u64).map(leaf).collect();
            let tree = MerkleTree::new(leaves.clone());
            let root = tree.root().unwrap();
            for i in 0..count {
                let proof = tree.get_proof(i).unwrap();
                assert!(
                    verify_merkle_proof(&leaves[i], i, count, &proof, &root),
                    "honest opening failed at {}/{}",
                    i,
                    count
                );
            }
        }
    }

    #[test]
    fn test_wrong_leaf_rejected() {
        let leaves: Vec<Vec<u8>> = (1..=4).map(leaf).collect();
        let tree = MerkleTree::new(leaves);
        let root = tree.root().unwrap();

        let proof = tree.get_proof(0).unwrap();
        // A different leaf value at position 0 must not verify.
        assert!(!verify_merkle_proof(&leaf(99), 0, 4, &proof, &root));
    }

    /// The exploit that motivated index binding: leaf `j` with leaf `j`'s
    /// authentication path must not verify at any position other than `j`.
    /// Without this, a prover can answer any query with any committed value and
    /// every FRI/DEEP query check becomes vacuous.
    #[test]
    fn index_substitution_is_rejected() {
        let leaves: Vec<Vec<u8>> = (1..=8).map(leaf).collect();
        let tree = MerkleTree::new(leaves.clone());
        let root = tree.root().unwrap();

        let proof = tree.get_proof(5).unwrap();
        assert!(verify_merkle_proof(&leaves[5], 5, 8, &proof, &root));
        for claimed in 0..8 {
            if claimed != 5 {
                assert!(
                    !verify_merkle_proof(&leaves[5], claimed, 8, &proof, &root),
                    "leaf 5 was accepted at position {}",
                    claimed
                );
            }
        }
    }

    /// The depth is pinned by the leaf count, so a truncated or padded path
    /// cannot verify.
    #[test]
    fn wrong_path_length_is_rejected() {
        let leaves: Vec<Vec<u8>> = (1..=8).map(leaf).collect();
        let tree = MerkleTree::new(leaves.clone());
        let root = tree.root().unwrap();

        let mut short = tree.get_proof(0).unwrap();
        short.path.pop();
        assert!(!verify_merkle_proof(&leaves[0], 0, 8, &short, &root));

        let mut long = tree.get_proof(0).unwrap();
        long.path.push(vec![0u8; 32]);
        assert!(!verify_merkle_proof(&leaves[0], 0, 8, &long, &root));
    }

    #[test]
    fn test_leaf_node_domain_separation() {
        // A two-leaf root is hash_node(hash_leaf(a), hash_leaf(b)). Because
        // leaves are tagged 0x00 and nodes 0x01, that node hash cannot be
        // reinterpreted as a leaf: committing to it as a single leaf yields a
        // different root, so an internal node can never masquerade as a leaf.
        let tree = MerkleTree::new(vec![leaf(1), leaf(2)]);
        let node_root = tree.root().unwrap();

        let masquerade = MerkleTree::new(vec![node_root.clone()]);
        assert_ne!(masquerade.root().unwrap(), node_root);
    }
}
