use sha2::{Digest, Sha256};

#[derive(Debug, Clone)]
pub struct MerkleProof {
    pub path: Vec<Vec<u8>>,
    pub position: Vec<bool>,
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
        let mut position = Vec::new();
        let mut current_index = index;

        // Start from the leaf level
        for level in &self.levels[..self.levels.len() - 1] {
            let sibling_index = if current_index % 2 == 0 {
                current_index + 1
            } else {
                current_index - 1
            };

            // If we're at the last node in an odd-sized level, use the node itself as sibling
            if sibling_index >= level.len() {
                path.push(level.get(current_index).unwrap().clone());
                position.push(true); // Treat as if sibling is on the right
            } else {
                path.push(level.get(sibling_index).unwrap().clone());
                position.push(current_index % 2 == 1);
            }

            current_index /= 2;
        }

        Some(MerkleProof { path, position })
    }

    pub fn root(&self) -> Option<Vec<u8>> {
        self.levels.last().unwrap().first().cloned().to_owned()
    }
}

pub fn verify_merkle_proof(leaf: Vec<u8>, proof: &MerkleProof, root: &Vec<u8>) -> bool {
    // Mirror build_tree: the supplied leaf is first domain-separated as a leaf,
    // then combined upward as internal nodes.
    let mut current_hash = hash_leaf(&leaf);

    for (sibling, is_right) in proof.path.iter().zip(proof.position.iter()) {
        current_hash = if *is_right {
            hash_node(sibling, &current_hash)
        } else {
            hash_node(&current_hash, sibling)
        };
    }

    current_hash == *root
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
        let leaves: Vec<Vec<u8>> = (1..=4).map(leaf).collect();
        let tree = MerkleTree::new(leaves.clone());
        let root = tree.root().unwrap();

        for i in 0..4 {
            let proof = tree.get_proof(i).unwrap();
            assert!(verify_merkle_proof(leaves[i].clone(), &proof, &root));
        }
    }

    #[test]
    fn test_merkle_proof_odd_leaves() {
        let leaves: Vec<Vec<u8>> = (1..=3).map(leaf).collect();
        let tree = MerkleTree::new(leaves.clone());
        let root = tree.root().unwrap();

        for i in 0..3 {
            let proof = tree.get_proof(i).unwrap();
            assert!(verify_merkle_proof(leaves[i].clone(), &proof, &root));
        }
    }

    #[test]
    fn test_merkle_proof_single_leaf() {
        let leaves = vec![leaf(1)];
        let tree = MerkleTree::new(leaves.clone());
        let root = tree.root().unwrap();

        let proof = tree.get_proof(0).unwrap();
        assert!(verify_merkle_proof(leaves[0].clone(), &proof, &root));
    }

    #[test]
    fn test_wrong_leaf_rejected() {
        let leaves: Vec<Vec<u8>> = (1..=4).map(leaf).collect();
        let tree = MerkleTree::new(leaves);
        let root = tree.root().unwrap();

        let proof = tree.get_proof(0).unwrap();
        // A different leaf value at position 0 must not verify.
        assert!(!verify_merkle_proof(leaf(99), &proof, &root));
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
