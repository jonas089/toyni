use sha2::{Digest, Sha256};

#[derive(Debug)]
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
        let mut current_level = self.leaves.clone();
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

                let mut combined = Vec::new();
                combined.extend_from_slice(&left);
                combined.extend_from_slice(&right);
                next_level.push(sha_digest(&combined));
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

pub fn verify_merkle_proof(leaf: &Vec<u8>, proof: &MerkleProof, root: &Vec<u8>) -> bool {
    let mut current_hash = leaf.clone();

    for (sibling, is_right) in proof.path.iter().zip(proof.position.iter()) {
        let mut combined = Vec::new();
        if *is_right {
            combined.extend_from_slice(sibling);
            combined.extend_from_slice(&current_hash);
        } else {
            combined.extend_from_slice(&current_hash);
            combined.extend_from_slice(sibling);
        }
        current_hash = sha_digest(&combined);
    }

    current_hash == *root
}

fn sha_digest(data: &[u8]) -> Vec<u8> {
    let mut hasher = Sha256::new();
    hasher.update(data);
    hasher.finalize().to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merkle_proof_verification() {
        // Create leaves
        let leaves = vec![
            sha_digest(&1u64.to_le_bytes()),
            sha_digest(&2u64.to_le_bytes()),
            sha_digest(&3u64.to_le_bytes()),
            sha_digest(&4u64.to_le_bytes()),
        ];

        let tree = MerkleTree::new(leaves);
        let root = tree.root().unwrap();

        // Test proof for each leaf
        for i in 0..4 {
            let proof = tree.get_proof(i).unwrap();
            let leaf = sha_digest(&(i as u64 + 1).to_le_bytes());
            assert!(verify_merkle_proof(&leaf, &proof, &root));
        }
    }

    #[test]
    fn test_merkle_proof_odd_leaves() {
        // Test with odd number of leaves
        let leaves = vec![
            sha_digest(&1u64.to_le_bytes()),
            sha_digest(&2u64.to_le_bytes()),
            sha_digest(&3u64.to_le_bytes()),
        ];

        let tree = MerkleTree::new(leaves);
        let root = tree.root().unwrap();

        // Test proof for each leaf
        for i in 0..3 {
            let proof = tree.get_proof(i).unwrap();
            let leaf = sha_digest(&(i as u64 + 1).to_le_bytes());
            assert!(verify_merkle_proof(&leaf, &proof, &root));
        }
    }

    #[test]
    fn test_merkle_proof_single_leaf() {
        // Test with single leaf
        let leaves = vec![sha_digest(&1u64.to_le_bytes())];
        let tree = MerkleTree::new(leaves);
        let root = tree.root().unwrap();

        let proof = tree.get_proof(0).unwrap();
        let leaf = sha_digest(&1u64.to_le_bytes());
        assert!(verify_merkle_proof(&leaf, &proof, &root));
    }
}
