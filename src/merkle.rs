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

/// One opening for a whole set of leaves.
///
/// Queries into the same tree share almost every internal node, so sending a
/// full path per leaf sends the upper levels over and over. This sends each
/// needed node once, in the order both sides walk the tree, which is where most
/// of a FRI proof's size goes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MerkleMultiProof {
    /// Sibling nodes, bottom level first and ascending by index within a level.
    pub nodes: Vec<Vec<u8>>,
}

impl MerkleTree {
    /// Open several leaves at once. Indices need not be sorted or distinct.
    pub fn multi_proof(&self, indices: &[usize]) -> Option<MerkleMultiProof> {
        let mut known = sorted_distinct(indices);
        if known.last().is_some_and(|&i| i >= self.leaves.len()) {
            return None;
        }
        if known.is_empty() {
            return Some(MerkleMultiProof { nodes: Vec::new() });
        }

        let mut nodes = Vec::new();
        let mut level_size = self.leaves.len();
        for level in 0..depth_for(self.leaves.len()) {
            let mut parents = Vec::with_capacity(known.len());
            let mut i = 0;
            while i < known.len() {
                let idx = known[i];
                if idx == level_size - 1 && level_size % 2 == 1 {
                    // Last node of an odd level is its own sibling.
                } else if idx % 2 == 0 && known.get(i + 1) == Some(&(idx + 1)) {
                    // Both halves are known, so no node is needed.
                    i += 1;
                } else {
                    nodes.push(self.levels[level][idx ^ 1].clone());
                }
                parents.push(idx / 2);
                i += 1;
            }
            known = parents;
            level_size = (level_size + 1) / 2;
        }
        Some(MerkleMultiProof { nodes })
    }
}

/// Verify a batch opening.
///
/// `leaves` are the leaf preimages with the index each sits at. They may be in
/// any order and may repeat; a repeated index must carry the same leaf.
pub fn verify_multi_proof(
    leaves: &[(usize, Vec<u8>)],
    num_leaves: usize,
    proof: &MerkleMultiProof,
    root: &[u8],
) -> bool {
    if num_leaves == 0 || leaves.is_empty() {
        return false;
    }
    let mut sorted: Vec<(usize, Vec<u8>)> = leaves.to_vec();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));
    // A repeated index must agree with itself, or a prover could claim two
    // values for one position.
    sorted.dedup_by(|a, b| a.0 == b.0 && a.1 == b.1);
    if sorted.windows(2).any(|w| w[0].0 == w[1].0) {
        return false;
    }
    if sorted.last().is_some_and(|(i, _)| *i >= num_leaves) {
        return false;
    }

    let mut known: Vec<(usize, Vec<u8>)> =
        sorted.into_iter().map(|(i, leaf)| (i, hash_leaf(&leaf))).collect();
    let mut level_size = num_leaves;
    let mut feed = proof.nodes.iter();

    for _ in 0..depth_for(num_leaves) {
        let mut parents: Vec<(usize, Vec<u8>)> = Vec::with_capacity(known.len());
        let mut i = 0;
        while i < known.len() {
            let (idx, ref node) = known[i];
            let (left, right) = if idx == level_size - 1 && level_size % 2 == 1 {
                (node.clone(), node.clone())
            } else if idx % 2 == 0 && known.get(i + 1).map(|(j, _)| *j) == Some(idx + 1) {
                let pair = known[i + 1].1.clone();
                i += 1;
                (node.clone(), pair)
            } else {
                let Some(sibling) = feed.next() else {
                    return false;
                };
                if sibling.len() != 32 {
                    return false;
                }
                if idx % 2 == 0 {
                    (node.clone(), sibling.clone())
                } else {
                    (sibling.clone(), node.clone())
                }
            };
            parents.push((idx / 2, hash_node(&left, &right)));
            i += 1;
        }
        known = parents;
        level_size = (level_size + 1) / 2;
    }

    // Every supplied node must have been used, so a proof cannot carry slack.
    if feed.next().is_some() {
        return false;
    }
    known.len() == 1 && known[0].1 == *root
}

fn sorted_distinct(indices: &[usize]) -> Vec<usize> {
    let mut out = indices.to_vec();
    out.sort_unstable();
    out.dedup();
    out
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

    /// A batch opening must accept exactly what the individual openings do,
    /// and must be smaller.
    #[test]
    fn multi_proofs_agree_with_single_proofs() {
        for count in [1usize, 2, 3, 5, 8, 9, 16, 64] {
            let ls: Vec<Vec<u8>> = (1..=count as u64).map(leaf).collect();
            let tree = MerkleTree::new(ls.clone());
            let root = tree.root().unwrap();

            for step in [1usize, 2, 3, 7] {
                let indices: Vec<usize> = (0..count).step_by(step).collect();
                let batch = tree.multi_proof(&indices).unwrap();
                let opened: Vec<(usize, Vec<u8>)> =
                    indices.iter().map(|&i| (i, ls[i].clone())).collect();
                assert!(
                    verify_multi_proof(&opened, count, &batch, &root),
                    "count={count} step={step}"
                );

                // Never larger than sending a path each, and smaller as soon as
                // the paths overlap.
                let single: usize = indices
                    .iter()
                    .map(|&i| tree.get_proof(i).unwrap().path.len())
                    .sum();
                assert!(batch.nodes.len() <= single, "count={count} step={step}");
                if indices.len() > 1 && count > 2 {
                    assert!(batch.nodes.len() < single, "no saving at count={count}");
                }
            }
        }
    }

    #[test]
    fn batch_openings_are_bound_to_their_indices() {
        let ls: Vec<Vec<u8>> = (1..=16).map(leaf).collect();
        let tree = MerkleTree::new(ls.clone());
        let root = tree.root().unwrap();
        let indices = [1usize, 4, 5, 11];
        let batch = tree.multi_proof(&indices).unwrap();
        let good: Vec<(usize, Vec<u8>)> =
            indices.iter().map(|&i| (i, ls[i].clone())).collect();
        assert!(verify_multi_proof(&good, 16, &batch, &root));

        // A leaf moved to another of the opened positions must fail.
        let mut moved = good.clone();
        moved[0].1 = ls[4].clone();
        assert!(!verify_multi_proof(&moved, 16, &batch, &root));

        // Dropping or adding a leaf changes the shape and must fail.
        assert!(!verify_multi_proof(&good[..3], 16, &batch, &root));
        let mut extra = good.clone();
        extra.push((7, ls[7].clone()));
        assert!(!verify_multi_proof(&extra, 16, &batch, &root));

        // Two different values for one index must fail.
        let mut doubled = good.clone();
        doubled.push((1, ls[2].clone()));
        assert!(!verify_multi_proof(&doubled, 16, &batch, &root));

        // A proof carrying slack must fail.
        let mut padded = batch.clone();
        padded.nodes.push(vec![0u8; 32]);
        assert!(!verify_multi_proof(&good, 16, &padded, &root));
        let mut short = batch;
        short.nodes.pop();
        assert!(!verify_multi_proof(&good, 16, &short, &root));
    }

    /// Order must not matter to the caller, since query indices arrive in
    /// transcript order.
    #[test]
    fn batch_openings_ignore_input_order() {
        let ls: Vec<Vec<u8>> = (1..=32).map(leaf).collect();
        let tree = MerkleTree::new(ls.clone());
        let root = tree.root().unwrap();
        let batch = tree.multi_proof(&[9, 2, 30, 2]).unwrap();
        let shuffled = vec![
            (30usize, ls[30].clone()),
            (2, ls[2].clone()),
            (9, ls[9].clone()),
        ];
        assert!(verify_multi_proof(&shuffled, 32, &batch, &root));
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
