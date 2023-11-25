//! A combination of an interval tree and the precomputed edges for SIMD operations.
/*
    Combine;
        - Precomputation for edge checks
        - SIMD operations
        - Interval tree to determine what checks.

    Since we can check 4 edges for the roughly the cost of one, we can terminate the tree whenever
    there's just 4 entries left.

    If at the Imid ranges we have less than four edges, we can lift ANY of the edges from lower in
    the tree up to ensure Imid is populated with multiples of four.

    Idealy
        - Strong data locality; everything in a single vector?
        - Cheap construction... O(n log n) is what literature claims.

*/

// Instead of 4 indepent vectors, we still interleave the individual chunks, that way the
// values we are loading into simd vectors are likely in the cache after the first iy read.
#[repr(C)]
#[derive(Default, Debug, Clone, Copy)]
struct EdgeVector {
    iy: [f64; 4],
    jy: [f64; 4],
    sub: [f64; 4],
    slope: [f64; 4],
}
// Will at least be 4 * 4 * 4 = 4 * 16 = 64 bytes long.

#[repr(C)]
#[derive(Default, Debug, Clone, Copy)]
struct Edge {
    lower: f64,
    upper: f64,
    sub: f64,
    slope: f64,
}
// At least 32 long.

impl Edge {
    pub fn new(vj: &(f64, f64), vi: &(f64, f64)) -> Self {
        let slope = (vj.0 - vi.0) / (vj.1 - vi.1);
        // Min and max here to account for edges that are rising & lowering,
        // avoids a second compare later.
        let lower = vi.1.min(vj.1);
        let upper = vj.1.max(vi.1);
        Edge {
            lower,
            upper,
            sub: -vi.1 * slope + vi.0,
            slope,
        }
    }
}

use std::num::NonZeroUsize;

#[derive(Debug, Clone, Default)]
struct Branch {
    /// Index to elements below of this pivot.
    left: Option<NonZeroUsize>,

    /// The actual center, the median of the values below this.
    pivot: f64,

    /// Index to elements above or equal to this pivot.
    right: Option<NonZeroUsize>,

    /// Number of imid elements
    imid_count: usize,

    /// Index where Imid starts in the main vector.
    imid_index: Option<NonZeroUsize>,
}
// at least 8 * 5 = 40 long

#[derive(Debug, Clone, Default)]
enum Node {
    #[default]
    Placeholder,
    Branch(Branch),
    Vector(EdgeVector),
}

struct EdgeTree {
    nodes: Vec<Node>,
}

impl EdgeTree {
    pub fn new(vertices: &[(f64, f64)]) -> Self {

        fn get_median(intervals: &[(f64, f64)]) -> f64 {
            let mut z = intervals
                .iter()
                .map(|se| [se.0, se.1])
                .flatten()
                .collect::<Vec<_>>();
            z.sort_by(|&a, &b| a.partial_cmp(&b).unwrap());
            z[z.len() / 2]
        }
        fn get_intervals(
            intervals: &[(f64, f64)],
            v: f64,
        ) -> Vec<(f64, f64)> {
            intervals
                .iter()
                .filter(|se| se.0 <= v && v <= se.1)
                .cloned()
                .collect()
        }
        fn get_left(
            intervals: &[(f64, f64)],
            v: f64,
        ) -> Vec<(f64, f64)> {
            intervals.iter().filter(|se| se.1 < v).cloned().collect()
        }
        fn get_right(
            intervals: &[(f64, f64)],
            v: f64,
        ) -> Vec<(f64, f64)> {
            intervals.iter().filter(|se| v < se.0).cloned().collect()
        }

        fn sort_imid(intervals: &[(f64, f64)]) -> Vec<f64> {
            let mut left = intervals.iter().map(|z| z.0).collect::<Vec<_>>();
            left.sort_by(|&a, &b| a.partial_cmp(&b).unwrap());
            let mut right = intervals.iter().map(|z| z.1).collect::<Vec<_>>();
            right.sort_by(|&a, &b| b.partial_cmp(&a).unwrap());
            left.iter().chain(right.iter()).cloned().collect()
        }


        let mut nodes = vec![];
        // If there's no work to do, don't do work.
        if vertices.is_empty() {
            nodes.push(Node::Branch(Branch::default()));
            return EdgeTree { nodes };
        }

        // First, convert the edges to vertices.
        let mut edges = vec![];
        let mut i = 0;
        let mut j = vertices.len() - 1;
        while i < vertices.len() {
            edges.push(Edge::new(&vertices[j], &vertices[i]));
            j = i;
            i += 1;
        }


        /*
        #[derive(Debug)]
        /// Helper struct for the nodes yet to be processed
        struct ProcessNode {
            intervals: Vec<Edge>,
            precursor: usize,
        }


        // We use a deque, such that we can insert in the rear and pop from the front.
        // This ensures that we don't get a depth first tree.
        nodes.push(Node::default()); // push the first placeholder node
                                     // Push the list of ondices to work on.
        use std::collections::VecDeque;
        let mut to_process: VecDeque<ProcessNode> = VecDeque::new();
        to_process.push_back(ProcessNode {
            intervals: intervals.to_vec(),
            precursor: 0,
        });

        while let Some(v) = to_process.pop_front() {
            // Determine the three interval options.
            let pivot = get_median(&v.intervals);
            let i_mid = get_intervals(&v.intervals, pivot);
            let i_left = get_left(&v.intervals, pivot);
            let i_right = get_right(&v.intervals, pivot);
            let count = i_mid.len();
            let sorted_imid = sort_imid(&i_mid);
            // println!("");
            // println!("pivot: {pivot}");
            // println!("i_mid: {i_mid:?}");
            // println!("i_left: {i_left:?}");
            // println!("i_right: {i_right:?}");

            // Determine what to do with left.
            let left = if i_left.is_empty() {
                None
            } else {
                let left = nodes.len();
                nodes.push(Node::default()); // left node
                to_process.push_back(ProcessNode {
                    intervals: i_left,
                    precursor: left,
                });
                NonZeroUsize::new(left)
            };
            let right = if i_right.is_empty() {
                None
            } else {
                let right = nodes.len();
                nodes.push(Node::default()); // left node
                to_process.push_back(ProcessNode {
                    intervals: i_right,
                    precursor: right,
                });
                NonZeroUsize::new(right)
            };

            // Finally, update the precursor, replacing the placeholder with a split pointing to
            // the correct nodes.
            nodes[v.precursor] = Node {
                pivot,
                left,
                right,
                count,
                sorted_imid,
            };
        }
        */





        EdgeTree{nodes}
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_edge_tree_creation() {
    }
}
