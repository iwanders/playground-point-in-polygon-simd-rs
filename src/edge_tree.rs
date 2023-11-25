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
    lower: [f64; 4],
    upper: [f64; 4],
    sub: [f64; 4],
    slope: [f64; 4],
}
// Will at least be 4 * 4 * 4 = 4 * 16 = 64 bytes long.

impl EdgeVector {
    fn combine(edges: &[&Edge]) -> Vec<EdgeVector> {
        let mut edge_v: Vec<EdgeVector> = vec![];
        for i in 0..edges.len() {
            if i % 4 == 0 {
                edge_v.push(Default::default());
            }
            edge_v[i / 4].lower[i % 4] = edges[i].lower;
            edge_v[i / 4].upper[i % 4] = edges[i].upper;
            edge_v[i / 4].sub[i % 4] = edges[i].sub;
            edge_v[i / 4].slope[i % 4] = edges[i].slope;
        }
        edge_v
    }
}

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

    // Amount of references here is kinda ugly, but it basically means we operate on indice instead
    // of copying values around.

    fn get_median(edges: &[&Edge]) -> f64 {
        let mut z = edges
            .iter()
            .map(|se| [se.lower, se.upper])
            .flatten()
            .collect::<Vec<_>>();
        z.sort_by(|&a, &b| a.partial_cmp(&b).unwrap());
        z[z.len() / 2]
    }

    fn get_overlapping<'a>(edges: &[&'a Edge], v: f64) -> Vec<&'a Edge> {
        edges
            .iter()
            .filter(|e| e.lower <= v && v <= e.upper)
            .copied()
            .collect()
    }
    fn get_left<'a>(edges: &[&'a Edge], v: f64) -> Vec<&'a Edge> {
        edges.iter().filter(|e| e.upper < v).copied().collect()
    }
    fn get_right<'a>(edges: &[&'a Edge], v: f64) -> Vec<&'a Edge> {
        edges.iter().filter(|e| v < e.lower).copied().collect()
    }

    fn sort_imid<'a>(edges: &[&'a Edge]) -> Vec<&'a Edge> {
        let mut left = edges.iter().map(|z| (z.lower, z)).collect::<Vec<_>>();
        left.sort_by(|&a, &b| a.0.partial_cmp(&b.0).unwrap());
        let mut right = edges.iter().map(|z| (z.upper, z)).collect::<Vec<_>>();
        right.sort_by(|&a, &b| b.0.partial_cmp(&a.0).unwrap());
        left.iter()
            .copied()
            .chain(right.iter().copied())
            .map(|(_, a)| *a)
            .collect()
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
    imid_index: usize,
}
// at least 8 * 5 = 40 long

#[derive(Debug, Clone, Default)]
enum Node {
    #[default]
    Placeholder,
    Branch(Branch),
    Vector(EdgeVector),
    Edges(Vec<Edge>),
}

#[derive(Debug)]
pub struct EdgeTree {
    nodes: Vec<Node>,
}

impl EdgeTree {
    pub fn new(vertices: &[(f64, f64)]) -> Self {
        let mut nodes = vec![];

        // If there's no work to do, don't do work.
        if vertices.is_empty() {
            nodes.push(Node::Branch(Branch::default()));
            return EdgeTree { nodes };
        }

        // First, convert the vertices to edges.
        let mut edges = vec![];
        let mut i = 0;
        let mut j = vertices.len() - 1;
        while i < vertices.len() {
            edges.push(Edge::new(&vertices[j], &vertices[i]));
            j = i;
            i += 1;
        }

        // Now we have a set of edges which we can convert into a tree.
        #[derive(Debug)]
        /// Helper struct for the nodes yet to be processed
        struct ProcessNode<'a> {
            intervals: Vec<&'a Edge>,
            precursor: usize,
        }

        // We use a deque, such that we can insert in the rear and pop from the front.
        // This ensures that we don't get a depth first tree.
        nodes.push(Node::default()); // push the first placeholder node
                                     // Push the list of ondices to work on.
        use std::collections::VecDeque;
        let mut to_process: VecDeque<ProcessNode> = VecDeque::new();
        to_process.push_back(ProcessNode {
            intervals: edges.iter().collect::<Vec<_>>(),
            precursor: 0,
        });

        while let Some(v) = to_process.pop_front() {
            // Determine the three interval options.
            let pivot = Edge::get_median(&v.intervals);
            let i_mid = Edge::get_overlapping(&v.intervals, pivot);
            let i_left = Edge::get_left(&v.intervals, pivot);
            let i_right = Edge::get_right(&v.intervals, pivot);
            let imid_count = i_mid.len();
            let sorted_imid = Edge::sort_imid(&i_mid);

            /*
            let (imid_count, imid_index) = if imid_count != 0 {
                let combined = EdgeVector::combine(&sorted_imid);
                let imid_start = nodes.len();
                nodes.extend(combined.iter().map(|v| Node::Vector(*v)));
                (combined.len(), imid_start)
            } else {
                (imid_count, 0)
            };*/
            let imid_index = if imid_count != 0 {
                let imid_vec = nodes.len();
                nodes.push(Node::Edges(sorted_imid.iter().cloned().cloned().collect()));
                imid_vec
            } else {
                0
            };

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
            let branch = Branch {
                pivot,
                left,
                right,
                imid_count,
                imid_index,
            };
            nodes[v.precursor] = Node::Branch(branch);
        }

        EdgeTree { nodes }
    }

    pub fn inside(&self, p: &(f64, f64)) -> bool {
        use std::arch::x86_64::*;
        struct RecurseState {
            px: f64,
            py: f64,
            tx: __m256d,
            ty: __m256d,
            crossings_totals: __m256i,
            crossings_integer: usize,
        };

        fn recurser<'a>(o: &mut RecurseState, index: usize, nodes: &'a [Node]) {
            match &nodes[index] {
                Node::Branch(Branch {
                    pivot,
                    left,
                    right,
                    imid_count,
                    imid_index,
                }) => {
                    if o.py < *pivot {
                        // Search the left side of i mid up to left endpoint > v
                        if *imid_count != 0 {
                            if let Node::Edges(e) = &nodes[*imid_index] {
                                o.crossings_integer +=
                                    (e[0..*imid_count].iter().take_while(|e| e.lower <= o.py))
                                        .map(|edge| {
                                            if (edge.lower <= o.py)
                                                && (o.py < edge.upper)
                                                && o.px < (o.py * edge.slope + edge.sub)
                                            {
                                                1
                                            } else {
                                                0
                                            }
                                        })
                                        .sum::<usize>()
                            } else {
                                panic!("boo");
                            }
                        };

                        if let Some(left_index) = left {
                            recurser(o, left_index.get(), nodes)
                        }
                    } else {
                        // Search the right side of i mid up to right endpoint < v
                        if *imid_count != 0 {
                            if let Node::Edges(e) = &nodes[*imid_index] {
                                o.crossings_integer +=
                                    (e[*imid_count..].iter().take_while(|e| e.upper >= o.py))
                                        .map(|edge| {
                                            if (edge.lower <= o.py)
                                                && (o.py < edge.upper)
                                                && o.px < (o.py * edge.slope + edge.sub)
                                            {
                                                1
                                            } else {
                                                0
                                            }
                                        })
                                        .sum::<usize>()
                            } else {
                                panic!("boo");
                            }
                        };
                        if let Some(right_index) = right {
                            recurser(o, right_index.get(), nodes)
                        }
                    }
                }
                _ => unimplemented!(),
            }
        }

        let crossings_totals = unsafe { _mm256_set_epi64x(0, 0, 0, 0) };
        let ty = unsafe { _mm256_set1_pd(p.1) };
        let tx = unsafe { _mm256_set1_pd(p.0) };

        let mut o = RecurseState {
            px: p.0,
            py: p.1,
            ty,
            tx,
            crossings_totals,
            crossings_integer: 0,
        };
        recurser(&mut o, 0, &self.nodes);
        // println!("r: {r:?}");

        o.crossings_integer % 2 == 1
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn inside_edge_tree(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
        let z = EdgeTree::new(vertices);
        // println!("z: {z:?}");
        z.inside(test)
    }
    #[test]
    fn test_edge_tree_creation() {
        let triangle = vec![(0.0, 0.0), (0.0, 4.0), (6.0, 0.0), (0.0, 0.0)];
        for f in [inside_edge_tree] {
            assert_eq!(f(&triangle, &(1.0, 1.0)), true);
            assert_eq!(f(&triangle, &(3.0, 3.0)), false);

            assert_eq!(f(&triangle, &(0.0, 0.0)), true);
            assert_eq!(f(&triangle, &(0.0, 4.0)), false);
            assert_eq!(f(&triangle, &(5.0, 0.0)), true);

            assert_eq!(f(&triangle, &(0.0, 2.0)), true);
            assert_eq!(f(&triangle, &(3.0, 2.0)), false);
            assert_eq!(f(&triangle, &(3.0, 0.0)), true);
        }
    }

    #[test]
    fn test_edge_tree_square() {
        for f in [inside_edge_tree] {
            let square = vec![(0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0), (0.0, 0.0)];
            assert_eq!(f(&square, &(1.0, 1.0)), true);
            assert_eq!(f(&square, &(3.0, 3.0)), false);

            assert_eq!(f(&square, &(0.0, 0.0)), true);
            assert_eq!(f(&square, &(0.0, 2.0)), false);
            assert_eq!(f(&square, &(2.0, 2.0)), false);
            assert_eq!(f(&square, &(2.0, 0.0)), false);

            assert_eq!(f(&square, &(0.0, 1.0)), true);
            assert_eq!(f(&square, &(1.0, 2.0)), false);
            assert_eq!(f(&square, &(2.0, 1.0)), false);
            assert_eq!(f(&square, &(1.0, 0.0)), true);
        }
    }

    // Something to create a polygon from polar coordinates, that way it always a valid polygon.
    fn create_circle_parts(segments: &[(f64, f64)]) -> Vec<(f64, f64)> {
        let mut v = Vec::<(f64, f64)>::new();
        for segment in segments.iter() {
            let theta_c = segment.0 * std::f64::consts::PI * 2.0;
            let r_c = segment.1;
            v.push((theta_c.cos() * r_c, theta_c.sin() * r_c));
        }
        // Closing point.
        v.push(v[0]);
        v
    }

    #[test]
    fn test_edge_tree_circular() {
        use rand::prelude::*;
        use rand_xoshiro::rand_core::SeedableRng;
        use rand_xoshiro::Xoshiro256PlusPlus;
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);
        for _ in 0..100 {
            let mut distances: Vec<(f64, f64)> = vec![];
            let r = rng.gen::<f64>() * 100.0;
            let mut s = 0.0;
            distances.push((s, r * rng.gen::<f64>()));
            let segments: usize = rng.gen_range(3..100);
            for _ in 0..segments {
                let sa = (1.0 / (segments as f64)) * rng.gen::<f64>();
                s += sa;
                if s > 1.0 {
                    break;
                }
                distances.push((s, r * rng.gen::<f64>()));
            }
            let poly = create_circle_parts(&distances);
            // println!("Poly: {poly:?}");
            for _ in 0..1000 {
                let x = r * 1.25 * rng.gen::<f64>();
                let y = r * 1.25 * rng.gen::<f64>();
                let point = (x, y);
                let expected = crate::inside(&poly, &point);
                // println!("Point: {point:?}");
                assert_eq!(inside_edge_tree(&poly, &point), expected);
            }
        }
    }
}
