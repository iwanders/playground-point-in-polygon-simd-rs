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
#[derive(Debug, Clone, Copy)]
pub struct EdgeVector {
    lower: [f64; 4],
    upper: [f64; 4],
    sub: [f64; 4],
    slope: [f64; 4],
}

impl Default for EdgeVector {
    fn default() -> Self {
        let mut n: EdgeVector = EdgeVector{
            lower: [0.0, 0.0, 0.0, 0.0],
            upper: [0.0, 0.0, 0.0, 0.0],
            sub: [0.0, 0.0, 0.0, 0.0],
            slope: [0.0, 0.0, 0.0, 0.0],
        };
        for vi in 0..4 {
            n.upper[vi] = f64::NEG_INFINITY;
            n.lower[vi] = f64::INFINITY;
        }
        n
    }
}

// Will at least be 4 * 4 * 4 = 4 * 16 = 64 bytes long.
use std::arch::x86_64::{__m256d, __m256i};
impl EdgeVector {
    pub fn combine(edges: &[&Edge]) -> Vec<EdgeVector> {
        let mut edge_v: Vec<EdgeVector> = vec![];
        for i in 0..edges.len() {
            if i % 4 == 0 {
                let n: EdgeVector = Default::default();
                edge_v.push(n);
            }
            edge_v[i / 4].lower[i % 4] = edges[i].lower;
            edge_v[i / 4].upper[i % 4] = edges[i].upper;
            edge_v[i / 4].sub[i % 4] = edges[i].sub;
            edge_v[i / 4].slope[i % 4] = edges[i].slope;
        }
        edge_v
    }

    /// Returns a vector of [-1, -1, -1, -1] for everything that crosses.
    pub fn calculate_crossings(&self, crossings_totals: &mut __m256i, tx: &__m256d, ty: &__m256d) {
        unsafe {
            use std::arch::x86_64::*;
            let iy = _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&self.lower[0]));
            let jy = _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&self.upper[0]));
            let sub = _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&self.sub[0]));
            let slope = _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&self.slope[0]));

            //  edge.iy <= test.1
            let above_lower = _mm256_cmp_pd(iy, *ty, _CMP_LE_OQ);
            // trace!("above_lower {}", pd(&above_lower));
            // test.1 <= edge.jy
            let below_upper = _mm256_cmp_pd(*ty, jy, _CMP_LT_OQ);
            // trace!("below_upper {}", pd(&below_upper));

            let in_range = _mm256_and_pd(above_lower, below_upper);

            // Or the other way around.
            // edge.iy <= test.1 <= edge.jy

            // let c = test.0 < (((test.1 * edge.slope + edge.sub) ));
            let right = _mm256_fmadd_pd(*ty, slope, sub);
            let t_l_right = _mm256_cmp_pd(*tx, right, _CMP_LT_OQ);
            // trace!("t_l_right {}", pd(&t_l_right));
            // println!("jskldfjsd");

            // Now, we mask all three together.
            let crosses = _mm256_and_pd(in_range, t_l_right);

            // Cast this to integers, which gets is 0 for not true, -1 for true;

            let crossess_i = _mm256_castpd_si256(crosses);
            *crossings_totals = _mm256_sub_epi64(*crossings_totals, crossess_i);
        }
    }

    pub fn calculate_crossings_range(
        crossings_totals: &mut __m256i,
        tx: &__m256d,
        ty: &__m256d,
        edges: &[EdgeVector],
    ) {
        // Search the left side of i mid up to left endpoint > v
        for v in edges.iter() {
            unsafe {
                let _f64_bitsset = f64::from_ne_bytes(0xFFFFFFFF_FFFFFFFFu64.to_ne_bytes());
                let _mask_1111 =
                    _mm256_set_pd(_f64_bitsset, _f64_bitsset, _f64_bitsset, _f64_bitsset);

                use std::arch::x86_64::*;
                let iy = _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&v.lower[0]));
                let jy = _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&v.upper[0]));
                let sub = _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&v.sub[0]));
                let slope = _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&v.slope[0]));

                //  edge.iy <= test.1
                let above_lower = _mm256_cmp_pd(iy, *ty, _CMP_LE_OQ);
                // trace!("above_lower {}", pd(&above_lower));
                // test.1 <= edge.jy
                let below_upper = _mm256_cmp_pd(*ty, jy, _CMP_LT_OQ);
                // trace!("below_upper {}", pd(&below_upper));

                let in_range = _mm256_and_pd(above_lower, below_upper);
                // let in_range = above_lower;

                if true {
                    if _mm256_testz_pd(in_range, _mask_1111) != 0 {
                        break;
                    }
                }
                // Or the other way around.
                // edge.iy <= test.1 <= edge.jy

                // let c = test.0 < (((test.1 * edge.slope + edge.sub) ));
                let right = _mm256_fmadd_pd(*ty, slope, sub);
                let t_l_right = _mm256_cmp_pd(*tx, right, _CMP_LT_OQ);
                // trace!("t_l_right {}", pd(&t_l_right));
                // println!("jskldfjsd");

                // Now, we mask all three together.
                let crosses = _mm256_and_pd(in_range, t_l_right);

                // Cast this to integers, which gets is 0 for not true, -1 for true;

                let crossess_i = _mm256_castpd_si256(crosses);
                *crossings_totals = _mm256_sub_epi64(*crossings_totals, crossess_i);
            }
        }
    }
}

#[repr(C)]
#[derive(Default, Debug, Clone, Copy)]
pub struct Edge {
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

    #[allow(dead_code)]
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

    fn sort_imid_left<'a>(edges: &[&'a Edge]) -> Vec<&'a Edge> {
        let mut left = edges.iter().map(|z| (z.lower, z)).collect::<Vec<_>>();
        left.sort_by(|&a, &b| a.0.partial_cmp(&b.0).unwrap());
        left.iter().map(|(_, a)| *a).copied().collect()
    }

    fn sort_imid_right<'a>(edges: &[&'a Edge]) -> Vec<&'a Edge> {
        let mut right = edges.iter().map(|z| (z.upper, z)).collect::<Vec<_>>();
        right.sort_by(|&a, &b| b.0.partial_cmp(&a.0).unwrap());
        right.iter().map(|(_, a)| *a).copied().collect()
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

    _padding: usize,
}
// at least 8 * 6 = 40 long

#[derive(Debug)]
pub struct EdgeTree {
    nodes: Vec<Branch>,
    edges_vector: Vec<EdgeVector>,
}

impl EdgeTree {
    pub fn new(vertices: &[(f64, f64)]) -> Self {
        let mut nodes = vec![];
        let mut edges_vector = vec![];

        // If there's no work to do, don't do work.
        if vertices.is_empty() {
            nodes.push(Branch::default());
            return EdgeTree {
                nodes,
                edges_vector,
            };
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
        nodes.push(Branch::default()); // push the first placeholder node
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
            // let sorted_imid = Edge::sort_imid(&i_mid);

            let sorted_imid_left = Edge::sort_imid_left(&i_mid);
            let sorted_imid_right = Edge::sort_imid_right(&i_mid);

            // Shortcut if there's four or less intervals, make a non-branching node and
            // shove them all into a vector.
            // This is not actually faster? O_o... no this costs 10% why!?
            if true && v.intervals.len() <= 4 {
                let imid_index = edges_vector.len();
                let imid_count = 1;
                let e = EdgeVector::combine(&v.intervals).pop().unwrap();
                edges_vector.push(e);
                edges_vector.push(e);
                let branch = Branch {
                    pivot,
                    left: None,
                    right: None,
                    imid_count,
                    imid_index,
                    _padding: 0,
                };
                nodes[v.precursor] = branch;
                continue;
            }

            let (imid_count, imid_index) = if imid_count != 0 {
                let imid_start = edges_vector.len();

                // Length is always the same...
                let combined_left = EdgeVector::combine(&sorted_imid_left);
                edges_vector.extend(combined_left.iter().map(|v| *v));

                let combined_right = EdgeVector::combine(&sorted_imid_right);
                edges_vector.extend(combined_right.iter().map(|v| *v));

                (combined_left.len(), imid_start)
            } else {
                (imid_count, 0)
            };

            // Determine what to do with left.
            let left = if i_left.is_empty() {
                None
            } else {
                let left = nodes.len();
                nodes.push(Branch::default()); // left node
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
                nodes.push(Branch::default()); // left node
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
                _padding: 0,
            };
            nodes[v.precursor] = branch;
        }

        EdgeTree {
            nodes,
            edges_vector,
        }
    }

    pub fn inside(&self, p: &(f64, f64)) -> bool {
        use std::arch::x86_64::*;
        struct RecurseState {
            // px: f64,
            tx: __m256d,
            ty: __m256d,
            crossings_totals: __m256i,
            py: f64,
        }

        fn recurser<'a>(
            o: &mut RecurseState,
            index: usize,
            nodes: &'a [Branch],
            edges_vector: &'a [EdgeVector],
        ) {
            let Branch {
                pivot,
                left,
                right,
                imid_count,
                imid_index,
                _padding,
            } = &nodes[index];
            if o.py < *pivot {
                // Search the left side of i mid up to left endpoint > v
                if *imid_count != 0 {
                    EdgeVector::calculate_crossings_range(
                        &mut o.crossings_totals,
                        &o.tx,
                        &o.ty,
                        &edges_vector[*imid_index..*imid_index + *imid_count],
                    );
                }

                if let Some(left_index) = left {
                    recurser(o, left_index.get(), nodes, edges_vector)
                }
            } else {
                // Search the right side of i mid up to right endpoint < v
                if *imid_count != 0 {
                    EdgeVector::calculate_crossings_range(
                        &mut o.crossings_totals,
                        &o.tx,
                        &o.ty,
                        &edges_vector[*imid_index + imid_count..*imid_index + 2 * *imid_count],
                    );
                }

                if let Some(right_index) = right {
                    recurser(o, right_index.get(), nodes, edges_vector)
                }
            }
        }

        let crossings_totals = unsafe { _mm256_set_epi64x(0, 0, 0, 0) };
        let ty = unsafe { _mm256_set1_pd(p.1) };
        let tx = unsafe { _mm256_set1_pd(p.0) };

        let mut o = RecurseState {
            py: p.1,
            ty,
            tx,
            crossings_totals,
        };
        recurser(&mut o, 0, &self.nodes, &self.edges_vector);
        // println!("r: {r:?}");

        let mut cross_normal = [0i64; 4];
        unsafe {
            _mm256_storeu_si256(
                std::mem::transmute::<_, *mut __m256i>(&mut cross_normal[0]),
                o.crossings_totals,
            )
        };
        let t: i64 = cross_normal.iter().sum();
        t.abs() % 2 == 1
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
