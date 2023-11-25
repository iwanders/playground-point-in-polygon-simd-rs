#[derive(Debug, Clone, PartialEq, Ord, PartialOrd, Eq, Copy)]
pub struct IntervalId(pub usize);

// We can use non zero indices to ensure the option is zero cost, this is guaranteed because
// no node will ever refer to the root node, which is the starting point for the search.
use std::num::NonZeroUsize;

#[derive(Debug, Clone, Default)]
/// A node in our interval tree.
struct Node {
    /// Index to elements below of this pivot.
    left: Option<NonZeroUsize>,

    /// The actual center, the median of the values below this.
    pivot: f64,

    /// Index to elements above or equal to this pivot.
    right: Option<NonZeroUsize>,

    /// Count of entries in Imid.
    count: usize,

    /// The next vector holds both left boundaries as well as right boundaries.
    /// 0..count is increasing by left endpoint.
    /// count..end is decreasing by right endpoint.
    sorted_imid: Vec<(f64, IntervalId)>,
}

#[derive(Debug, Clone)]
pub struct IntervalTree {
    nodes: Vec<Node>,
}

// Construction is very similar to a KDTree, and some of the same concepts are applicable here
// see https://github.com/iwanders/playground-registration_icp_2d/blob/65dc9351fce30906fe7be1fadc57f95920ecd7f1/src/kdtree.rs

impl IntervalTree {
    pub fn new(intervals: &[((f64, f64), IntervalId)]) -> Self {
        // This creation can be cleaned up significantly, currently there's a lof copying around.

        fn get_median(intervals: &[((f64, f64), IntervalId)]) -> f64 {
            let mut z = intervals
                .iter()
                .map(|(se, _)| [se.0, se.1])
                .flatten()
                .collect::<Vec<_>>();
            z.sort_by(|&a, &b| a.partial_cmp(&b).unwrap());
            z[z.len() / 2]
        }
        fn get_intervals(
            intervals: &[((f64, f64), IntervalId)],
            v: f64,
        ) -> Vec<((f64, f64), IntervalId)> {
            intervals
                .iter()
                .filter(|a| a.0 .0 <= v && v <= a.0 .1)
                .cloned()
                .collect()
        }
        fn get_left(
            intervals: &[((f64, f64), IntervalId)],
            v: f64,
        ) -> Vec<((f64, f64), IntervalId)> {
            intervals.iter().filter(|a| a.0 .1 < v).cloned().collect()
        }
        fn get_right(
            intervals: &[((f64, f64), IntervalId)],
            v: f64,
        ) -> Vec<((f64, f64), IntervalId)> {
            intervals.iter().filter(|a| v < a.0 .0).cloned().collect()
        }

        fn sort_imid(intervals: &[((f64, f64), IntervalId)]) -> Vec<(f64, IntervalId)> {
            let mut left = intervals.iter().map(|z| (z.0 .0, z.1)).collect::<Vec<_>>();
            left.sort_by(|&a, &b| a.0.partial_cmp(&b.0).unwrap());
            let mut right = intervals.iter().map(|z| (z.0 .1, z.1)).collect::<Vec<_>>();
            right.sort_by(|&a, &b| b.0.partial_cmp(&a.0).unwrap());
            left.iter().chain(right.iter()).cloned().collect()
        }

        let mut nodes = vec![];

        #[derive(Debug)]
        /// Helper struct for the nodes yet to be processed
        struct ProcessNode {
            intervals: Vec<((f64, f64), IntervalId)>,
            precursor: usize,
        }

        if intervals.is_empty() {
            nodes.push(Node::default());
            return IntervalTree { nodes };
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

        IntervalTree { nodes }
    }

    pub fn intervals(&self, v: f64) -> Vec<IntervalId> {
        let mut r = vec![];

        fn recurser(o: &mut Vec<IntervalId>, v: f64, index: usize, nodes: &[Node]) {
            let Node {
                pivot,
                left,
                right,
                count,
                sorted_imid,
            } = &nodes[index];
            if v < *pivot {
                // Search the left side of i mid up to left endpoint > v
                o.extend(
                    sorted_imid[0..*count]
                        .iter()
                        .take_while(|(le, _)| le <= &v)
                        .map(|(_, i)| *i),
                );
                if let Some(left_index) = left {
                    recurser(o, v, left_index.get(), nodes)
                }
            } else {
                // Search the right side of i mid up to right endpoint < v
                o.extend(
                    sorted_imid[*count..]
                        .iter()
                        .take_while(|(re, _)| re >= &v)
                        .map(|(_, i)| *i),
                );
                if let Some(right_index) = right {
                    recurser(o, v, right_index.get(), nodes)
                }
            }
        }

        recurser(&mut r, v, 0, &self.nodes);

        r
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn get_interval_ids(intervals: &[((f64, f64), IntervalId)], v: f64) -> Vec<IntervalId> {
        intervals
            .iter()
            .filter_map(
                |((s, e), i)| {
                    if *s <= v && v <= *e {
                        Some(i)
                    } else {
                        None
                    }
                },
            )
            .cloned()
            .collect()
    }

    fn assert_intervals(a: &[IntervalId], b: &[IntervalId]) {
        let mut a = a.to_vec();
        let mut b = b.to_vec();
        a.sort();
        b.sort();
        assert_eq!(a, b);
    }
    fn range(start: f64, end: f64, interval: f64) -> Vec<f64> {
        let mut r = vec![];
        let mut v = start;
        r.push(start);
        while v < end {
            r.push(v);
            v += interval;
        }
        r.push(end);
        r
    }

    #[test]
    fn test_interval_tree_basic() {
        let intervals = [
            ((-1.0, -0.5), IntervalId(0)),
            ((-0.75, 0.0), IntervalId(1)),
            ((-0.5, 0.5), IntervalId(2)),
            ((0.25, 0.5), IntervalId(3)),
            ((0.4, 0.5), IntervalId(4)),
            ((0.0, 0.3), IntervalId(5)),
        ];
        let t = IntervalTree::new(&intervals);
        println!("t: {t:?}");

        for v in range(-1.5, 1.0, 0.01) {
            assert_intervals(&t.intervals(v), &get_interval_ids(&intervals, v));
        }
    }

    #[test]
    fn test_random_intervals() {
        use rand::prelude::*;
        use rand_xoshiro::rand_core::SeedableRng;
        use rand_xoshiro::Xoshiro256PlusPlus;
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(1);
        for _ in 0..100 {
            let interval_count: usize = rng.gen_range(0..100);

            let mut intervals = vec![];
            for i in 0..interval_count {
                let a = rng.gen::<f64>() * 100.0;
                let b = rng.gen::<f64>() * 100.0;
                let left = a.min(b);
                let right = a.max(b);

                intervals.push(((left, right), IntervalId(i)));
            }
            
            let t = IntervalTree::new(&intervals);
            // println!("t: {t:?}");

            let points = intervals.iter().map(|(a,_)| a.0).chain(intervals.iter().map(|(a,_)| a.1)).collect::<Vec<_>>();
            for v in points {
                assert_intervals(&t.intervals(v), &get_interval_ids(&intervals, v));
            }
            for v in range(0.0, 100.0, 0.01) {
                assert_intervals(&t.intervals(v), &get_interval_ids(&intervals, v));
            }
        }
    }
}
