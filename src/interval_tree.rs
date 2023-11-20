
#[derive(Debug, Clone, PartialEq, Ord, PartialOrd, Eq)]
pub struct IntervalId(pub usize);

#[derive(Debug, Clone)]
/// A node in our interval tree.
enum Node {
    Split {
        /// Index to elements below of this pivot.
        left: usize,
        /// The actual center value.
        pivot: f64,
        /// Index to elements above or equal to this pivot.
        right: usize,

        /// Storage of Imid
        data: Vec<((f64, f64), IntervalId)>
    },
    Leaf ,
    /// A placeholder, used during construction.
    Placeholder,
}

#[derive(Debug, Clone)]
pub struct IntervalTree {
    nodes: Vec<Node>,
}


// Construction is very similar to a KDTree, and some of the same concepts are applicable here
// see https://github.com/iwanders/playground-registration_icp_2d/blob/65dc9351fce30906fe7be1fadc57f95920ecd7f1/src/kdtree.rs


impl IntervalTree {
    pub fn new(intervals: &[((f64, f64), IntervalId)]) -> Self {
        // collect the endpoints into a vector;

        fn get_median(intervals: &[((f64, f64), IntervalId)]) -> f64 {
            let mut z = intervals.iter().map(|(se, _)| [se.0, se.1]).flatten().collect::<Vec<_>>();
            z.sort_by(|&a, &b| a.partial_cmp(&b).unwrap());
            z[z.len() / 2]
        }
        fn get_intervals(intervals: &[((f64, f64), IntervalId)], v: f64) -> Vec<((f64, f64), IntervalId)> {
            intervals.iter().filter(|a| a.0.0 <= v && v <= a.0.1).cloned().collect()
        }
        fn get_left(intervals: &[((f64, f64), IntervalId)], v: f64) -> Vec<((f64, f64), IntervalId)> {
            intervals.iter().filter(|a| a.0.1 < v).cloned().collect()
        }
        fn get_right(intervals: &[((f64, f64), IntervalId)], v: f64) -> Vec<((f64, f64), IntervalId)> {
            intervals.iter().filter(|a| v < a.0.0).cloned().collect()
        }
            
        let mut nodes = vec![];
    
        #[derive(Debug)]
        /// Helper struct for the nodes yet to be processed
        struct ProcessNode {
            intervals: Vec<((f64, f64), IntervalId)>,
            precursor: usize,
        }


        // We use a deque, such that we can insert in the rear and pop from the front.
        // This ensures that we don't get a depth first tree.
        nodes.push(Node::Placeholder); // push the first placeholder node
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
            println!("");
            println!("pivot: {pivot}");
            println!("i_mid: {i_mid:?}");
            println!("i_left: {i_left:?}");
            println!("i_right: {i_right:?}");

            // Determine what to do with left.
            let left = if i_left.is_empty() {
                // Drop the points into a container.
                nodes.push(Node::Leaf);
                nodes.len() - 1
            } else {
                let left = nodes.len();
                nodes.push(Node::Placeholder); // left node
                to_process.push_back(ProcessNode {
                    intervals: i_left,
                    precursor: left,
                });
                left
            };
            let right = if i_right.is_empty() {
                // Drop the points into a container.
                nodes.push(Node::Leaf);
                nodes.len() - 1
            } else {
                let right = nodes.len();
                nodes.push(Node::Placeholder); // left node
                to_process.push_back(ProcessNode {
                    intervals: i_right,
                    precursor: right,
                });
                right
            };

            // Finally, update the precursor, replacing the placeholder with a split pointing to
            // the correct nodes.
            nodes[v.precursor] = Node::Split {
                pivot,
                left,
                right,
                data: i_mid,
            };
        }

        IntervalTree{nodes}
    }

    pub fn intervals(&self, v: f64) -> Vec<IntervalId> {
        let mut r = vec![];

        fn recurser(o: &mut Vec<IntervalId>, v: f64, index: usize, nodes: &[Node]) {
            match &nodes[index] {
                Node::Placeholder => panic!("placeholder encountered during search"),
                Node::Split{
                    left, pivot, right, data
                } => {
                    if v <= *pivot {
                        // check imid, and traverse left.
                        o.extend(data.iter().filter_map(|((s, _e), z)| {if s <= &v {Some(z)} else {None}}).cloned());
                        recurser(o, v, *left, nodes)
                    } else {
                        o.extend(data.iter().filter_map(|((_s, e), z)| {if e >= &v {Some(z)} else {None}}).cloned());
                        recurser(o, v, *right, nodes)
                    }
                },
                Node::Leaf => {},
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
        intervals.iter().filter_map(|((s, e), i)| {
            if *s <= v && v <= *e {
                Some(i)
            } else {
                None
            }
        }).cloned().collect()
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

        for v in range(-1.5, 1.0, 0.01) {
            assert_intervals(&t.intervals(v), &get_interval_ids(&intervals, v));
        }
    }
}
