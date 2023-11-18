// https://wrfranklin.org/Research/Short_Notes/pnpoly.html

/*

Argument	Meaning
nvert	Number of vertices in the polygon. Whether to repeat the first vertex at the end is discussed below.
vertx, verty	Arrays containing the x- and y-coordinates of the polygon's vertices.
testx, testy	X- and y-coordinate of the test point.

int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}
*/

/// Purest form of Franklin, exactly from the original logic.
pub fn inside(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
    let mut i = 0;
    let mut j = vertices.len() - 1;
    let mut inside = false;
    while i < vertices.len() {
        if ((vertices[i].1 > test.1) != (vertices[j].1 > test.1))
            && (test.0
                < (vertices[j].0 - vertices[i].0) * (test.1 - vertices[i].1)
                    / (vertices[j].1 - vertices[i].1)
                    + vertices[i].0)
        {
            inside = !inside;
        }
        j = i;
        i += 1;
    }
    inside
}

// maybe
// https://en.wikipedia.org/wiki/Interval_tree
// https://en.wikipedia.org/wiki/Segment_tree
#[repr(C)]
struct Edge {
    /// Lower y coordinate of current slope.
    iy: f64,
    /// Upper y coordinate of current slope.
    jy: f64,
    /// Subtraction for right hand compare.
    sub: f64,
    /// Slope for the right hand multiplication with y coordinate.
    slope: f64,
}

/// Precomputed version of Franklin
pub struct Precomputed {
    edges: Vec<Edge>,
}
impl Precomputed {
    pub fn new(vertices: &[(f64, f64)]) -> Self {
        let mut edges: Vec<Edge> = vec![];
        let mut i = 0;
        let mut j = vertices.len() - 1;
        while i < vertices.len() {
            let slope = (vertices[j].0 - vertices[i].0) / (vertices[j].1 - vertices[i].1);
            let lower_y = vertices[i].1;//.min(vertices[j].1);
            let upper_y = vertices[j].1;//.max(vertices[i].1);
            edges.push(Edge {
                // Min and max here to account for edges that are rising & lowering,
                // avoids a second compare later.
                iy: lower_y,
                jy: upper_y,
                sub: -vertices[i].1 + vertices[i].0,
                slope,
            });
            j = i;
            i += 1;
        }
        Self { edges }
    }

    /// Simple inside check, iterating through the edges.
    pub fn inside(self, test: &(f64, f64)) -> bool {
        let mut inside = false;
        for edge in self.edges.iter() {
            if ((edge.iy > test.1) != (edge.jy > test.1)) {
            // if (edge.iy <= test.1) && (test.1 < edge.jy) {
            // if (edge.iy <= test.1) && (test.1 < edge.jy) || ((edge.jy <= test.1) && (test.1 < edge.iy)) {
                let c = test.0 < (test.1 * edge.slope + edge.sub);
                if c {
                    inside = !inside
                }
            }
        }
        inside
    }

    /// SIMD inside check, iterating through edges in steps of four.
    pub fn inside_simd(self, test: &(f64, f64)) -> bool {
        use crate::print::pd;
        use crate::trace;
        use std::arch::x86_64::*;
        // Now, we think.
        /*
            if (edge.iy <= test.1) && (test.1 < edge.jy)
            let c = test.0 < (((test.1 + edge.sub) * edge.slope));


            // Also ok;
            if (edge.iy <= test.1) && (test.1 <= edge.jy)
            let c = test.0 < (((test.1 + edge.sub) * edge.slope));


            // then rewrite as fmadd; _mm256_fmadd_pd
            if (edge.iy <= test.1) && (test.1 <= edge.jy)
                let c = test.0 < (((test.1 * edge.slope + edge.sub) ));

            edge;
                iy: f64,
                jy: f64,
                sub: f64,
                slope: f64,

            // So now we have:
                edge.iy <= test.1
                test.1 <= edge.jy
                let c = test.0 < (((test.1 * edge.slope + edge.sub) ));
        */

        unsafe {
            let mut inside = false;
            let mut i = 0;

            // Step in fours;
            let ty = _mm256_set1_pd(test.1);
            let tx = _mm256_set1_pd(test.0);
            while i + 4 < self.edges.len() {
                let e = &self.edges[i..i + 4];
                // Could do a gather here, but lets get this working and then just make the struct
                // axis first instead of vertex first.
                let iy = _mm256_set_pd(e[3].iy, e[2].iy, e[1].iy, e[0].iy);
                let jy = _mm256_set_pd(e[3].jy, e[2].jy, e[1].jy, e[0].jy);
                let sub = _mm256_set_pd(e[3].sub, e[2].sub, e[1].sub, e[0].sub);
                let slope = _mm256_set_pd(e[3].slope, e[2].slope, e[1].slope, e[0].slope);

                //  edge.iy <= test.1
                let above_lower = _mm256_cmp_pd(iy, ty, _CMP_LE_OQ);
                // trace!("above_lower {}", pd(&above_lower));
                // test.1 <= edge.jy
                let below_upper = _mm256_cmp_pd(ty, jy, _CMP_LT_OQ);
                // trace!("below_upper {}", pd(&below_upper));

                // Or the other way around.
                // edge.iy <= test.1 <= edge.jy
                // 
                let inverted_coords_above_lower = _mm256_cmp_pd(jy, ty, _CMP_LT_OQ);
                let inverted_coords_below_upper = _mm256_cmp_pd(ty, iy, _CMP_LT_OQ);

                // let c = test.0 < (((test.1 * edge.slope + edge.sub) ));
                let right = _mm256_fmadd_pd(ty, slope, sub);
                let t_l_right = _mm256_cmp_pd(tx, right, _CMP_LT_OQ);
                // trace!("t_l_right {}", pd(&t_l_right));
                // println!("jskldfjsd");

                // Now, we mask all three together.
                let in_range = _mm256_and_pd(above_lower, below_upper);
                let inverted_in_range = _mm256_and_pd(inverted_coords_above_lower, inverted_coords_below_upper);
                let in_range = _mm256_or_pd(in_range, inverted_in_range);
                let crosses = _mm256_and_pd(in_range, t_l_right);

                // This section could be reduced to a vertical addition, followed by a
                // horizontal addition at the end.
                let bits = _mm256_movemask_pd(crosses);
                // trace!("bits: {:?}", bits);

                let count = _popcnt64(bits as i64);
                // trace!("count: {:?}", count);
                for _ in 0..count {
                    inside = !inside
                }

                i += 4;
            }

            // Finish the tail if not a multiple of four.
            while i < self.edges.len() {
                let edge = &self.edges[i];
                if (edge.iy <= test.1) && (test.1 < edge.jy) {
                    let c = test.0 < (test.1 * edge.slope + edge.sub);
                    if c {
                        inside = !inside
                    }
                }
                i += 1;
            }

            inside
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    // https://github.com/boostorg/geometry/blob/3f5c044abcc813a36d6af83465a9c086f9728a2f/test/strategies/franklin.cpp

    fn inside_precomputed(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
        let z = Precomputed::new(vertices);
        z.inside(test)
    }
    fn inside_precomputed_simd(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
        let z = Precomputed::new(vertices);
        z.inside_simd(test)
    }

    #[test]
    fn test_triangle() {
        for f in [inside, inside_precomputed, inside_precomputed_simd] {
            let triangle = vec![(0.0, 0.0), (0.0, 4.0), (6.0, 0.0), (0.0, 0.0)];
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
    fn test_square() {
        for f in [inside, inside_precomputed, inside_precomputed_simd] {
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
        for i in 0..segments.len() - 1 {
            let c = &segments[i];
            let n = &segments[i + 1];
            let theta_c = c.0 * std::f64::consts::PI * 2.0;
            let theta_n = n.0 * std::f64::consts::PI * 2.0;
            let r_c = c.1;
            let r_n = n.1;
            v.push((theta_c.cos() * r_c, theta_c.sin() * r_c));
            v.push((theta_n.cos() * r_n, theta_n.sin() * r_n));
        }
        // Closing point.
        v.push((segments[0].0.cos() * segments[0].1, segments[0].0.sin() * segments[0].1));
        // v.reverse();
        v
    }

    #[test]
    fn test_circle_parts() {
        let distances = [(0.0, 0.5),  (0.25, 0.5), (0.5, 1.0), (0.75, 1.0)];
        let poly = create_circle_parts(&distances);
        use rand::prelude::*;
        let mut rng = rand_xorshift::XorShiftRng::seed_from_u64(1);
        for _ in 0..10 {
            let x: f64 = rng.gen();
            let y: f64 = rng.gen();
            let point = (x, y);
            let expected = inside(&poly, &point);
            crate::trace!("Point: {:?}, exp: {}", point, expected);

            assert_eq!(inside_precomputed(&poly, &point), expected);
            assert_eq!(inside_precomputed_simd(&poly, &point), expected);            
        }
        
        
    }
}
