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

pub fn inside_assume_closing(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
    let mut i = 1;
    let mut j = 0;
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
// Repr c to ensure we can read the first four doubles as simd vectors.
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

    // Only for debugging
    vi: (f64, f64),
    vj: (f64, f64),
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
            let lower_y = vertices[i].1; //.min(vertices[j].1);
            let upper_y = vertices[j].1; //.max(vertices[i].1);
            edges.push(Edge {
                // Min and max here to account for edges that are rising & lowering,
                // avoids a second compare later.
                iy: lower_y,
                jy: upper_y,
                sub: -vertices[i].1 * slope + vertices[i].0,
                slope,
                vi: vertices[i],
                vj: vertices[j],
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
            if (edge.iy > test.1) != (edge.jy > test.1) {
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
            while i + 4 <= self.edges.len() {
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
                let inverted_in_range =
                    _mm256_and_pd(inverted_coords_above_lower, inverted_coords_below_upper);
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

    #[test]
    fn test_single_edge() {
        let vi = (3.061616997868383e-17, 0.5);
        let vj = (0.5, 0.0);
        let test = (0.1449590873681229, 0.05982148134111476);

        let a = vi.1 > test.1;
        let b = vj.1 > test.1;
        let cr = (vj.0 - vi.0) * (test.1 - vi.1) / (vj.1 - vi.1) + vi.0;
        // crate::trace!("cr: {cr:?}");
        let c = test.0 < cr;
        assert_eq!(a, true);
        assert_eq!(b, false);
        assert_eq!(c, true);

        // now test the slope decomposed flavour.
        /*
            let cr = (vj.0 - vi.0) * (test.1 - vi.1) / (vj.1 - vi.1) + vi.0;
            let cr = (jx - ix) * (ty - iy) / (jy - iy) + ix;
        */

        let slope = (vj.0 - vi.0) / (vj.1 - vi.1);
        // crate::trace!("slope: {slope:?}");
        let lower_y = vi.1;
        let upper_y = vj.1;
        let sub = -vi.1 * slope + vi.0;
        // crate::trace!("sub: {sub:?}");
        let iy = lower_y;
        let jy = upper_y;

        let a = iy > test.1;
        let b = jy > test.1;
        let c = test.0 < (test.1 * slope + sub);
        assert_eq!(a, true);
        assert_eq!(b, false);
        assert_eq!(c, true);
    }
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
        for f in [
            inside,
            inside_precomputed,
            inside_precomputed_simd,
            inside_assume_closing,
        ] {
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
        for f in [
            inside,
            inside_precomputed,
            inside_precomputed_simd,
            inside_assume_closing,
        ] {
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

    #[test]
    fn test_polygon_found_with_closing_removal() {
        let poly = vec![
            (46.737189787317114, 0.0),
            (-5.487807116126276, 16.463760906592245),
            (-1.9544451225548618, 1.027559503579469),
            (46.737189787317114, 0.0),
        ];
        let point = (54.57356785726406, 85.87822568825256);
        assert_eq!(inside(&poly, &point), false);
        assert_eq!(inside_precomputed(&poly, &point), false);
        assert_eq!(inside_precomputed_simd(&poly, &point), false);
        assert_eq!(inside_assume_closing(&poly, &point), false);
        let point = (12.230435569131824, 3.1868743396643913);
        assert_eq!(inside(&poly, &point), true);
        assert_eq!(inside_precomputed(&poly, &point), true);
        assert_eq!(inside_precomputed_simd(&poly, &point), true);
        assert_eq!(inside_assume_closing(&poly, &point), true);
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
        // v.reverse();
        v
    }

    #[test]
    fn test_circle_parts() {
        // return;
        use rand::prelude::*;
        let mut rng = rand_xorshift::XorShiftRng::seed_from_u64(1);
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
                let expected = inside(&poly, &point);
                // println!("Point: {point:?}");
                assert_eq!(inside_precomputed(&poly, &point), expected);
                assert_eq!(inside_precomputed_simd(&poly, &point), expected);
                assert_eq!(inside_assume_closing(&poly, &point), expected);
            }
        }
    }
}
