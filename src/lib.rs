pub mod print {
    #[allow(dead_code)]
    pub const DO_PRINTS: bool = false;
    use std::arch::x86_64::{__m256, __m256d, __m256i};
    #[allow(dead_code)]
    /// Print a vector of m128 type.
    pub fn pi(input: &__m256i) -> String {
        use std::arch::x86_64::_mm256_storeu_si256;
        let v: [u8; 32] = [0; 32];
        unsafe { _mm256_storeu_si256(v.as_ptr() as *mut _, *input) };
        format!("{:02X?}", v)
    }
    pub fn pli(input: &__m256i) -> String {
        use std::arch::x86_64::_mm256_storeu_si256;
        let v: [i64; 4] = [0; 4];
        unsafe { _mm256_storeu_si256(v.as_ptr() as *mut _, *input) };
        format!("{:02?}", v)
    }


    #[allow(dead_code)]
    /// Print a vector of m128d type.
    pub fn pf(input: &__m256) -> String {
        use std::arch::x86_64::_mm256_storeu_ps;
        let v: [f32; 8] = [0.0; 8];
        unsafe { _mm256_storeu_ps(v.as_ptr() as *mut _, *input) };
        format!("{:?}", v)
    }
    #[allow(dead_code)]
    /// Print a vector of m128d type.
    pub fn pd(input: &__m256d) -> String {
        use std::arch::x86_64::_mm256_storeu_pd;
        let v: [f64; 4] = [0.0; 4];
        unsafe { _mm256_storeu_pd(v.as_ptr() as *mut _, *input) };
        format!(
            "[{:?} ({:x?}), {:?} ({:x?}), {:?} ({:x?}), {:?} ({:x?})]",
            v[0],
            v[0].to_bits(),
            v[1],
            v[1].to_bits(),
            v[2],
            v[2].to_bits(),
            v[3],
            v[3].to_bits(),
        )
    }
}

#[allow(unused_macros)]
/// Helper print macro that can be enabled or disabled.
macro_rules! trace {
    () => (if crate::print::DO_PRINTS {println!("\n");});
    ($($arg:tt)*) => {
        if crate::print::DO_PRINTS {
            println!($($arg)*);
        }
    }
}

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

/// Simplify the start / end condition since we always have a closing point.
pub fn inside_assume_closed(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
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

/// Vectorized form without precomputation. Assumes closed.
pub fn inside_simd(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
    unsafe {
        #[allow(unused_imports)]
        use print::*;
        use std::arch::x86_64::*;
        let mut i = 0;

        // Step in fours;
        let ty = _mm256_set1_pd(test.1);
        let tx = _mm256_set1_pd(test.0);

        let mut crossings_totals = _mm256_set_epi64x(0, 0, 0, 0);
        while i + (4 + 1) <= vertices.len() {
            let base_x01 = std::mem::transmute::<_, *const f64>(&vertices[i].0);
            let base_y01 = std::mem::transmute::<_, *const f64>(&vertices[i].1);
            let base_x23 = std::mem::transmute::<_, *const f64>(&vertices[i + 1].1);
            let base_y23 = std::mem::transmute::<_, *const f64>(&vertices[i + 2].0);

            let base_x4 = std::mem::transmute::<_, *const f64>(&vertices[i + 4].0);
            let base_y4 = std::mem::transmute::<_, *const f64>(&vertices[i + 4].1);

            // vertices;
            // x0, y0, x1, y1| x2, y2, x3, y3| x4, y4
            // need:
            // x0, x1, x2, x3
            //     x1, x2, x3, x4
            let ix;
            let iy;
            let jx;
            let jy;
            const USE_MOVES: bool = true;
            if USE_MOVES {
                let mask_1010 = _mm256_set_epi64x(0, -1, 0, -1);
                let mask_0101 = _mm256_set_epi64x(-1, 0, -1, 0);
                let mask_0001 = _mm256_set_epi64x(0, 0, 0, -1);

                let f64_bitsset = f64::from_ne_bytes(0xFFFFFFFF_FFFFFFFFu64.to_ne_bytes());
                let mask_1110 = _mm256_set_pd(f64_bitsset, f64_bitsset, f64_bitsset, 0.0);
                // Order doesn't actually matter, as long as it is consistent.
                // __m256d _mm256_maskload_pd (double const * mem_addr, __m256i mask)
                let ix_1010 = _mm256_maskload_pd(base_x01, mask_1010); // holds x0, 0, x1, 0
                let iy_1010 = _mm256_maskload_pd(base_y01, mask_1010); // holds y0, 0, y1, 0
                let ix_0101 = _mm256_maskload_pd(base_x23, mask_0101); // holds 0, x2, 0, x3
                let iy_0101 = _mm256_maskload_pd(base_y23, mask_0101); // holds 0, y2, 0, y3
                                                                       // trace!("iy_1010: {}", pd(&iy_1010));
                                                                       // trace!("iy_0101: {}", pd(&iy_0101));

                ix = _mm256_or_pd(ix_1010, ix_0101);
                iy = _mm256_or_pd(iy_1010, iy_0101);

                // These vectors are      x0, x2, x1, x3
                // To get the +1, we need to make:
                // here, we need to make; x1, x3, x2, x4
                // So that's shuffle the lanes.
                // First, replace 0 with 4, then we can do a permutate to cross lane boundaries.
                let jx_0001 = _mm256_maskload_pd(base_x4, mask_0001); // holds x4, 0, 0, 0
                let jy_0001 = _mm256_maskload_pd(base_y4, mask_0001); // holds y4, 0, 0, 0

                let ix_1110 = _mm256_and_pd(ix, mask_1110);
                let iy_1110 = _mm256_and_pd(iy, mask_1110);

                let jx_wrong = _mm256_or_pd(jx_0001, ix_1110);
                let jy_wrong = _mm256_or_pd(jy_0001, iy_1110);
                // trace!("jx_wrong: {}", pd(&jx_wrong));
                // trace!("jy_wrong: {}", pd(&jy_wrong));
                // Now, all that remains is a permutate to cross some lane boundaries.
                //                        11  10  01  00
                // holds                  x4, x2, x1, x3
                // here, we need to make; x1, x3, x2, x4
                jx = _mm256_permute4x64_pd(jx_wrong, 0b00_01_11_10);
                jy = _mm256_permute4x64_pd(jy_wrong, 0b00_01_11_10);
            } else {
                ix = _mm256_set_pd(vertices[i+3].0, vertices[i+2].0, vertices[i+1].0, vertices[i].0);
                jx = _mm256_set_pd(vertices[i+3 + 1].0, vertices[i+2 + 1].0, vertices[i+1 + 1].0, vertices[i + 1].0);
                iy = _mm256_set_pd(vertices[i+3].1, vertices[i+2].1, vertices[i+1].1, vertices[i].1);
                jy = _mm256_set_pd(vertices[i+3 + 1].1, vertices[i+2 + 1].1, vertices[i+1 + 1].1, vertices[i + 1].1);
            }

            /*
                if ((vertices[i].1 > test.1) != (vertices[j].1 > test.1))
                    && (test.0
                        < (vertices[j].0 - vertices[i].0) * (test.1 - vertices[i].1)
                            / (vertices[j].1 - vertices[i].1)
                            + vertices[i].0)
            */

            //  edge.iy <= test.1
            let above_lower = _mm256_cmp_pd(iy, ty, _CMP_LE_OQ);
            // trace!("above_lower {}", pd(&above_lower));
            // test.1 <= edge.jy
            let below_upper = _mm256_cmp_pd(ty, jy, _CMP_LT_OQ);
            // trace!("below_upper {}", pd(&below_upper));

            // Or the other way around.
            // edge.iy <= test.1 <= edge.jy
            //
            let inverted_coords_above_lower = _mm256_cmp_pd(jy, ty, _CMP_LE_OQ);
            let inverted_coords_below_upper = _mm256_cmp_pd(ty, iy, _CMP_LT_OQ);

            // let c = test.0 < (((test.1 * edge.slope + edge.sub) ));

            // let right = _mm256_fmadd_pd(ty, slope, sub);
            /*
                (vertices[j].0 - vertices[i].0) *
                    (test.1 - vertices[i].1) / (vertices[j].1 - vertices[i].1) + vertices[i].0)
                a * b / c + z
            */

            let a = _mm256_sub_pd(jx, ix);
            let b = _mm256_sub_pd(ty, iy);
            let c = _mm256_sub_pd(jy, iy);
            let z = ix;
            let bc = _mm256_div_pd(b, c);
            let right = _mm256_fmadd_pd(a, bc, z);

            let t_l_right = _mm256_cmp_pd(tx, right, _CMP_LT_OQ);
            // trace!("t_l_right {}", pd(&t_l_right));
            // println!("jskldfjsd");

            // Now, we mask all three together.
            let in_range = _mm256_and_pd(above_lower, below_upper);
            let inverted_in_range =
                _mm256_and_pd(inverted_coords_above_lower, inverted_coords_below_upper);
            let in_range = _mm256_or_pd(in_range, inverted_in_range);
            let crosses = _mm256_and_pd(in_range, t_l_right);

            // Cast this to integers, which gets is 0 for not true, -1 for true;
            let crossess_i = _mm256_castpd_si256(crosses);
            crossings_totals = _mm256_sub_epi64(crossings_totals, crossess_i);
            i += 4;
        }

        // Finish the tail if not a multiple of four.
        let mut cross_normal = [0i64; 4];
        _mm256_storeu_si256(std::mem::transmute::<_, *mut __m256i>(&cross_normal[0]), crossings_totals);


        while i < vertices.len() - 1 {
            let j = i + 1;
            if ((vertices[i].1 > test.1) != (vertices[j].1 > test.1))
                && (test.0
                    < (vertices[j].0 - vertices[i].0) * (test.1 - vertices[i].1)
                        / (vertices[j].1 - vertices[i].1)
                        + vertices[i].0)
            {
                cross_normal[0] -= 1;
            }
            i += 1;
        }

        let t: i64 = cross_normal.iter().sum();
        t.abs() % 2 == 1
    }
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

// Instead of 4 indepent vectors, we still interleave the individual chunks, that way the
// values we are loading into simd vectors are likely in the cache after the first iy read.
#[repr(C)]
#[derive(Default)]
struct EdgeVector {
    iy: [f64; 4],
    jy: [f64; 4],
    sub: [f64; 4],
    slope: [f64; 4],
}

/// Precomputed version of Franklin
pub struct Precomputed {
    edges: Vec<Edge>,
    edge_v: Vec<EdgeVector>,
}
impl Precomputed {
    pub fn new(vertices: &[(f64, f64)]) -> Self {
        let mut edges: Vec<Edge> = vec![];
        let mut edge_v: Vec<EdgeVector> = vec![];
        edge_v.push(Default::default());

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
            edge_v[i / 4].iy[i % 4] = lower_y;
            edge_v[i / 4].jy[i % 4] = upper_y;
            edge_v[i / 4].sub[i % 4] = -vertices[i].1 * slope + vertices[i].0;
            edge_v[i / 4].slope[i % 4] = slope;
            j = i;
            i += 1;
            if i % 4 == 0 {
                edge_v.push(Default::default());
            }
        }
    

        Self { edges, edge_v }
    }

    /// Simple inside check, iterating through the edges.
    pub fn inside(&self, test: &(f64, f64)) -> bool {
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
    pub fn inside_simd(&self, test: &(f64, f64)) -> bool {
        // use print::pd;
        // use trace;
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
            let mut i = 0;

            // Step in fours;
            let ty = _mm256_set1_pd(test.1);
            let tx = _mm256_set1_pd(test.0);
            let mut crossings_totals = _mm256_set_epi64x(0, 0, 0, 0);
            while i + 4 <= self.edges.len() {
                let e = &self.edge_v[i / 4];
                // Could do a gather here, but lets get this working and then just make the struct
                // axis first instead of vertex first.
                let iy = _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&e.iy[0]));
                let jy = _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&e.jy[0]));
                let sub =  _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&e.sub[0]));
                let slope =  _mm256_loadu_pd(std::mem::transmute::<_, *const f64>(&e.slope[0]));

                //  edge.iy <= test.1
                let above_lower = _mm256_cmp_pd(iy, ty, _CMP_LE_OQ);
                // trace!("above_lower {}", pd(&above_lower));
                // test.1 <= edge.jy
                let below_upper = _mm256_cmp_pd(ty, jy, _CMP_LT_OQ);
                // trace!("below_upper {}", pd(&below_upper));

                // Or the other way around.
                // edge.iy <= test.1 <= edge.jy
                //
                let inverted_coords_above_lower = _mm256_cmp_pd(jy, ty, _CMP_LE_OQ);
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

                // Cast this to integers, which gets is 0 for not true, -1 for true;
                let crossess_i = _mm256_castpd_si256(crosses);
                crossings_totals = _mm256_sub_epi64(crossings_totals, crossess_i);

                i += 4;
            }

            // Finish the tail if not a multiple of four.
            let mut cross_normal = [0i64; 4];
            _mm256_storeu_si256(std::mem::transmute::<_, *mut __m256i>(&cross_normal[0]), crossings_totals);

            while i < self.edges.len() {
                let edge = &self.edges[i];
                if (edge.iy <= test.1) && (test.1 < edge.jy) {
                    let c = test.0 < (test.1 * edge.slope + edge.sub);
                    if c {
                        cross_normal[0] -= 1;
                    }
                }
                i += 1;
            }

            let t: i64 = cross_normal.iter().sum();
            t.abs() % 2 == 1
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    // Square and triangle from:
    // https://github.com/boostorg/geometry/blob/3f5c044abcc813a36d6af83465a9c086f9728a2f/test/strategies/franklin.cpp

    // https://www.desmos.com/calculator is helpful, rust debug prints can be pasted there.

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
            inside_assume_closed,
            inside_simd,
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
            inside_assume_closed,
            inside_simd,
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
        assert_eq!(inside_assume_closed(&poly, &point), false);
        assert_eq!(inside_simd(&poly, &point), false);
        let point = (12.230435569131824, 3.1868743396643913);
        assert_eq!(inside(&poly, &point), true);
        assert_eq!(inside_precomputed(&poly, &point), true);
        assert_eq!(inside_precomputed_simd(&poly, &point), true);
        assert_eq!(inside_assume_closed(&poly, &point), true);
        assert_eq!(inside_simd(&poly, &point), true);
    }

    #[test]
    fn test_polygon_increasing() {
        // Nice properties of the vertices incrementing.
        let polygon = vec![
            (1.0, 2.0),
            (3.0, 4.0),
            (5.0, 6.0),
            (7.0, 8.0),
            (9.0, 10.0),
            (0.0, 10.0),
            (0.0, 1.0),
            (1.0, 2.0),
        ];
        let test_points = [
            ((2.0, 4.0), true),
            ((0.0, 0.0), false),
            ((4.0, 8.0), true),
            ((6.0, 0.0), false),
        ];
        for f in [
            inside,
            inside_precomputed,
            inside_precomputed_simd,
            inside_assume_closed,
            inside_simd,
        ] {
            for (p, expected) in test_points.iter() {
                assert_eq!(f(&polygon, p), *expected);
            }
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
    fn test_circle_parts() {
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
                let expected = inside(&poly, &point);
                // println!("Point: {point:?}");
                assert_eq!(inside_precomputed(&poly, &point), expected);
                assert_eq!(inside_precomputed_simd(&poly, &point), expected);
                assert_eq!(inside_assume_closed(&poly, &point), expected);
                assert_eq!(inside_simd(&poly, &point), expected);
            }
        }
    }
}
