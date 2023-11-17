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

pub fn inside(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
    let mut i = 0;
    let mut j = vertices.len() - 1;
    let mut inside = false;
    while i < vertices.len() {
        if ((vertices[i].1 > test.1) != (vertices[j].1 > test.1)) &&
            // (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
            (test.0 < (vertices[j].0 - vertices[i].0) * (test.1 - vertices[i].1) / (vertices[j].1 - vertices[i].1) + vertices[i].0)
        {
            inside = !inside;
        }
        j = i;
        i += 1;
    }
    inside
}

pub fn inside_simd(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
    use crate::print::pd;
    use crate::trace;
    use std::arch::x86_64::*;
    unsafe {
        let mut i = 0;
        let mut j = vertices.len() - 1;
        // __m128d
        let mut inside = false;
        let mask = _mm_set1_epi8(-1);

        let nans_v = [f64::NAN, f64::NAN];
        let testv = _mm_maskload_pd(std::mem::transmute::<_, *const f64>(&test.0), mask);
        while i < vertices.len() {
            trace!();
            let curr = _mm_maskload_pd(std::mem::transmute::<_, *const f64>(&vertices[i].0), mask);
            let next = _mm_maskload_pd(std::mem::transmute::<_, *const f64>(&vertices[j].0), mask);
            let nanv = _mm_maskload_pd(std::mem::transmute::<_, *const f64>(&nans_v[0]), mask);
            trace!("testv {}", pd(&testv));
            trace!("curr {}", pd(&curr));
            trace!("next {}", pd(&next));

            let a = (vertices[i].1 > test.1);
            let b = (vertices[j].1 > test.1);
            let cmpa = _mm_cmp_pd(curr, testv, _CMP_GE_OQ);
            let cmpb = _mm_cmp_pd(next, testv, _CMP_GE_OQ);
            trace!("cmpa {}, a: {}", pd(&cmpa), a);
            trace!("cmpb {}, b: {}", pd(&cmpb), b);
            // note, only care about upper of cmpa & cmpb

            // The harder part is the actual comparison that follows:

            // let c1 = (test.0 < (vertices[j].0 - vertices[i].0) * (test.1 - vertices[i].1) / (vertices[j].1 - vertices[i].1) + vertices[i].0);
            let c1 = (test.0 - vertices[i].0)
                < (vertices[j].0 - vertices[i].0) * (test.1 - vertices[i].1)
                    / (vertices[j].1 - vertices[i].1);

            // Split this;

            let c1_left = (test.0 - vertices[i].0);
            let c1_right = (vertices[j].0 - vertices[i].0) * (test.1 - vertices[i].1) / (vertices[j].1 - vertices[i].1);
            trace!("c1_left {}, c1_right: {}", c1_left, c1_right);

            let c1_l_m = (test.0 - vertices[i].0) / (vertices[j].0 - vertices[i].0);
            let c1_r_m = (test.1 - vertices[i].1) / (vertices[j].1 - vertices[i].1);

            trace!("c1_l_m {}, c1_r_m: {}", c1_l_m, c1_r_m);
            let c1_from_split = c1_l_m < c1_r_m;
            trace!("c1_l_m < c1_r_m: {}", c1_from_split);
            trace!("c1: {}", c1);
            if (c1 != c1_from_split) {
                panic!("c1 and c1_from_split don't agree.");
            }

            let t_min_curr = _mm_sub_pd(testv, curr);
            trace!("t_min_curr : {}", pd(&t_min_curr));

            let v_diff = _mm_sub_pd(next, curr);
            trace!("v_diff : {}", pd(&v_diff));
            let cc_diff = _mm_div_pd(t_min_curr, v_diff); // is a division by zero, that's why this all breaks down.
                                                 // c1 = lower < upper (both of cc_diff);
            trace!("cc_diff : {}", pd(&cc_diff));
            // Compare left with right.

            // let cc_cmp_nan = _mm_cmp_sd(cc_diff, nanv, _CMP_UNORD_Q);
            // let cc_is_nan = _mm_testc_pd(cc_cmp_nan, cc_cmp_nan);

            let upper_at_both = _mm_permute_pd(cc_diff, 0b11);
            let lower_at_both = _mm_permute_pd(cc_diff, 0b00);
            trace!("upper_at_left : {}", pd(&upper_at_both));
            trace!("lower_at_both : {}", pd(&lower_at_both));
            let c1_s = _mm_cmp_sd(lower_at_both, upper_at_both, _CMP_LT_OQ);
            trace!("c1_s : {}", pd(&c1_s)); // only care about lower in this one.
            let c1_at_both = _mm_permute_pd(c1_s, 0b00);
            trace!("c1_at_both : {}", pd(&c1_at_both));

            // use std::arch::x86_64::_mm_storeu_pd;
            // let v: [f64; 2] = [0.0; 2];
            // unsafe { _mm_storeu_pd(v.as_ptr() as *mut _, c1_at_both) };
            let c1_f = (_mm_movemask_pd(c1_s) & 0b10) != 0;
            // let c1_f = _mm_testc_pd(c1_at_both, mask) != 0;
            // let c1_f = v[0] == 0.0;
            trace!("c1_f : {}, c1: {}", c1_f, c1);


            if ((a != b) && (c1_f != c1)) {
                panic!();
            }

            // a != b condition checks if test is between a or b's y coordinate.
            // Then the 'c' condition checks on which side we are of the line.

            if (a != b) && c1_f {
                inside = !inside;
            }
            j = i;
            i += 1;
        }
        inside
    }
}

#[cfg(test)]
mod test {
    use super::*;
    // https://github.com/boostorg/geometry/blob/3f5c044abcc813a36d6af83465a9c086f9728a2f/test/strategies/franklin.cpp

    #[test]
    fn test_triangle() {
        let triangle = vec![(0.0, 0.0), (0.0, 4.0), (6.0, 0.0), (0.0, 0.0)];
        assert_eq!(inside(&triangle, &(1.0, 1.0)), true);
        assert_eq!(inside(&triangle, &(3.0, 3.0)), false);

        assert_eq!(inside(&triangle, &(0.0, 0.0)), true);
        assert_eq!(inside(&triangle, &(0.0, 4.0)), false);
        assert_eq!(inside(&triangle, &(5.0, 0.0)), true);

        assert_eq!(inside(&triangle, &(0.0, 2.0)), true);
        assert_eq!(inside(&triangle, &(3.0, 2.0)), false);
        assert_eq!(inside(&triangle, &(3.0, 0.0)), true);
    }

    #[test]
    fn test_triangle_simd() {
        let triangle = vec![(0.0, 0.0), (0.0, 4.0), (6.0, 0.0), (0.0, 0.0)];
        assert_eq!(inside_simd(&triangle, &(1.0, 1.0)), true);
        assert_eq!(inside_simd(&triangle, &(3.0, 3.0)), false);

        assert_eq!(inside_simd(&triangle, &(0.0, 0.0)), true);
        assert_eq!(inside_simd(&triangle, &(0.0, 4.0)), false);
        assert_eq!(inside_simd(&triangle, &(5.0, 0.0)), true);

        assert_eq!(inside_simd(&triangle, &(0.0, 2.0)), true);
        assert_eq!(inside_simd(&triangle, &(3.0, 2.0)), false);
        assert_eq!(inside_simd(&triangle, &(3.0, 0.0)), true);
    }

    #[test]
    fn test_square() {
        let square = vec![(0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0), (0.0, 0.0)];
        assert_eq!(inside(&square, &(1.0, 1.0)), true);
        assert_eq!(inside(&square, &(3.0, 3.0)), false);

        assert_eq!(inside(&square, &(0.0, 0.0)), true);
        assert_eq!(inside(&square, &(0.0, 2.0)), false);
        assert_eq!(inside(&square, &(2.0, 2.0)), false);
        assert_eq!(inside(&square, &(2.0, 0.0)), false);

        assert_eq!(inside(&square, &(0.0, 1.0)), true);
        assert_eq!(inside(&square, &(1.0, 2.0)), false);
        assert_eq!(inside(&square, &(2.0, 1.0)), false);
        assert_eq!(inside(&square, &(1.0, 0.0)), true);
    }

    #[test]
    fn test_square_simd() {
        let square = vec![(0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0), (0.0, 0.0)];
        assert_eq!(inside_simd(&square, &(1.0, 1.0)), true);
        assert_eq!(inside_simd(&square, &(3.0, 3.0)), false);

        assert_eq!(inside_simd(&square, &(0.0, 0.0)), true);
        assert_eq!(inside_simd(&square, &(0.0, 2.0)), false);
        assert_eq!(inside_simd(&square, &(2.0, 2.0)), false);
        assert_eq!(inside_simd(&square, &(2.0, 0.0)), false);

        assert_eq!(inside_simd(&square, &(0.0, 1.0)), true);
        assert_eq!(inside_simd(&square, &(1.0, 2.0)), false);
        assert_eq!(inside_simd(&square, &(2.0, 1.0)), false);
        assert_eq!(inside_simd(&square, &(1.0, 0.0)), true);
    }
}
