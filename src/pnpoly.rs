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
        let testv = _mm_maskload_pd(std::mem::transmute::<_, *const f64>(&test.0), mask);
        while i < vertices.len() {
            let curr = _mm_maskload_pd(std::mem::transmute::<_, *const f64>(&vertices[i].0), mask);
            let next = _mm_maskload_pd(std::mem::transmute::<_, *const f64>(&vertices[j].0), mask);
            trace!("curr {}", pd(&curr));
            trace!("next {}", pd(&next));

            let a = (vertices[i].1 > test.1);
            let b = (vertices[j].1 > test.1);
            let cmpa = _mm_cmp_pd(curr, testv, _CMP_GE_OQ);
            let cmpb = _mm_cmp_pd(next, testv, _CMP_GE_OQ);
            trace!("cmpa {}, a: {}", pd(&cmpa), a);
            trace!("cmpb {}, b: {}", pd(&cmpb), b);
            // note, only care about upper of cmpa & cmpb
            // let c1 = (test.0 < (vertices[j].0 - vertices[i].0) * (test.1 - vertices[i].1) / (vertices[j].1 - vertices[i].1) + vertices[i].0);
            let c1 = (test.0 - vertices[i].0)
                < (vertices[j].0 - vertices[i].0) * (test.1 - vertices[i].1)
                    / (vertices[j].1 - vertices[i].1);
            let cc = _mm_sub_pd(testv, curr);
            let vdiff = _mm_sub_pd(next, curr);
            let cc_diff = _mm_div_pd(cc, vdiff); // is a division by zero, that's why this all breaks down.
                                                 // c1 = lower < upper (both of cc_diff);
            let upper_at_left = _mm_permute_pd(cc_diff, 0b11);
            let c1_s = _mm_cmp_sd(cc_diff, upper_at_left, _CMP_LE_OQ);
            let c1 = _mm_testc_pd(c1_s, c1_s) != 0;

            trace!("c1_s {}, c1: {}", pd(&c1_s), c1);

            // a != b condition checks if test is between a or b's y coordinate.
            // Then the 'c' condition checks on which side we are of the line.

            if (a != b) && c1 {
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
