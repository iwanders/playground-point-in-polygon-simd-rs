

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

pub fn inside(vertices: &[(f32, f32)], test: &(f32, f32)) -> bool {
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

#[cfg(test)]
mod test{
    use super::*;
    // https://github.com/boostorg/geometry/blob/3f5c044abcc813a36d6af83465a9c086f9728a2f/test/strategies/franklin.cpp

    #[test]
    fn test_triangle() {
        let triangle = vec![(0.0, 0.0), (0.0,  4.0), (6.0,  0.0), (0.0, 0.0)];
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
    fn test_square() {
        let square = vec![(0.0, 0.0), (0.0,  2.0), (2.0, 2.0), (2.0,  0.0), (0.0, 0.0)];
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
}