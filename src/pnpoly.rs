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

// maybe
// https://en.wikipedia.org/wiki/Interval_tree
// https://en.wikipedia.org/wiki/Segment_tree
struct Edge {
    iy: f64,
    jy: f64,
    sub: f64,
    slope: f64,
}
pub struct Minimal {
    edges: Vec<Edge>,
}
impl Minimal {
    pub fn new(vertices: &[(f64, f64)]) -> Self {
        let mut edges: Vec<Edge> = vec![];
        let mut i = 0;
        let mut j = vertices.len() - 1;
        while i < vertices.len() {
            let slope = (vertices[j].0 - vertices[i].0) /  (vertices[j].1 - vertices[i].1);
            edges.push(Edge {
                iy: vertices[i].1,
                jy: vertices[j].1,
                sub: -vertices[i].1 + vertices[i].0 / slope,
                slope,
            });
            j = i;
            i += 1;
        }
        Minimal{edges}
    }
    pub fn inside(self, test: &(f64, f64)) -> bool {
        let mut inside = false;
        for edge in self.edges.iter() {
            if (edge.iy <= test.1) && (test.1 < edge.jy) {
                let c = test.0 < (((test.1 + edge.sub) * edge.slope));
                if c {
                    inside = !inside
                }
            }
        }
        inside
    }
}


pub fn inside_simd(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
    use crate::print::pd;
    use crate::trace;
    use std::arch::x86_64::*;
    unsafe {
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
}

#[cfg(test)]
mod test {
    use super::*;
    // https://github.com/boostorg/geometry/blob/3f5c044abcc813a36d6af83465a9c086f9728a2f/test/strategies/franklin.cpp

    fn inside_minimal(vertices: &[(f64, f64)], test: &(f64, f64)) -> bool {
        let z = Minimal::new(vertices);
        z.inside(test)
    }

    #[test]
    fn test_triangle() {
        for f in [inside, inside_minimal] {
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
        for f in [inside, inside_minimal] {
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
    }

}
