# Playground point in polygon simd

Repo holding some explorations aimed at writing the fastest possible point in polygon checks using AVX2.

## Notes

Some notes from explorations so far.

### Direct
Starting point for all algorithms is [Franklin's](https://wrfranklin.org/Research/Short_Notes/pnpoly.html), this is a seemingly simple C function;

```
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
```

The naive rust implementation of the exact same logic is optimised amazingly well by the rustc compiler (1.73), using `-C target-feature=+avx2 -C target-feature=+fma`. When trying to tackle 4 edges in one step (`__mm256d` holds 4 doubles), a lot of performance is lost on loading the interleaved `x` and `y` coordinates, requiring either masked loads and permutates or `_mm256_set_pd`.

So far, `inside_simd` does *not* beat the naive `inside` function because the compiler can optimise it so well.


### Precomputed

If we are allowed to prepare the vertices and precompute data before performing a large number of point in polygon checks, this changes things.

In the original code, there's roughly two conditions;
```
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
```

The first, the `( ((verty[i]>testy) != (verty[j]>testy))` only checks if the current test y is straddled by the edge under consideration. Then the second condition determines whether the test point is left or right of the line.

We can restructure the second condition, by splitting on the `(testy-verty[i])` term;
```
testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]
testx < ((vertx[j]-vertx[i]) * testy / (verty[j]-verty[i])) - (vertx[j]-vertx[i]) * verty[i] / (verty[j]-verty[i]) + vertx[i]
//  Then, introduce an intermediate variable;
let slope = (vertx[j]-vertx[i]) / (verty[j]-verty[i]);
// It becomes;
(testx < (testy * slope) - ( verty[i] * slope + vertx[i])
(testx < (testy * slope -  verty[i] * slope - vertx[i])
// Since the subtraction doesn't depend on the test coordinate, precalculate that;
let sub = -verty[i] * slope - vertx[i];
(testx < (testy * slope + sub)
```

Which leaves `testy * slope + sub`, which can be written with `_mm256_fmadd_pd` in a single SIMD instruction.

Next step is storing all this data in a convenient format that allows for checking 4 edges in one step;
```
struct EdgeVector {
    iy: [f64; 4],
    jy: [f64; 4],
    sub: [f64; 4],
    slope: [f64; 4],
}
pub struct Precomputed {
    edges: Vec<Edge>, // for the tail, could be combined into edge_v.
    edge_v: Vec<EdgeVector>,
}
```

Where the various values are (roughly):
```
while i < vertices.len() {
    let slope = (vertices[j].0 - vertices[i].0) / (vertices[j].1 - vertices[i].1);
    let iy = vertices[i].1;
    let jy = vertices[j].1;
    let sub = -vertices[i].1 * slope + vertices[i].0;
    // Store things.
    j = i;
    i += 1;
}
```

This algorithm still iterates linearly over the vector, but benchmarks do indeed show that we are now iterating through this list of precomputed values four times as fast (7cf6cbe2d2fc8295f2f187a36c5cf46aea49971b). So `Precomputed::inside_simd` is now 4 times as fast as `Precomputed::inside`, which still uses the slope and sub precalculated values.

### Next
Probably a range tree or something, that allows checking only relevant edges, then terminate the nodes there base on the `EdgeVector` datastructure, and put the entire tree in a single vector.

## License
License is `BSD-3-Clause`.