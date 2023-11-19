use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use point_in_polygon::{inside, inside_simd};

use rand::prelude::*;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

fn make_polygon(points: usize, seed: u64, radius: f64) -> Vec<(f64, f64)> {
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
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);

    let mut distances: Vec<(f64, f64)> = vec![];
    let r = rng.gen::<f64>() * radius;
    let mut s = 0.0;
    distances.push((s, r * rng.gen::<f64>()));
    let segments: usize = points;
    for _ in 0..segments {
        let sa = (1.0 / (segments as f64)) * rng.gen::<f64>();
        s += sa;
        if s > 1.0 {
            break;
        }
        distances.push((s, r * rng.gen::<f64>()));
    }
    create_circle_parts(&distances)
}

fn make_points(points: usize, seed: u64, radius: f64) -> Vec<(f64, f64)> {
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut r = vec![];
    for _ in 0..points {
        let x = radius * 1.25 * rng.gen::<f64>();
        let y = radius * 1.25 * rng.gen::<f64>();
        r.push((x, y));
    }
    r
}

#[allow(dead_code)]
fn z() {
    let r = 100.0;
    let ppoints = 1000;
    let polygon = make_polygon(ppoints, 1, r);
    let points = make_points(10000, 2, r);
    for _ in 0..10 {
        for point in points.iter() {
            let r = inside_simd(&polygon, point);
            black_box(r);
        }
    }
}

pub fn criterion_benchmark(c: &mut Criterion) {
    // return z();
    let mut group = c.benchmark_group("poly_size");
    // for poly_points in [10, 100, 500, 1000] {
    for poly_points in [500] {
        group.bench_with_input(
            BenchmarkId::new("inside", poly_points),
            &poly_points,
            |b, &ppoints| {
                let r = 100.0;
                let polygon = make_polygon(ppoints, 1, r);
                let points = make_points(10000, 2, r);
                b.iter(|| {
                    for point in points.iter() {
                        let r = inside(&polygon, point);
                        black_box(r);
                    }
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("inside_simd", poly_points),
            &poly_points,
            |b, &ppoints| {
                let r = 100.0;
                let polygon = make_polygon(ppoints, 1, r);
                let points = make_points(10000, 2, r);
                b.iter(|| {
                    for point in points.iter() {
                        let r = inside_simd(&polygon, point);
                        black_box(r);
                    }
                });
            },
        );
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
