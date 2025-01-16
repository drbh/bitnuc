use bitnuc::{as_2bit, encode_alloc, hdist, hdist_scalar};
use criterion::{criterion_group, criterion_main, Criterion};

fn generate_sequence(length: usize, modulo: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    (0..length).map(|i| bases[i % modulo]).collect()
}

#[inline]
fn naive_hamming_distance(s1: &[u8], s2: &[u8]) -> usize {
    s1.iter().zip(s2.iter()).filter(|(a, b)| a != b).count()
}

fn bench_hdist_scalar(c: &mut Criterion) {
    let mut group = c.benchmark_group("hdist_scalar");

    let l = 32;
    let s1 = generate_sequence(l, 4);
    let s2 = generate_sequence(l, 3);

    let expected_diff = naive_hamming_distance(&s1, &s2);

    let es1 = as_2bit(&s1).unwrap();
    let es2 = as_2bit(&s2).unwrap();

    group.bench_function("naive_hdist", |b| {
        b.iter(|| {
            let dist = naive_hamming_distance(&s1, &s2);
            assert_eq!(dist, expected_diff);
        });
    });

    group.bench_function("bitnuc_hdist_scalar", |b| {
        b.iter(|| {
            let dist = hdist_scalar(es1, es2, l).unwrap();
            assert_eq!(dist, expected_diff as u32);
        });
    });

    group.finish();
}

fn bench_hdist_multi(c: &mut Criterion) {
    let mut group = c.benchmark_group("hdist_multi");

    let impl_type = if cfg!(feature = "nosimd") {
        "nosimd"
    } else {
        "simd"
    };

    let l = 512;
    let s1 = generate_sequence(l, 4);
    let s2 = generate_sequence(l, 3);

    let expected_diff = naive_hamming_distance(&s1, &s2);

    let es1 = encode_alloc(&s1).unwrap();
    let es2 = encode_alloc(&s2).unwrap();

    group.bench_function("naive_hdist", |b| {
        b.iter(|| {
            let dist = naive_hamming_distance(&s1, &s2);
            assert_eq!(dist, expected_diff);
        });
    });

    group.bench_function(format!("bitnuc_hdist_multi_{}", impl_type), |b| {
        b.iter(|| {
            let dist = hdist(&es1, &es2, l).unwrap();
            assert_eq!(dist, expected_diff as u32);
        });
    });

    group.finish();
}

criterion_group!(benches, bench_hdist_multi, bench_hdist_scalar);
criterion_main!(benches);
