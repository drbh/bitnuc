use bitnuc::{as_2bit, from_2bit};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

fn generate_sequence(length: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    (0..length).map(|i| bases[i % 4]).collect()
}

fn bench_simd_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("simd_comparison");

    let impl_type = if cfg!(feature = "nosimd") {
        "nosimd"
    } else {
        "simd"
    };

    // Test different sequence lengths
    for size in [4, 8, 16, 24, 32].iter() {
        let seq = generate_sequence(*size);

        group.bench_with_input(
            BenchmarkId::new(format!("packing_{}", impl_type), size),
            &seq,
            |b, seq| b.iter(|| as_2bit(seq)),
        );
    }

    group.finish();
}

fn bench_unpacking(c: &mut Criterion) {
    let mut group = c.benchmark_group("unpacking");

    let impl_type = if cfg!(feature = "nosimd") {
        "no_simd"
    } else {
        "simd"
    };

    // Test different sequence lengths
    for size in [4, 8, 16, 24, 32].iter() {
        let seq = generate_sequence(*size);
        let packed = as_2bit(&seq).unwrap();

        group.bench_with_input(
            BenchmarkId::new(format!("unpacking_{}", impl_type), size),
            &packed,
            |b, packed| b.iter(|| from_2bit(*packed, *size, &mut Vec::new()).unwrap()),
        );
    }

    group.finish();
}

criterion_group!(benches, bench_simd_comparison, bench_unpacking);
criterion_main!(benches);
