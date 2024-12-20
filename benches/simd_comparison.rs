use bitnuc::PackedSequence;
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
            |b, seq| b.iter(|| PackedSequence::new(seq)),
        );
    }

    group.finish();
}

criterion_group!(benches, bench_simd_comparison);
criterion_main!(benches);
