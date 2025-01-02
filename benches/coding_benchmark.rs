use bitnuc::{decode, encode};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

fn generate_sequence(length: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    (0..length).map(|i| bases[i % 4]).collect()
}

fn bench_decode(c: &mut Criterion) {
    let mut group = c.benchmark_group("decode");

    let impl_type = if cfg!(feature = "nosimd") {
        "nosimd"
    } else {
        "simd"
    };

    // test different sequence lengths
    for size in [30, 80, 140, 280].iter() {
        let seq = generate_sequence(*size);
        let mut ebuf = Vec::new();
        encode(&seq, &mut ebuf).unwrap();

        group.bench_with_input(BenchmarkId::new(impl_type, size), size, |b, _| {
            b.iter(|| {
                let mut dbuf = Vec::with_capacity(*size);
                decode(&ebuf, *size, &mut dbuf).unwrap();
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_decode);
criterion_main!(benches);
