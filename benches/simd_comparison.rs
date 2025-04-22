use bitnuc::{as_2bit, decode, encode_alloc, from_2bit};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

fn generate_sequence(length: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    (0..length).map(|i| bases[i % 4]).collect()
}

fn bench_packing(c: &mut Criterion) {
    let mut group = c.benchmark_group("packing");

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

fn bench_encoding(c: &mut Criterion) {
    let mut group = c.benchmark_group("packing");

    let impl_type = if cfg!(feature = "nosimd") {
        "nosimd"
    } else {
        "simd"
    };

    // Test different sequence lengths
    for size in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024].iter() {
        let seq = generate_sequence(*size);

        group.bench_with_input(
            BenchmarkId::new(format!("encoding_{}", impl_type), size),
            &seq,
            |b, seq| b.iter(|| encode_alloc(seq).unwrap()),
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

fn bench_decoding(c: &mut Criterion) {
    let mut group = c.benchmark_group("unpacking");

    let impl_type = if cfg!(feature = "nosimd") {
        "no_simd"
    } else {
        "simd"
    };

    // Test different sequence lengths
    for size in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024].iter() {
        let seq = generate_sequence(*size);
        let packed = encode_alloc(&seq).unwrap();

        group.bench_with_input(
            BenchmarkId::new(format!("decoding_{}", impl_type), size),
            &packed,
            |b, packed| b.iter(|| decode(packed, *size, &mut Vec::new()).unwrap()),
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_packing,
    bench_encoding,
    bench_unpacking,
    bench_decoding
);
criterion_main!(benches);
