// benches/packing_benchmarks.rs
use bitnuc::PackedSequence;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

fn generate_sequence(length: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    (0..length).map(|i| bases[i % 4]).collect()
}

fn bench_packing(c: &mut Criterion) {
    let mut group = c.benchmark_group("packing");

    // Test different sequence lengths
    for size in [4, 8, 16, 24, 32].iter() {
        let seq = generate_sequence(*size);

        group.bench_with_input(BenchmarkId::new("pack", size), &seq, |b, seq| {
            b.iter(|| PackedSequence::new(black_box(seq)))
        });
    }

    group.finish();
}

fn bench_unpacking(c: &mut Criterion) {
    let mut group = c.benchmark_group("unpacking");

    // Test different sequence lengths
    for size in [4, 8, 16, 24, 32].iter() {
        let seq = generate_sequence(*size);
        let packed = PackedSequence::new(&seq).unwrap();

        group.bench_with_input(BenchmarkId::new("unpack", size), &packed, |b, packed| {
            b.iter(|| black_box(packed).to_vec())
        });
    }

    group.finish();
}

fn bench_roundtrip(c: &mut Criterion) {
    let mut group = c.benchmark_group("roundtrip");

    // Test different sequence lengths
    for size in [4, 8, 16, 24, 32].iter() {
        let seq = generate_sequence(*size);

        group.bench_with_input(BenchmarkId::new("pack_unpack", size), &seq, |b, seq| {
            b.iter(|| {
                let packed = PackedSequence::new(black_box(seq)).unwrap();
                black_box(packed).to_vec()
            })
        });
    }

    group.finish();
}

fn bench_sequence_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequence_patterns");

    // Test different sequence patterns
    let patterns = vec![
        ("homopolymer", vec![b'A'; 32]),
        (
            "alternating",
            (0..32)
                .map(|i| if i % 2 == 0 { b'A' } else { b'T' })
                .collect(),
        ),
        (
            "repeating_4",
            b"ACGT".iter().cycle().take(32).cloned().collect(),
        ),
        (
            "gc_rich",
            b"GCGC".iter().cycle().take(32).cloned().collect(),
        ),
        (
            "at_rich",
            b"ATAT".iter().cycle().take(32).cloned().collect(),
        ),
    ];

    for (name, seq) in patterns {
        group.bench_with_input(BenchmarkId::new("pattern", name), &seq, |b, seq| {
            b.iter(|| {
                let packed = PackedSequence::new(black_box(seq)).unwrap();
                black_box(packed).to_vec()
            })
        });
    }

    group.finish();
}

fn bench_access_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("access_patterns");

    let seq = PackedSequence::new(&generate_sequence(32)).unwrap();

    group.bench_function("sequential_access", |b| {
        b.iter(|| {
            for i in 0..32 {
                black_box(seq.get(i)).unwrap();
            }
        })
    });

    group.bench_function("random_access", |b| {
        b.iter(|| {
            // Access some random positions (but deterministic for benchmark)
            for &i in &[7, 31, 0, 15, 23, 4, 12, 28] {
                black_box(seq.get(i)).unwrap();
            }
        })
    });

    group.bench_function("slice_small", |b| {
        b.iter(|| black_box(seq.slice(0..8)).unwrap())
    });

    group.bench_function("slice_large", |b| {
        b.iter(|| black_box(seq.slice(0..24)).unwrap())
    });

    group.finish();
}

fn bench_fast_packing_long_sequences(c: &mut Criterion) {
    let mut group = c.benchmark_group("packing");

    let impl_type = if cfg!(feature = "nosimd") {
        "nosimd"
    } else {
        "simd"
    };

    // Test different sequence lengths
    for size in [32_000, 64_000, 128_000, 256_000, 512_000].iter() {
        let seq = generate_sequence(*size);

        group.bench_with_input(
            BenchmarkId::new(format!("pack_{}", impl_type), size),
            &seq,
            |b, seq| b.iter(|| PackedSequence::new(black_box(seq))),
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_packing,
    bench_unpacking,
    bench_roundtrip,
    bench_sequence_patterns,
    bench_access_patterns,
    bench_fast_packing_long_sequences
);
criterion_main!(benches);
