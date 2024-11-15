use bitnuc::{GCContent, PackedSequence};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

pub fn sequence_benchmark(c: &mut Criterion) {
    c.bench_function("pack sequence", |b| {
        b.iter(|| PackedSequence::new(black_box(b"ACGTACGTACGTACGT")).unwrap())
    });

    let seq = PackedSequence::new(b"ACGTACGTACGTACGT").unwrap();

    c.bench_function("gc content", |b| b.iter(|| black_box(&seq).gc_content()));
}

criterion_group!(benches, sequence_benchmark);
criterion_main!(benches);
