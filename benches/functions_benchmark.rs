use bitnuc::{decode, encode, split_packed};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

fn generate_sequence(length: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    (0..length).map(|i| bases[i % 4]).collect()
}

fn naive_split(
    ebuf: &[u64],
    slen: usize,
    pos: usize,
    dbuf: &mut Vec<u8>,
    lbuf: &mut Vec<u64>,
    rbuf: &mut Vec<u64>,
) {
    decode(ebuf, slen, dbuf).unwrap();
    let (l, r) = dbuf.split_at(pos);
    encode(l, lbuf).unwrap();
    encode(r, rbuf).unwrap();
}

fn bench_splitting_naive(c: &mut Criterion) {
    let mut naive_group = c.benchmark_group("naive_split");
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

        naive_group.bench_with_input(BenchmarkId::new(impl_type, size), size, |b, _| {
            b.iter(|| {
                let mut dbuf = Vec::with_capacity(*size);
                let mut lbuf = Vec::new();
                let mut rbuf = Vec::new();
                naive_split(&ebuf, *size, *size / 2, &mut dbuf, &mut lbuf, &mut rbuf);
            });
        });
    }

    naive_group.finish();
}

fn bench_splitting_explicit(c: &mut Criterion) {
    let mut explicit_group = c.benchmark_group("explicit_split");
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

        explicit_group.bench_with_input(BenchmarkId::new(impl_type, size), size, |b, _| {
            b.iter(|| {
                let mut lbuf = Vec::new();
                let mut rbuf = Vec::new();
                split_packed(&ebuf, *size, *size / 2, &mut lbuf, &mut rbuf).unwrap();
            });
        });
    }

    explicit_group.finish();
}

criterion_group!(benches, bench_splitting_naive, bench_splitting_explicit);
criterion_main!(benches);
