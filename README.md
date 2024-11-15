# bitnuc

[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE.md)
![actions status](https://github.com/noamteyssier/bitnuc/workflows/Rust/badge.svg)
[![Crates.io](https://img.shields.io/crates/d/bitnuc?color=orange&label=crates.io)](https://crates.io/crates/bitnuc)
[![docs.rs](https://img.shields.io/docsrs/bitnuc?color=green&label=docs.rs)](https://docs.rs/bitnuc/latest/bitnuc/)

A library for efficient nucleotide sequence manipulation using 2-bit encoding.

## Features

- 2-bit nucleotide encoding (A=00, C=01, G=10, T=11)
- Direct bit manipulation functions for custom implementations
- Higher-level sequence type with additional analysis features

## Low-Level Packing Functions

For direct bit manipulation, use the `as_2bit` and `from_2bit` functions:

```rust
use bitnuc::{as_2bit, from_2bit};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Pack a sequence into a u64
    let packed = as_2bit(b"ACGT")?;
    assert_eq!(packed, 0b11100100);

    // Unpack back to a sequence
    let unpacked = from_2bit(packed, 4)?;
    assert_eq!(&unpacked, b"ACGT");
    Ok(())
}
```

These functions are useful when you need to:
- Implement custom sequence storage
- Manipulate sequences at the bit level
- Integrate with other bioinformatics tools
- Copy sequences more efficiently
- Hash sequences more efficiently

For example, packing multiple short sequences:

```rust
use bitnuc::{as_2bit, from_2bit};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Pack multiple 4-mers into u64s
    let kmers = [b"ACGT", b"TGCA", b"GGCC"];
    let packed: Vec<u64> = kmers
        .into_iter()
        .map(|kmer| as_2bit(kmer))
        .collect::<Result<_, _>>()?;

    // Unpack when needed
    let first_kmer = from_2bit(packed[0], 4)?;
    assert_eq!(&first_kmer, b"ACGT");
    Ok(())
}
```

## High-Level Sequence Type

For more complex sequence manipulation, use the [`PackedSequence`] type:

```rust
use bitnuc::{PackedSequence, GCContent, BaseCount};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let seq = PackedSequence::new(b"ACGTACGT")?;

    // Sequence analysis
    println!("GC Content: {}%", seq.gc_content());
    let [a_count, c_count, g_count, t_count] = seq.base_counts();

    // Slicing
    let subseq = seq.slice(1..5)?;
    assert_eq!(&subseq, b"CGTA");
    Ok(())
}
```

## Memory Usage

The 2-bit encoding provides significant memory savings:

```text
Standard encoding: 1 byte per base
ACGT = 4 bytes = 32 bits

2-bit encoding: 2 bits per base
ACGT = 8 bits
```

This means you can store 4 times as many sequences in the same amount of memory.

## Error Handling

All operations that could fail return a [`Result`] with [`NucleotideError`]:

```rust
use bitnuc::{as_2bit, NucleotideError};

// Invalid nucleotide
let err = as_2bit(b"ACGN").unwrap_err();
assert!(matches!(err, NucleotideError::InvalidBase(b'N')));

// Sequence too long
let long_seq = vec![b'A'; 33];
let err = as_2bit(&long_seq).unwrap_err();
assert!(matches!(err, NucleotideError::SequenceTooLong(33)));
```

## Performance Considerations

When working with many short sequences (like k-mers), using `as_2bit` and `from_2bit`
directly can be more efficient than creating [`PackedSequence`] instances:

```rust
use bitnuc::{as_2bit, from_2bit};
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Efficient k-mer counting
    let mut kmer_counts = HashMap::new();

    // Pack k-mers directly into u64s
    let sequence = b"ACGTACGT";
    for window in sequence.windows(4) {
        let packed = as_2bit(window)?;
        *kmer_counts.entry(packed).or_insert(0) += 1;
    }

    // Count of "ACGT"
    let acgt_packed = as_2bit(b"ACGT")?;
    assert_eq!(kmer_counts.get(&acgt_packed), Some(&2));
    Ok(())
}
```

See the documentation for [`as_2bit`] and [`from_2bit`] for more details on
working with packed sequences directly.
