//! # bitnuc
//!
//! A library for efficient nucleotide sequence manipulation using 2-bit encoding.
//!
//! ## Features
//!
//! - 2-bit nucleotide encoding (A=00, C=01, G=10, T=11)
//! - Direct bit manipulation functions for custom implementations
//! - Higher-level sequence type with additional analysis features
//!
//! ## Low-Level Packing Functions
//!
//! For direct bit manipulation, use the `as_2bit` and `from_2bit` functions:
//!
//! ```rust
//! use bitnuc::{as_2bit, from_2bit, from_2bit_alloc};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Pack a sequence into a u64
//!     let packed = as_2bit(b"ACGT")?;
//!     assert_eq!(packed, 0b11100100);
//!
//!     // Unpack back to a sequence using a reusable buffer
//!     let mut unpacked = Vec::new();
//!     from_2bit(packed, 4, &mut unpacked)?;
//!     assert_eq!(&unpacked, b"ACGT");
//!     unpacked.clear();
//!
//!     // Unpack back to a sequence with a reallocation
//!     let unpacked = from_2bit_alloc(packed, 4)?;
//!     assert_eq!(&unpacked, b"ACGT");
//!
//!     Ok(())
//! }
//! ```
//!
//! These functions are useful when you need to:
//! - Implement custom sequence storage
//! - Manipulate sequences at the bit level
//! - Integrate with other bioinformatics tools
//! - Copy sequences more efficiently
//! - Hash sequences more efficiently
//!
//! For example, packing multiple short sequences:
//!
//! ```rust
//! use bitnuc::{as_2bit, from_2bit};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Pack multiple 4-mers into u64s
//!     let kmers = [b"ACGT", b"TGCA", b"GGCC"];
//!     let packed: Vec<u64> = kmers
//!         .into_iter()
//!         .map(|kmer| as_2bit(kmer))
//!         .collect::<Result<_, _>>()?;
//!
//!     // Unpack when needed
//!     let mut unpacked = Vec::new();
//!     from_2bit(packed[0], 4, &mut unpacked)?;
//!     assert_eq!(&unpacked, b"ACGT");
//!     Ok(())
//! }
//! ```
//!
//! ## High-Level Sequence Type
//!
//! For more complex sequence manipulation, use the [`PackedSequence`] type:
//!
//! ```rust
//! use bitnuc::{PackedSequence, GCContent, BaseCount};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let seq = PackedSequence::new(b"ACGTACGT")?;
//!
//!     // Sequence analysis
//!     println!("GC Content: {}%", seq.gc_content());
//!     let [a_count, c_count, g_count, t_count] = seq.base_counts();
//!
//!     // Slicing
//!     let subseq = seq.slice(1..5)?;
//!     assert_eq!(&subseq, b"CGTA");
//!     Ok(())
//! }
//! ```
//!
//! ## Memory Usage
//!
//! The 2-bit encoding provides significant memory savings:
//!
//! ```text
//! Standard encoding: 1 byte per base
//! ACGT = 4 bytes = 32 bits
//!
//! 2-bit encoding: 2 bits per base
//! ACGT = 8 bits
//! ```
//!
//! This means you can store 4 times as many sequences in the same amount of memory.
//!
//! ## Error Handling
//!
//! All operations that could fail return a [`Result`] with [`NucleotideError`]:
//!
//! ```rust
//! use bitnuc::{as_2bit, NucleotideError};
//!
//! // Invalid nucleotide
//! let err = as_2bit(b"ACGN").unwrap_err();
//! assert!(matches!(err, NucleotideError::InvalidBase(b'N')));
//!
//! // Sequence too long
//! let long_seq = vec![b'A'; 33];
//! let err = as_2bit(&long_seq).unwrap_err();
//! assert!(matches!(err, NucleotideError::SequenceTooLong(33)));
//! ```
//!
//! ## Performance Considerations
//!
//! When working with many short sequences (like k-mers), using `as_2bit` and `from_2bit`
//! directly can be more efficient than creating [`PackedSequence`] instances:
//!
//! ```rust
//! use bitnuc::{as_2bit, from_2bit};
//! use std::collections::HashMap;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Efficient k-mer counting
//!     let mut kmer_counts = HashMap::new();
//!
//!     // Pack k-mers directly into u64s
//!     let sequence = b"ACGTACGT";
//!     for window in sequence.windows(4) {
//!         let packed = as_2bit(window)?;
//!         *kmer_counts.entry(packed).or_insert(0) += 1;
//!     }
//!
//!     // Count of "ACGT"
//!     let acgt_packed = as_2bit(b"ACGT")?;
//!     assert_eq!(kmer_counts.get(&acgt_packed), Some(&2));
//!     Ok(())
//! }
//! ```
//!
//! If you are unpacking many sequences, consider reusing a buffer to avoid reallocations:
//!
//! ```rust
//! use bitnuc::{as_2bit, from_2bit};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!
//!     // Pack a sequence
//!     let packed = as_2bit(b"ACGT")?;
//!
//!     // Reusable buffer for unpacking
//!     let mut unpacked = Vec::new();
//!     from_2bit(packed, 4, &mut unpacked)?;
//!     assert_eq!(&unpacked, b"ACGT");
//!     unpacked.clear();
//!
//!     // Pack another sequence
//!     let packed = as_2bit(b"TGCA")?;
//!     from_2bit(packed, 4, &mut unpacked)?;
//!     assert_eq!(&unpacked, b"TGCA");
//!     Ok(())
//! }
//! ```
//!
//!
//! See the documentation for [`as_2bit`] and [`from_2bit`] for more details on
//! working with packed sequences directly.

mod error;
mod sequence;
mod utils;

pub use error::NucleotideError;
pub use sequence::PackedSequence;
pub use utils::{
    analysis::{BaseCount, GCContent},
    as_2bit, decode, encode, from_2bit, from_2bit_alloc,
};

#[cfg(test)]
mod testing {
    use crate::{BaseCount, GCContent, PackedSequence};

    #[test]
    fn test_sequence_creation_and_analysis() {
        let seq = PackedSequence::new(b"ACGTACGT").unwrap();

        // Test basic properties
        assert_eq!(seq.len(), 8);
        assert!(!seq.is_empty());

        // Test sequence content
        assert_eq!(seq.to_vec().unwrap(), b"ACGTACGT");

        // Test analysis
        assert_eq!(seq.gc_content(), 50.0);
        assert_eq!(seq.base_counts(), [2, 2, 2, 2]);
    }

    #[test]
    fn test_sequence_mutations() {
        let seq = PackedSequence::new(b"ACGTACGT").unwrap();

        // Test slicing
        let slice = seq.slice(2..6).unwrap();
        assert_eq!(slice, b"GTAC");

        // Test individual base access
        assert_eq!(seq.get(0).unwrap(), b'A');
        assert_eq!(seq.get(7).unwrap(), b'T');
    }

    #[test]
    fn test_error_handling() {
        // Test invalid sequence
        assert!(PackedSequence::new(b"ACGN").is_err());

        // Test out of bounds
        let seq = PackedSequence::new(b"ACGT").unwrap();
        assert!(seq.get(4).is_err());
        assert!(seq.slice(2..5).is_err());
    }
}
