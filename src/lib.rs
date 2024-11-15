//! A library for efficient nucleotide sequence manipulation using 2-bit encoding
//!
//! This library provides tools for working with DNA sequences in a memory-efficient
//! manner by storing nucleotides in a 2-bit packed format.

mod error;
mod sequence;
mod utils;

pub use error::NucleotideError;
pub use sequence::PackedSequence;
pub use utils::{
    analysis::{BaseCount, GCContent},
    as_2bit, from_2bit,
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
