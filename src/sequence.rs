use crate::error::NucleotideError;
use crate::utils::packing::as_2bit;
use std::ops::Range;

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct PackedSequence {
    data: Vec<u64>,
    length: usize,
}

impl PackedSequence {
    pub fn new(seq: &[u8]) -> Result<Self, NucleotideError> {
        let chunks = (seq.len() + 31) / 32;
        let mut data = Vec::with_capacity(chunks);

        for chunk in seq.chunks(32) {
            data.push(as_2bit(chunk)?);
        }

        Ok(Self {
            data,
            length: seq.len(),
        })
    }

    pub fn len(&self) -> usize {
        self.length
    }

    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    pub fn get(&self, index: usize) -> Result<u8, NucleotideError> {
        if index >= self.length {
            return Err(NucleotideError::IndexOutOfBounds {
                index,
                length: self.length,
            });
        }

        let chunk_idx = index / 32;
        let bit_idx = (index % 32) * 2;
        let bits = (self.data[chunk_idx] >> bit_idx) & 0b11;

        Ok(match bits {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => unreachable!(),
        })
    }

    pub fn slice(&self, range: Range<usize>) -> Result<Vec<u8>, NucleotideError> {
        if range.start > range.end || range.end > self.length {
            return Err(NucleotideError::InvalidRange {
                start: range.start,
                end: range.end,
                length: self.length,
            });
        }

        let mut result = Vec::with_capacity(range.end - range.start);
        for i in range {
            result.push(self.get(i)?);
        }
        Ok(result)
    }

    pub fn to_vec(&self) -> Result<Vec<u8>, NucleotideError> {
        self.slice(0..self.length)
    }
}
#[cfg(test)]
mod tests {
    use crate::error::NucleotideError;
    use crate::sequence::PackedSequence;
    use std::collections::HashSet;

    #[test]
    fn test_new_sequence() {
        let seq = PackedSequence::new(b"ACGT").unwrap();
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.to_vec().unwrap(), b"ACGT");
    }

    #[test]
    fn test_sequence_get() {
        let seq = PackedSequence::new(b"ACGT").unwrap();
        assert_eq!(seq.get(0).unwrap(), b'A');
        assert_eq!(seq.get(1).unwrap(), b'C');
        assert_eq!(seq.get(2).unwrap(), b'G');
        assert_eq!(seq.get(3).unwrap(), b'T');
    }

    #[test]
    fn test_sequence_get_out_of_bounds() {
        let seq = PackedSequence::new(b"ACGT").unwrap();
        assert!(matches!(
            seq.get(4),
            Err(NucleotideError::IndexOutOfBounds {
                index: 4,
                length: 4
            })
        ));
    }

    #[test]
    fn test_sequence_slice() {
        let seq = PackedSequence::new(b"ACGTACGT").unwrap();
        assert_eq!(seq.slice(1..5).unwrap(), b"CGTA");
    }

    #[test]
    fn test_sequence_invalid_slice() {
        let seq = PackedSequence::new(b"ACGT").unwrap();
        assert!(matches!(
            seq.slice(3..2),
            Err(NucleotideError::InvalidRange {
                start: 3,
                end: 2,
                length: 4
            })
        ));
    }

    #[test]
    fn test_sequence_equality() {
        let seq1 = PackedSequence::new(b"ACGT").unwrap();
        let seq2 = PackedSequence::new(b"ACGT").unwrap();
        let seq3 = PackedSequence::new(b"TGCA").unwrap();

        assert_eq!(seq1, seq2);
        assert_ne!(seq1, seq3);
    }

    #[test]
    fn test_sequence_hashability() {
        let mut set = HashSet::new();
        let seq1 = PackedSequence::new(b"ACGT").unwrap();
        let seq2 = PackedSequence::new(b"ACGT").unwrap();
        let seq3 = PackedSequence::new(b"TGCA").unwrap();

        set.insert(seq1.clone());
        assert!(set.contains(&seq2));
        assert!(!set.contains(&seq3));
    }
}
