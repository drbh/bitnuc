use crate::sequence::PackedSequence;

pub trait GCContent {
    fn gc_content(&self) -> f64;
}

impl GCContent for PackedSequence {
    fn gc_content(&self) -> f64 {
        let seq = self.to_vec().unwrap_or_default();
        if seq.is_empty() {
            0.0
        } else {
            let gc_count = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
            (gc_count as f64 / self.len() as f64) * 100.0
        }
    }
}

pub trait BaseCount {
    fn base_counts(&self) -> [usize; 4];
}

impl BaseCount for PackedSequence {
    fn base_counts(&self) -> [usize; 4] {
        let seq = self.to_vec().unwrap_or_default();
        let mut counts = [0; 4];
        for &base in &seq {
            let idx = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => continue,
            };
            counts[idx] += 1;
        }
        counts
    }
}

#[cfg(test)]
mod tests {
    use crate::sequence::PackedSequence;
    use crate::utils::analysis::{BaseCount, GCContent};

    #[test]
    fn test_gc_content() {
        let tests = vec![
            (b"ACGT".as_slice(), 50.0),
            (b"AAAA".as_slice(), 0.0),
            (b"CCCC".as_slice(), 100.0),
            (b"AACG".as_slice(), 50.0),
            (b"ACGTA".as_slice(), 40.0),
        ];

        for (seq, expected) in tests {
            let packed = PackedSequence::new(seq).unwrap();
            assert_eq!(packed.gc_content(), expected);
        }
    }

    #[test]
    fn test_base_counts() {
        let tests = vec![
            (b"ACGT".as_slice(), [1, 1, 1, 1]),
            (b"AAAA".as_slice(), [4, 0, 0, 0]),
            (b"CCCC".as_slice(), [0, 4, 0, 0]),
            (b"AACG".as_slice(), [2, 1, 1, 0]),
            (b"ACGTA".as_slice(), [2, 1, 1, 1]),
        ];

        for (seq, expected) in tests {
            let packed = PackedSequence::new(seq).unwrap();
            assert_eq!(packed.base_counts(), expected);
        }
    }

    #[test]
    fn test_empty_sequence_analysis() {
        let seq = PackedSequence::new(b"").unwrap();
        assert_eq!(seq.gc_content(), 0.0);
        assert_eq!(seq.base_counts(), [0, 0, 0, 0]);
    }
}
