use crate::NucleotideError;

#[cfg(target_arch = "aarch64")]
mod aarch64;
#[cfg(target_arch = "x86_64")]
mod avx;
mod naive;

pub fn as_2bit(seq: &[u8]) -> Result<u64, NucleotideError> {
    #[cfg(target_arch = "aarch64")]
    return aarch64::as_2bit(seq);

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx")]
    return avx::as_2bit(seq);

    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    return naive::as_2bit(seq);
}

#[cfg(test)]
mod testing {
    use super::*;

    #[test]
    fn test_as_2bit_valid_sequence() {
        let tests = vec![
            (b"ACGT", 0b11100100),
            (b"AAAA", 0b00000000),
            (b"TTTT", 0b11111111),
            (b"GGGG", 0b10101010),
            (b"CCCC", 0b01010101),
        ];

        for (input, expected) in tests {
            assert_eq!(as_2bit(input).unwrap(), expected);
        }
    }

    #[test]
    fn test_as_2bit_alignments() {
        let tests = vec![(b"ACTGGAAAATTTTAAGG", 0b1010000011111111000000001010110100)];
        for (input, expected) in tests {
            assert_eq!(as_2bit(input).unwrap(), expected);
        }
    }

    #[test]
    fn test_as_2bit_lowercase() {
        assert_eq!(as_2bit(b"acgt").unwrap(), as_2bit(b"ACGT").unwrap());
    }

    #[test]
    fn test_as_2bit_invalid_base() {
        let result = as_2bit(b"ACGN");
        assert!(matches!(result, Err(NucleotideError::InvalidBase(b'N'))));
    }

    #[test]
    fn test_as_2bit_sequence_too_long() {
        let long_seq = vec![b'A'; 33];
        assert!(matches!(
            as_2bit(&long_seq),
            Err(NucleotideError::SequenceTooLong(33))
        ));
    }
}
