use crate::error::NucleotideError;

pub fn as_2bit(seq: &[u8]) -> Result<u64, NucleotideError> {
    if seq.len() > 32 {
        return Err(NucleotideError::SequenceTooLong(seq.len()));
    }

    let mut packed = 0u64;

    for (i, &base) in seq.iter().enumerate() {
        let bits = match base {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            invalid => return Err(NucleotideError::InvalidBase(invalid)),
        };

        packed |= (bits as u64) << (i * 2);
    }

    Ok(packed)
}

pub fn from_2bit(packed: u64, expected_size: usize) -> Result<Vec<u8>, NucleotideError> {
    if expected_size > 32 {
        return Err(NucleotideError::InvalidLength(expected_size));
    }

    let mut sequence = Vec::with_capacity(expected_size);

    for i in 0..expected_size {
        let bits = (packed >> (i * 2)) & 0b11;
        let base = match bits {
            0b00 => b'A',
            0b01 => b'C',
            0b10 => b'G',
            0b11 => b'T',
            _ => unreachable!(),
        };
        sequence.push(base);
    }

    Ok(sequence)
}

#[cfg(test)]
mod tests {
    use crate::error::NucleotideError;
    use crate::utils::packing::{as_2bit, from_2bit};

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

    #[test]
    fn test_from_2bit_valid_sequence() {
        let tests = vec![
            (0b11100100, 4, b"ACGT"),
            (0b00000000, 4, b"AAAA"),
            (0b11111111, 4, b"TTTT"),
        ];

        for (input, size, expected) in tests {
            assert_eq!(from_2bit(input, size).unwrap(), expected);
        }
    }
}
