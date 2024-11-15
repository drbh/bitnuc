use crate::error::NucleotideError;

/// Converts a nucleotide sequence into a 2-bit packed representation.
///
/// Each nucleotide is encoded using 2 bits:
/// - A/a = 00
/// - C/c = 01
/// - G/g = 10
/// - T/t = 11
///
/// The bases are packed from least significant to most significant bits.
/// For example, "ACGT" becomes 0b11100100.
///
/// # Arguments
///
/// * `seq` - A byte slice containing ASCII nucleotides (A,C,G,T, case insensitive)
///
/// # Returns
///
/// Returns a `u64` containing the packed representation.
///
/// # Errors
///
/// Returns `NucleotideError::InvalidBase` if the sequence contains any characters
/// other than A,C,G,T (case insensitive).
///
/// Returns `NucleotideError::SequenceTooLong` if the input sequence is longer
/// than 32 bases (as a u64 can only store 32 * 2 bits).
///
/// # Examples
///
/// Basic packing:
/// ```rust
/// use bitnuc::as_2bit;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let packed = as_2bit(b"ACGT")?;
/// assert_eq!(packed, 0b11100100);
/// # Ok(())
/// # }
/// ```
///
/// Case insensitivity:
/// ```rust
/// use bitnuc::as_2bit;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// assert_eq!(as_2bit(b"acgt")?, as_2bit(b"ACGT")?);
/// # Ok(())
/// # }
/// ```
///
/// Error handling:
/// ```rust
/// use bitnuc::{as_2bit, NucleotideError};
///
/// # fn main() {
/// // Invalid base
/// assert!(matches!(
///     as_2bit(b"ACGN"),
///     Err(NucleotideError::InvalidBase(b'N'))
/// ));
///
/// // Sequence too long
/// let long_seq = vec![b'A'; 33];
/// assert!(matches!(
///     as_2bit(&long_seq),
///     Err(NucleotideError::SequenceTooLong(33))
/// ));
/// # }
/// ```
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

/// Converts a 2-bit packed representation back into a nucleotide sequence.
///
/// This function reverses the packing performed by `as_2bit`.
///
/// # Arguments
///
/// * `packed` - A u64 containing the 2-bit packed sequence
/// * `expected_size` - The number of bases to unpack
///
/// # Returns
///
/// Returns a `Vec<u8>` containing the ASCII sequence.
///
/// # Errors
///
/// Returns `NucleotideError::InvalidLength` if `expected_size` is greater than 32
/// (as a u64 can only store 32 * 2 bits).
///
/// # Examples
///
/// Basic unpacking:
/// ```rust
/// use bitnuc::{as_2bit, from_2bit};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// // Pack and unpack
/// let packed = as_2bit(b"ACGT")?;
/// let unpacked = from_2bit(packed, 4)?;
/// assert_eq!(&unpacked, b"ACGT");
///
/// // Partial unpacking
/// let partial = from_2bit(packed, 2)?;
/// assert_eq!(&partial, b"AC");
/// # Ok(())
/// # }
/// ```
///
/// Error handling:
/// ```rust
/// use bitnuc::{from_2bit, NucleotideError};
///
/// # fn main() {
/// // Length too long
/// assert!(matches!(
///     from_2bit(0, 33),
///     Err(NucleotideError::InvalidLength(33))
/// ));
/// # }
/// ```
///
/// # Implementation Details
///
/// The bases are unpacked from least significant to most significant bits:
/// ```rust
/// use bitnuc::{as_2bit, from_2bit};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let packed = 0b11100100; // "ACGT" in 2-bit encoding
/// let seq = from_2bit(packed, 4)?;
/// assert_eq!(seq[0], b'A'); // 00
/// assert_eq!(seq[1], b'C'); // 01
/// assert_eq!(seq[2], b'G'); // 10
/// assert_eq!(seq[3], b'T'); // 11
/// # Ok(())
/// # }
/// ```
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

    #[test]
    fn test_pack_unpack_roundtrip() {
        let test_cases = [
            b"A".as_slice(),
            b"C".as_slice(),
            b"G".as_slice(),
            b"T".as_slice(),
            b"AC".as_slice(),
            b"GT".as_slice(),
            b"ACG".as_slice(),
            b"TGC".as_slice(),
            b"ACGT".as_slice(),
            b"TGCA".as_slice(),
            b"ACGTACGT".as_slice(),
            b"AAAA".as_slice(),
            b"CCCC".as_slice(),
            b"GGGG".as_slice(),
            b"TTTT".as_slice(),
        ];

        for &seq in &test_cases {
            let packed = as_2bit(seq).unwrap();
            let unpacked = from_2bit(packed, seq.len()).unwrap();
            assert_eq!(seq, unpacked.as_slice());
        }
    }

    #[test]
    fn test_partial_unpack() {
        let packed = as_2bit(b"ACGT").unwrap();
        assert_eq!(from_2bit(packed, 2).unwrap(), b"AC");
        assert_eq!(from_2bit(packed, 3).unwrap(), b"ACG");
    }
}
