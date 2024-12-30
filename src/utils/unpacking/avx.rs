use crate::NucleotideError;
use std::arch::x86_64::*;

unsafe fn unpack_8_bases(packed: u64, base_lookup: __m256i) -> __m256i {
    // Extract bits for 8 bases all at once
    let mut values = [0u8; 8];
    for i in 0..8 {
        values[i] = ((packed >> (i * 2)) & 0b11) as u8;
    }

    // Load these values into a SIMD vector
    let indices = _mm256_set1_epi64x(
        ((values[7] as i64) << 56)
            | ((values[6] as i64) << 48)
            | ((values[5] as i64) << 40)
            | ((values[4] as i64) << 32)
            | ((values[3] as i64) << 24)
            | ((values[2] as i64) << 16)
            | ((values[1] as i64) << 8)
            | (values[0] as i64),
    );

    // Use shuffle to convert to bases
    _mm256_shuffle_epi8(base_lookup, indices)
}

pub unsafe fn from_2bit_simd(
    packed: u64,
    expected_size: usize,
    sequence: &mut Vec<u8>,
) -> Result<(), NucleotideError> {
    if expected_size > 32 {
        return Err(NucleotideError::InvalidLength(expected_size));
    }

    // Only handle full 8-base chunks with SIMD
    let simd_chunks = expected_size / 8;
    sequence.reserve(expected_size);

    // Create lookup table for each possible 2-bit value
    let base_lookup = _mm256_setr_epi8(
        b'A' as i8, b'C' as i8, b'G' as i8, b'T' as i8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        b'A' as i8, b'C' as i8, b'G' as i8, b'T' as i8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    );

    // Process 8 bases at a time
    for chunk in 0..simd_chunks {
        let result = unpack_8_bases(packed >> (chunk * 16), base_lookup);

        let mut temp = [0u8; 32];
        _mm256_storeu_si256(temp.as_mut_ptr() as *mut __m256i, result);
        sequence.extend_from_slice(&temp[..8]);
    }

    // Handle remaining bases
    let remaining_start = simd_chunks * 8;
    for i in remaining_start..expected_size {
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

    Ok(())
}

#[cfg(test)]
mod testing {
    use super::*;
    use crate::as_2bit;

    #[test]
    fn test_from_2bit_simd() {
        let expected = b"ACTGACTGACTGACTGACTG";
        let packed = as_2bit(expected).unwrap();
        let mut observed = Vec::new();
        unsafe {
            from_2bit_simd(packed, 20, &mut observed).unwrap();
        }
        assert_eq!(&observed, expected);
    }

    #[test]
    fn test_various_lengths() {
        for len in 1..=32 {
            let input = b"ACTGACTGACTGACTGACTGACTGACTGACTG";
            let packed = as_2bit(&input[..len]).unwrap();
            let mut observed = Vec::new();
            unsafe {
                from_2bit_simd(packed, len, &mut observed).unwrap();
            }
            assert_eq!(&observed, &input[..len], "Failed at length {}", len);
        }
    }

    #[test]
    fn test_append() {
        let sequence = b"ACTGACTGACTGACTGACTG";
        let packed = as_2bit(sequence).unwrap();
        let mut observed = Vec::new();
        unsafe {
            from_2bit_simd(packed, 10, &mut observed).unwrap();
            from_2bit_simd(packed, 10, &mut observed).unwrap();
        }
        let expected = b"ACTGACTGACACTGACTGAC"; // repeated out of phase
        assert_eq!(&observed, expected);
    }
}
