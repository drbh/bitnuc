use crate::NucleotideError;
use std::arch::aarch64::*;

/// Unpacks 8 bases at once using NEON SIMD instructions
unsafe fn unpack_8_bases(packed: u64, base_lookup: uint8x8_t) -> uint8x8_t {
    // Extract 2-bit values for 8 bases
    let mut values = [0u8; 8];

    for (idx, v) in values.iter_mut().enumerate() {
        *v = ((packed >> (idx * 2)) & 0b11) as u8;
    }

    // Load these values into a NEON vector
    // First create an 8-byte vector with our indices
    let indices = vld1_u8(values.as_ptr());

    // Use TBL (table lookup) to convert to bases
    // TBL is the NEON equivalent of x86's shuffle
    vtbl1_u8(base_lookup, indices)
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
    // NEONs table lookup works with 8-byte vectors
    let base_lookup = vld1_u8([b'A', b'C', b'G', b'T', 0, 0, 0, 0].as_ptr());

    // Process 8 bases at a time
    for chunk in 0..simd_chunks {
        let result = unpack_8_bases(packed >> (chunk * 16), base_lookup);

        // Store the 8 resulting bases
        let mut temp = [0u8; 8];
        vst1_u8(temp.as_mut_ptr(), result);
        sequence.extend_from_slice(&temp);
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
