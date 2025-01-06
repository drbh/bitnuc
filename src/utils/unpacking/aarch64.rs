use crate::NucleotideError;
use std::arch::aarch64::*;

#[inline(always)]
unsafe fn unpack_8_bases(packed: u64, lookup: uint8x8_t) -> uint8x8_t {
    // Create indices array for 8 bases
    let mut indices = [0u8; 8];
    for (i, v) in indices.iter_mut().enumerate() {
        *v = ((packed >> (i * 2)) & 0b11) as u8;
    }

    // Load indices and perform table lookup
    let index_vec = vld1_u8(indices.as_ptr());
    vtbl1_u8(lookup, index_vec)
}

#[inline(always)]
unsafe fn process_remainder(packed: u64, start: usize, end: usize, sequence: &mut Vec<u8>) {
    static LOOKUP: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let count = end - start;
    let old_len = sequence.len();
    sequence.reserve(count);

    let ptr = sequence.as_mut_ptr().add(old_len);
    for i in 0..count {
        let bits = (packed >> ((start + i) * 2)) & 0b11;
        *ptr.add(i) = LOOKUP[bits as usize];
    }
    sequence.set_len(old_len + count);
}

pub unsafe fn from_2bit_simd(
    packed: u64,
    expected_size: usize,
    sequence: &mut Vec<u8>,
) -> Result<(), NucleotideError> {
    if expected_size > 32 {
        return Err(NucleotideError::InvalidLength(expected_size));
    }

    sequence.reserve(expected_size);

    // Create lookup table for base conversion
    let lookup = vld1_u8([b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'].as_ptr());

    // Process in chunks of 8 bases
    let simd_chunks = expected_size / 8;
    for chunk in 0..simd_chunks {
        let chunk_data = packed >> (chunk * 16);
        let result = unpack_8_bases(chunk_data, lookup);
        let mut temp = [0u8; 8];
        vst1_u8(temp.as_mut_ptr(), result);
        sequence.extend_from_slice(&temp);
    }

    // Handle remaining bases
    let remaining_start = simd_chunks * 8;
    if remaining_start < expected_size {
        process_remainder(packed, remaining_start, expected_size, sequence);
    }

    Ok(())
}

pub unsafe fn from_2bit_multi_simd(
    ebuf: &[u64],
    n_bases: usize,
    sequence: &mut Vec<u8>,
) -> Result<(), NucleotideError> {
    sequence.reserve(n_bases);

    // Create lookup table
    let lookup = vld1_u8([b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'].as_ptr());

    // Process full u64 chunks (32 bases each)
    let full_chunks = n_bases / 32;
    for chunk in ebuf.iter().take(full_chunks) {
        // Process each 32-base chunk in 8-base segments
        for segment in 0..4 {
            let segment_data = chunk >> (segment * 16);
            let result = unpack_8_bases(segment_data, lookup);
            let mut temp = [0u8; 8];
            vst1_u8(temp.as_mut_ptr(), result);
            sequence.extend_from_slice(&temp);
        }
    }

    // Handle remaining bases in the last chunk
    let remaining_bases = n_bases % 32;
    if remaining_bases > 0 {
        let last_chunk = ebuf[full_chunks];

        // Process full 8-base segments in the last chunk
        let remaining_segments = remaining_bases / 8;
        for segment in 0..remaining_segments {
            let segment_data = last_chunk >> (segment * 16);
            let result = unpack_8_bases(segment_data, lookup);
            let mut temp = [0u8; 8];
            vst1_u8(temp.as_mut_ptr(), result);
            sequence.extend_from_slice(&temp);
        }

        // Handle any remaining bases
        let final_remaining = remaining_bases % 8;
        if final_remaining > 0 {
            let start = remaining_segments * 8;
            process_remainder(last_chunk, start, remaining_bases, sequence);
        }
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
        let expected = b"ACTGACTGACACTGACTGAC"; // Two copies of the first 10 bases
        assert_eq!(&observed, expected);
    }
}
