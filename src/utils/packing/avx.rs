use crate::NucleotideError;
use super::naive;
use std::arch::x86_64::*;

pub fn as_2bit(seq: &[u8]) -> Result<u64, NucleotideError> {
    if seq.len() > 32 {
        return Err(NucleotideError::SequenceTooLong(seq.len()));
    }

    // For very short sequences, use the simple version
    if seq.len() < 8 {
        return naive::as_2bit(seq);
    }

    // Pre-validate all bases
    if let Some(&invalid) = seq
        .iter()
        .find(|&&b| !matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'))
    {
        return Err(NucleotideError::InvalidBase(invalid));
    }

    let mut packed = 0u64;

    let len = seq.len();
    let simd_len = len - (len % 16); // Process 16 bytes at a time with AVX2
    unsafe {
        // Create constant vectors for comparison and result building
        let zeros = _mm256_set1_epi8(0i8);
        let ones = _mm256_set1_epi8(1i8);
        let twos = _mm256_set1_epi8(2i8);
        let threes = _mm256_set1_epi8(3i8);

        // Load uppercase and lowercase base values
        let a_upper = _mm256_set1_epi8(b'A' as i8);
        let c_upper = _mm256_set1_epi8(b'C' as i8);
        let g_upper = _mm256_set1_epi8(b'G' as i8);
        let t_upper = _mm256_set1_epi8(b'T' as i8);

        for chunk_idx in (0..simd_len).step_by(16) {
            // Load 16 bytes of sequence data
            let chunk = _mm256_loadu_si256(seq[chunk_idx..].as_ptr() as *const __m256i);

            // Convert lowercase to uppercase for simplified comparison
            let uppercased = _mm256_or_si256(
                chunk,
                _mm256_set1_epi8(0x20), // Set bit 5 to convert to lowercase
            );

            // Create masks for each nucleotide
            let c_mask = _mm256_cmpeq_epi8(uppercased, _mm256_set1_epi8(b'c' as i8));
            let g_mask = _mm256_cmpeq_epi8(uppercased, _mm256_set1_epi8(b'g' as i8));
            let t_mask = _mm256_cmpeq_epi8(uppercased, _mm256_set1_epi8(b't' as i8));

            // Build result vector
            let mut result = zeros;
            result = _mm256_blendv_epi8(result, ones, c_mask);
            result = _mm256_blendv_epi8(result, twos, g_mask);
            result = _mm256_blendv_epi8(result, threes, t_mask);

            // Store results
            let mut temp = [0u8; 32];
            _mm256_storeu_si256(temp.as_mut_ptr() as *mut __m256i, result);

            // Pack into u64
            for (i, &val) in temp.iter().take(16).enumerate() {
                packed |= (val as u64) << ((chunk_idx + i) * 2);
            }
        }

        // Handle remaining bases
        seq.iter()
            .skip(simd_len)
            .take(len)
            .enumerate()
            .for_each(|(i, &base)| {
                let bits = match base {
                    b'A' | b'a' => 0b00,
                    b'C' | b'c' => 0b01,
                    b'G' | b'g' => 0b10,
                    b'T' | b't' => 0b11,
                    _ => unreachable!(),
                };
                packed |= (bits as u64) << ((simd_len + i) * 2);
            });
    }

    Ok(packed)
}
