use super::naive;
use crate::error::NucleotideError;
use std::arch::aarch64::*;

/// Represents the 2-bit encoding for each nucleotide
#[repr(u8)]
enum NucleotideBits {
    A = 0b00,
    C = 0b01,
    G = 0b10,
    T = 0b11,
}

/// A reusable structure holding common SIMD constants
#[repr(align(16))] // Ensure proper alignment for SIMD
struct SimdConstants {
    zeros: uint8x8_t,
    ones: uint8x8_t,
    twos: uint8x8_t,
    threes: uint8x8_t,
}

impl SimdConstants {
    #[inline(always)]
    unsafe fn new() -> Self {
        Self {
            zeros: vdup_n_u8(NucleotideBits::A as u8),
            ones: vdup_n_u8(NucleotideBits::C as u8),
            twos: vdup_n_u8(NucleotideBits::G as u8),
            threes: vdup_n_u8(NucleotideBits::T as u8),
        }
    }
}

/// Creates a bitmask for matching both upper and lowercase versions of a nucleotide
///
/// This function combines equality comparisons for both cases of a nucleotide
/// to create a single mask where matching positions are set to all 1s.
#[inline(always)]
unsafe fn create_dual_pattern_mask(chunk: uint8x8_t, upper: u8, lower: u8) -> uint8x8_t {
    vorr_u8(
        vceq_u8(chunk, vdup_n_u8(upper)),
        vceq_u8(chunk, vdup_n_u8(lower)),
    )
}

/// Sets the appropriate 2-bit patterns based on nucleotide masks
#[inline(always)]
unsafe fn set_bits(
    c_mask: uint8x8_t,
    g_mask: uint8x8_t,
    t_mask: uint8x8_t,
    constants: &SimdConstants,
) -> uint8x8_t {
    let mut result = constants.zeros;
    // Use BSL (bit select) to set appropriate values based on masks
    result = vbsl_u8(c_mask, constants.ones, result);
    result = vbsl_u8(g_mask, constants.twos, result);
    result = vbsl_u8(t_mask, constants.threes, result);
    result
}

/// Processes a single SIMD chunk of 8 nucleotides
#[inline(always)]
unsafe fn process_simd_chunk(chunk: uint8x8_t, constants: &SimdConstants) -> uint8x8_t {
    let (c_mask, g_mask, t_mask) = (
        create_dual_pattern_mask(chunk, b'C', b'c'),
        create_dual_pattern_mask(chunk, b'G', b'g'),
        create_dual_pattern_mask(chunk, b'T', b't'),
    );
    set_bits(c_mask, g_mask, t_mask, constants)
}

#[cfg(target_arch = "aarch64")]
pub fn as_2bit(seq: &[u8]) -> Result<u64, NucleotideError> {
    if seq.len() > 32 {
        return Err(NucleotideError::SequenceTooLong(seq.len()));
    }

    if seq.len() < 8 {
        return naive::as_2bit(seq);
    }

    // Pre-validate all bases using SIMD when possible
    if let Some(&invalid) = seq
        .iter()
        .find(|&&b| !matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'))
    {
        return Err(NucleotideError::InvalidBase(invalid));
    }

    let mut packed = 0u64;
    let len = seq.len();
    let simd_len = len - (len % 8);

    unsafe {
        let constants = SimdConstants::new();

        // Process 8 nucleotides at a time using SIMD
        for chunk_idx in (0..simd_len).step_by(8) {
            let chunk = vld1_u8(seq[chunk_idx..].as_ptr());
            let result = process_simd_chunk(chunk, &constants);

            // Store and pack results
            let mut temp = [0u8; 8];
            vst1_u8(temp.as_mut_ptr(), result);

            for (i, &val) in temp.iter().enumerate() {
                packed |= (val as u64) << ((chunk_idx + i) * 2);
            }
        }

        // Handle remaining nucleotides
        for (i, &base) in seq.iter().skip(simd_len).enumerate() {
            let bits = match base {
                b'A' | b'a' => NucleotideBits::A as u64,
                b'C' | b'c' => NucleotideBits::C as u64,
                b'G' | b'g' => NucleotideBits::G as u64,
                b'T' | b't' => NucleotideBits::T as u64,
                _ => unreachable!(),
            };
            packed |= bits << ((simd_len + i) * 2);
        }
    }

    Ok(packed)
}
