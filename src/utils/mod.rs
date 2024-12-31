pub mod analysis;
pub mod packing;
pub mod unpacking;

pub use packing::as_2bit;
pub use unpacking::{from_2bit, from_2bit_alloc};

use crate::NucleotideError;

/// Encode a sequence into a buffer of 2-bit encoded nucleotides.
///
/// # Arguments
///
/// * `sequence` - The nucleotide sequence to encode.
/// * `ebuf` - The buffer to write the encoded nucleotides to.
///
/// # Errors
///
/// If the sequence cannot be encoded, an error is returned.
pub fn encode(sequence: &[u8], ebuf: &mut Vec<u64>) -> Result<(), NucleotideError> {
    // Clear the buffer
    ebuf.clear();

    // Calculate the number of chunks
    let n_chunks = sequence.len().div_ceil(32);

    let mut l_bounds = 0;
    for _ in 0..n_chunks - 1 {
        let r_bounds = l_bounds + 32;
        let chunk = &sequence[l_bounds..r_bounds];

        let bits = as_2bit(chunk)?;
        ebuf.push(bits);
        l_bounds = r_bounds;
    }

    let bits = as_2bit(&sequence[l_bounds..])?;
    ebuf.push(bits);

    Ok(())
}

/// Unpacks a 2-bit packed sequence into a nucleotide sequence.
///
/// The sequence is a collection of u64 values, where each u64 contains up to 32 nucleotides.
///
/// It is expected that the sequence is packed fully (32bp) until the final u64, which may contain
/// fewer than 32 nucleotides.
///
/// # Arguments
///
/// * `ebuf` - The buffer containing the packed nucleotides.
/// * `n_bases` - The number of nucleotides to unpack.
/// * `dbuf` - The buffer to write the unpacked nucleotides to.
///
/// # Errors
///
/// If the sequence cannot be unpacked, an error is returned.
pub fn decode(ebuf: &[u64], n_bases: usize, dbuf: &mut Vec<u8>) -> Result<(), NucleotideError> {
    // Calculate the number of chunks and the remainder
    let n_chunks = n_bases.div_ceil(32);
    let rem = match n_bases % 32 {
        0 => 32, // Full chunk
        rem => rem,
    };

    // Process all chunks except the last one
    ebuf.iter()
        .take(n_chunks - 1)
        .try_for_each(|component| from_2bit(*component, 32, dbuf))?;

    // Process the last one with the remainder
    ebuf.get(n_chunks - 1)
        .map_or(Err(NucleotideError::InvalidLength(n_bases)), |&component| {
            from_2bit(component, rem, dbuf)
        })?;

    Ok(())
}

#[cfg(test)]
mod testing {
    use super::*;
    use crate::NucleotideError;
    use nucgen::Sequence;

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
        let mut unpacked = Vec::new();

        for &seq in &test_cases {
            let packed = as_2bit(seq).unwrap();
            from_2bit(packed, seq.len(), &mut unpacked).unwrap();
            assert_eq!(seq, unpacked.as_slice());
            unpacked.clear();
        }
    }

    #[test]
    fn test_partial_unpack() {
        let packed = as_2bit(b"ACGT").unwrap();
        let mut unpacked = Vec::new();

        from_2bit(packed, 2, &mut unpacked).unwrap();
        assert_eq!(unpacked, b"AC");
        unpacked.clear();

        from_2bit(packed, 3, &mut unpacked).unwrap();
        assert_eq!(unpacked, b"ACG");
        unpacked.clear();
    }

    #[test]
    fn test_large_sequence_round_trip() -> Result<(), NucleotideError> {
        let mut rng = rand::thread_rng();
        let mut seq = Sequence::new();

        let sizes = 1..=1000;

        for len in sizes {
            seq.fill_buffer(&mut rng, len);

            let mut ebuf = Vec::new();
            encode(seq.bytes(), &mut ebuf)?;

            let mut unpacked = Vec::new();
            decode(&ebuf, len, &mut unpacked)?;

            assert_eq!(seq.bytes(), unpacked.as_slice());
        }

        Ok(())
    }
}
