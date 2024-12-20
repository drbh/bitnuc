pub mod analysis;
pub mod packing;
pub mod unpacking;

pub use packing::as_2bit;
pub use unpacking::from_2bit;

#[cfg(test)]
mod testing {
    use super::*;

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
}
