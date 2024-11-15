use std::fmt;

#[derive(Debug, PartialEq, Eq)]
pub enum NucleotideError {
    InvalidBase(u8),
    SequenceTooLong(usize),
    InvalidLength(usize),
    IndexOutOfBounds {
        index: usize,
        length: usize,
    },
    InvalidRange {
        start: usize,
        end: usize,
        length: usize,
    },
}

impl fmt::Display for NucleotideError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NucleotideError::InvalidBase(b) => write!(f, "Invalid nucleotide base: {}", b),
            NucleotideError::SequenceTooLong(len) => {
                write!(f, "Sequence length {} exceeds maximum", len)
            }
            NucleotideError::InvalidLength(len) => write!(f, "Invalid length: {}", len),
            NucleotideError::IndexOutOfBounds { index, length } => {
                write!(
                    f,
                    "Index {} out of bounds for sequence of length {}",
                    index, length
                )
            }
            NucleotideError::InvalidRange { start, end, length } => {
                write!(
                    f,
                    "Invalid range {}..{} for sequence of length {}",
                    start, end, length
                )
            }
        }
    }
}

impl std::error::Error for NucleotideError {}
