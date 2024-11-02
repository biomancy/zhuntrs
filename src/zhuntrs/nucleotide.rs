use std::fmt::Display;

use eyre::{eyre, Error, Result};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[repr(u8)]
pub enum Nucleotide {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4,
}

impl TryFrom<u8> for Nucleotide {
    type Error = Error;

    fn try_from(value: u8) -> Result<Self> {
        match value {
            b'A' | b'a' => Ok(Nucleotide::A),
            b'C' | b'c' => Ok(Nucleotide::C),
            b'G' | b'g' => Ok(Nucleotide::G),
            b'T' | b't' => Ok(Nucleotide::T),
            b'N' | b'n' => Ok(Nucleotide::N),
            _ => Err(eyre!("Invalid nucleotide: {value}")),
        }
    }
}

impl Display for Nucleotide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            Nucleotide::A => "A".to_string(),
            Nucleotide::C => "C".to_string(),
            Nucleotide::G => "G".to_string(),
            Nucleotide::T => "T".to_string(),
            Nucleotide::N => "N".to_string(),
        };
        write!(f, "{}", str)
    }
}
