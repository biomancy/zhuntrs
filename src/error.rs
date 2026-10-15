use std::{error, fmt, io};

/// Errors returned by the prediction engine and prediction streams.
#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    InvalidWindowRange {
        mindn: usize,
        maxdn: usize,
    },
    InvalidThreshold {
        threshold: f64,
    },
    WindowLengthOverflow {
        maxdn: usize,
    },
    InvalidNucleotide {
        value: u8,
    },
    /// The sequence and its required window context exceed addressable storage.
    InputTooLarge,
    ZeroStreamSize,
    StreamExhausted,
    /// A Rust panic interrupted a prediction while its stream state was locked.
    ///
    /// Rust marks the lock as "poisoned" because the panic may have left partial
    /// updates in that state. Later calls return this error instead of resuming
    /// with potentially inconsistent data; create a new stream to continue.
    /// Ordinary returned errors and threads waiting for the lock do not cause it.
    StreamPoisoned,
    Io(io::Error),
}

pub type Result<T> = std::result::Result<T, Error>;

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidWindowRange { mindn, maxdn } => write!(
                f,
                "Window lengths must satisfy 1 <= mindn <= maxdn (got mindn={mindn}, maxdn={maxdn})"
            ),
            Self::InvalidThreshold { threshold } => {
                write!(f, "Threshold must be >= 0 and not NaN (got {threshold})")
            }
            Self::WindowLengthOverflow { maxdn } => {
                write!(
                    f,
                    "Maximum window length is too large: {maxdn} dinucleotides"
                )
            }
            Self::InvalidNucleotide { value } => write!(f, "Invalid nucleotide byte: {value}"),
            Self::InputTooLarge => {
                f.write_str("Sequence and window context exceed the supported length")
            }
            Self::ZeroStreamSize => f.write_str("Stream size cannot be zero"),
            Self::StreamExhausted => f.write_str("No predictions left"),
            Self::StreamPoisoned => f.write_str(
                "Prediction stream cannot continue after an earlier Rust panic; create a new stream"
            ),
            Self::Io(source) => write!(f, "Unable to read sequence: {source}"),
        }
    }
}

impl error::Error for Error {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(source) => Some(source),
            _ => None,
        }
    }
}

impl From<io::Error> for Error {
    fn from(source: io::Error) -> Self {
        Self::Io(source)
    }
}
