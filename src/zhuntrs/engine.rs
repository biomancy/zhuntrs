use std::io::Read;
use std::num::NonZeroUsize;

use crate::{Error, Result};

use super::dbzed::DiNucleotide;
use super::input::Input;
use super::iterator::PredictionIterator;
use super::nucleotide::Nucleotide;
use super::predictor::Predictor;

/// Columns describing locally optimal Z-DNA windows.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Prediction {
    pub start: Vec<usize>,
    pub end: Vec<usize>,
    pub zhscore: Vec<f64>,
    pub sequence: Vec<Vec<Nucleotide>>,
    pub conformation: Vec<Vec<DiNucleotide>>,
}

/// Validated, immutable Z-Hunt parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Engine {
    mindn: NonZeroUsize,
    maxdn: NonZeroUsize,
    threshold: f64,
    window_size: NonZeroUsize,
    wrap: bool,
}

impl Engine {
    /// Validate and construct a new engine.
    ///
    /// Window lengths are in dinucleotides. `wrap` enables circular sequences;
    /// only windows whose score is at least `threshold` are returned.
    ///
    /// Window lengths must satisfy `1 <= mindn <= maxdn`, and the threshold
    /// must be nonnegative and not NaN. The nucleotide width `2 * maxdn`
    /// must also fit in `usize`.
    pub fn new(
        mindn: NonZeroUsize,
        maxdn: NonZeroUsize,
        threshold: f64,
        wrap: bool,
    ) -> Result<Self> {
        if mindn > maxdn {
            return Err(Error::InvalidWindowRange {
                mindn: mindn.get(),
                maxdn: maxdn.get(),
            });
        }
        if threshold.is_nan() || threshold < 0.0 {
            return Err(Error::InvalidThreshold { threshold });
        }
        if maxdn.get().checked_mul(2).is_none() {
            return Err(Error::WindowLengthOverflow { maxdn: maxdn.get() });
        }
        let window_size = NonZeroUsize::new(maxdn.get() * 2).unwrap();

        Ok(Self {
            mindn,
            maxdn,
            threshold,
            window_size,
            wrap,
        })
    }

    pub fn mindn(&self) -> NonZeroUsize {
        self.mindn
    }

    pub fn maxdn(&self) -> NonZeroUsize {
        self.maxdn
    }

    pub fn window_size(&self) -> NonZeroUsize {
        self.window_size
    }

    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    pub fn wrap(&self) -> bool {
        self.wrap
    }

    /// Predict all eligible windows in a sequence of ASCII A, C, G, T, or N.
    ///
    /// Input is case insensitive. Empty sequences produce empty predictions.
    pub fn predict<R: Read>(&self, sequence: R) -> Result<Prediction> {
        let input = Input::new(sequence, self.window_size(), self.wrap)?;
        match input.windows(0, NonZeroUsize::MAX) {
            Some(windows) => Ok(Predictor::new().predict(self, windows)),
            None => Ok(Prediction::default()),
        }
    }

    /// Stream predictions for up to `chunk` consecutive start positions per iteration.
    ///
    /// The source is read, validated, and extended for circular windows before
    /// returning. The iterator owns the prepared sequence and borrows this
    /// engine. Filtering may produce empty chunks; `chunk` must be nonzero.
    pub fn stream<R: Read>(
        &self,
        sequence: R,
        chunk: NonZeroUsize,
    ) -> Result<PredictionIterator<'_>> {
        PredictionIterator::new(self, sequence, chunk)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nz(value: usize) -> NonZeroUsize {
        NonZeroUsize::new(value).unwrap()
    }

    #[test]
    fn new_rejects_invalid_parameters() {
        for (mindn, maxdn, threshold) in [
            (3, 2, 0.0),
            (1, 2, -1.0),
            (1, 2, f64::NAN),
            (1, usize::MAX, 0.0),
        ] {
            let result = Engine::new(nz(mindn), nz(maxdn), threshold, false);
            match (mindn, maxdn) {
                (3, 2) => assert!(matches!(
                    result,
                    Err(Error::InvalidWindowRange { mindn: 3, maxdn: 2 })
                )),
                (_, usize::MAX) => assert!(matches!(
                    result,
                    Err(Error::WindowLengthOverflow { maxdn: usize::MAX })
                )),
                _ => assert!(matches!(result, Err(Error::InvalidThreshold { .. }))),
            }
        }
    }

    #[test]
    fn supports_concurrent_predictions() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<Engine>();
    }
}
