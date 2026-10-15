use std::io::Read;
use std::num::NonZeroUsize;

use crate::{Error, Result};

use super::input::Input;
use super::predictor::Predictor;
use super::{Engine, Prediction};

/// Mutable iterator state shared by the native iterator and Python bindings.
pub(crate) struct IteratorState {
    input: Input,
    predictor: Predictor,
    cursor: usize,
}

impl IteratorState {
    pub(crate) fn new<R: Read>(engine: &Engine, sequence: R) -> Result<Self> {
        Ok(Self {
            input: Input::new(sequence, engine.window_size(), engine.wrap())?,
            predictor: Predictor::new(),
            cursor: 0,
        })
    }

    pub(crate) fn next(&mut self, engine: &Engine, chunk: NonZeroUsize) -> Option<Prediction> {
        let windows = self.input.windows(self.cursor, chunk)?;
        let end = windows.starts.end;
        let output = self.predictor.predict(engine, windows);
        self.cursor = end;
        Some(output)
    }
}

/// An iterator owning its parsed sequence and reusing working memory between calls.
pub struct PredictionIterator<'a> {
    engine: &'a Engine,
    state: IteratorState,
    chunk: NonZeroUsize,
}

impl<'a> PredictionIterator<'a> {
    pub(crate) fn new<R: Read>(
        engine: &'a Engine,
        sequence: R,
        chunk: NonZeroUsize,
    ) -> Result<Self> {
        Ok(Self {
            engine,
            state: IteratorState::new(engine, sequence)?,
            chunk,
        })
    }

    /// Return all remaining predictions, or `Error::StreamExhausted`.
    pub fn consume(&mut self) -> Result<Prediction> {
        self.state
            .next(self.engine, NonZeroUsize::MAX)
            .ok_or(Error::StreamExhausted)
    }
}

impl Iterator for PredictionIterator<'_> {
    type Item = Result<Prediction>;

    fn next(&mut self) -> Option<Self::Item> {
        self.state.next(self.engine, self.chunk).map(Ok)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nz(value: usize) -> NonZeroUsize {
        NonZeroUsize::new(value).unwrap()
    }

    fn concatenate(chunks: impl IntoIterator<Item = Prediction>) -> Prediction {
        let mut output = Prediction::default();
        for chunk in chunks {
            output.start.extend(chunk.start);
            output.end.extend(chunk.end);
            output.zhscore.extend(chunk.zhscore);
            output.sequence.extend(chunk.sequence);
            output.conformation.extend(chunk.conformation);
        }
        output
    }

    #[test]
    fn iterators_have_independent_positions() -> Result<()> {
        let sequence = b"ACGTACGT";
        let engine = Engine::new(nz(1), nz(2), 0.0, true)?;
        let expected = engine.predict(sequence.as_slice())?;

        let mut first = engine.stream(sequence.as_slice(), nz(2))?;
        let mut second = engine.stream(sequence.as_slice(), nz(3))?;

        let head = first.next().unwrap()?;
        assert_eq!(second.consume()?, expected);

        let rest = first.consume()?;
        assert_eq!(concatenate([head, rest]), expected);
        assert!(first.next().is_none());
        assert!(matches!(first.consume(), Err(Error::StreamExhausted)));
        assert!(matches!(second.consume(), Err(Error::StreamExhausted)));
        Ok(())
    }
}
