use eyre::{eyre, WrapErr};
use pyo3::prelude::*;
use pyo3::types::PyBytes;

use crate::zhuntrs;

type ZHOutput = (Vec<usize>, Vec<usize>, Vec<f64>, Vec<String>, Vec<String>);

/// Estimate Z-DNA propensity of a given DNA sequence using the Z-Hunt algorithm.
///
/// # Arguments
/// * `sequence` - The sequence to be analyzed. Only A, C, G, T, and N are allowed. Case insensitive.
/// * `mindn` - Minimum length of the Z-DNA window, in dinucleotides.
/// * `maxdn` - Maximum length of the Z-DNA window, in dinucleotides.
/// * `threshold` - Only DNA windows with a ZH-score above this threshold will be returned.
/// * `wrap` - If true, the sequence will be wrapped around to estimate the ZH-scores of DNA windows ending outside the sequence.
///
/// # Returns
/// A tuple with 5 vectors describing locally ZH-optimal DNA windows:
/// * `start` - Start coordinate of the i-th DNA window, 0-based.
/// * `end` - End coordinate of the i-th DNA window, 0-based.
/// * `score` - Locally optimal ZH-score for the i-th DNA window.
/// * `sequence` - Sequence of the i-th DNA window.
/// * `conformation` - Conformation of the i-th DNA window.
///
/// Note: A locally optimal window is defined as a window with the maximum ZH-score among DNA windows starting at the i-th position.
#[pyfunction(signature = (sequence, mindn, maxdn, threshold, wrap = false))]
pub fn predict(
    sequence: &[u8],
    mindn: u8,
    maxdn: u8,
    threshold: f64,
    wrap: bool,
) -> PyResult<ZHOutput> {
    let wrap = match wrap {
        true => zhuntrs::TailExtension::WrapAround(maxdn as usize * 2),
        false => zhuntrs::TailExtension::None,
    };

    let mut engine = zhuntrs::Engine::new(None);
    let mut buffer = zhuntrs::Prediction::default();
    engine
        .predict(
            sequence,
            mindn as usize,
            maxdn as usize,
            threshold,
            wrap,
            &mut buffer,
        )
        .wrap_err("Z-Hunt prediction failed")?;

    let sequence = buffer
        .sequence
        .into_iter()
        .map(|x| x.into_iter().map(|n| n.to_string()).collect())
        .collect();

    let conformation = buffer
        .conformation
        .into_iter()
        .map(|x| x.into_iter().map(|dn| dn.to_string()).collect())
        .collect();

    Ok((
        buffer.start,
        buffer.end,
        buffer.zhscore,
        sequence,
        conformation,
    ))
}

/// Estimate Z-DNA propensity of a given DNA sequence using the Z-Hunt algorithm. Predictions are streamed in chunks of a given size.
///
/// # Arguments
/// * `sequence` - The sequence to be analyzed. Only A, C, G, T, and N are allowed. Case insensitive.
/// * `mindn` - Minimum length of the Z-DNA window, in dinucleotides.
/// * `maxdn` - Maximum length of the Z-DNA window, in dinucleotides.
/// * `threshold` - Only DNA windows with a ZH-score above this threshold will be returned.
/// * `wrap` - If true, the sequence will be wrapped around to estimate the ZH-scores of DNA windows ending outside the sequence.
/// * `ssize` - Size of the stream chunks, in nucleotides.
///
/// # Returns
/// An iterable object, where each iteration yields a tuple with 5 vectors describing locally ZH-optimal DNA windows (see `predict` for details).
#[pyfunction(signature = (sequence, mindn, maxdn, threshold, wrap = false, ssize = 1024))]
pub fn stream(
    sequence: Py<PyBytes>,
    mindn: usize,
    maxdn: usize,
    threshold: f64,
    wrap: bool,
    ssize: usize,
) -> PyResult<PyPredictionsStream> {
    if ssize == 0 {
        Err(eyre!("Stream size cannot be zero"))?;
    }

    Ok(PyPredictionsStream {
        engine: zhuntrs::Engine::new(None),
        buffer: zhuntrs::Prediction::default(),

        sequence,
        mindn,
        maxdn,
        threshold,
        wrap,

        ssize,
        index: 0,
    })
}

#[pyclass]
pub struct PyPredictionsStream {
    engine: zhuntrs::Engine,
    buffer: zhuntrs::Prediction,

    sequence: Py<PyBytes>,
    mindn: usize,
    maxdn: usize,
    threshold: f64,
    wrap: bool,

    ssize: usize,
    index: usize,
}

#[pymethods]
impl PyPredictionsStream {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<ZHOutput>> {
        let sequence = self.sequence.as_bytes(py);

        if self.index == sequence.len() {
            return Ok(None);
        }

        let start = self.index;
        let end = (self.index + self.ssize).min(sequence.len());
        debug_assert!(start < end);
        self.index = end;

        // Handle tail extension strategy
        let strategy = if end == sequence.len() {
            match self.wrap {
                true => zhuntrs::TailExtension::WrapAroundFrom(
                    &sequence[..self.maxdn * 2],
                    self.maxdn * 2,
                ),
                false => zhuntrs::TailExtension::None,
            }
        } else {
            zhuntrs::TailExtension::None
        };

        // Run & return predictions
        self.engine.predict(
            &sequence[start..end],
            self.mindn,
            self.maxdn,
            self.threshold,
            strategy,
            &mut self.buffer,
        )?;

        Ok(Some((
            self.buffer.start.clone(),
            self.buffer.end.clone(),
            self.buffer.zhscore.clone(),
            self.buffer
                .sequence
                .iter()
                .map(|x| x.iter().map(|n| n.to_string()).collect())
                .collect(),
            self.buffer
                .conformation
                .iter()
                .map(|x| x.iter().map(|dn| dn.to_string()).collect())
                .collect(),
        )))
    }

    fn consume(&mut self, py: Python<'_>) -> PyResult<ZHOutput> {
        let sequence = self.sequence.as_bytes(py);
        self.ssize = sequence.len() - self.index;
        match self.__next__(py)? {
            None => Err(eyre!("No predictions left"))?,
            Some(result) => Ok(result),
        }
    }
}

#[pymodule]
#[pyo3(name = "zhuntrs")]
pub fn _zhuntrs(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(predict, module)?)?;
    module.add_function(wrap_pyfunction!(stream, module)?)?;
    Ok(())
}
