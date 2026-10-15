use std::num::NonZeroUsize;
use std::sync::Mutex;

use pyo3::exceptions::{PyOSError, PyOverflowError, PyRuntimeError, PyStopIteration, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyTuple, PyType};

use crate::{Error, zhuntrs};

impl From<Error> for PyErr {
    fn from(error: Error) -> Self {
        match error {
            Error::InvalidWindowRange { .. }
            | Error::InvalidThreshold { .. }
            | Error::InvalidNucleotide { .. }
            | Error::ZeroStreamSize => PyValueError::new_err(error.to_string()),
            Error::WindowLengthOverflow { .. } | Error::InputTooLarge => {
                PyOverflowError::new_err(error.to_string())
            }
            Error::StreamExhausted => PyStopIteration::new_err(error.to_string()),
            Error::StreamPoisoned => PyRuntimeError::new_err(error.to_string()),
            Error::Io(source) => match source.raw_os_error() {
                // The fourth argument preserves winerror on Windows and is
                // ignored elsewhere. OSError selects the appropriate subclass.
                Some(code) => PyOSError::new_err((code, source.to_string(), None::<String>, code)),
                None => source.into(),
            },
        }
    }
}

/// Locally optimal Z-DNA subsequences.
#[pyclass(name = "Prediction", module = "zhuntrs", frozen, get_all)]
#[derive(Debug)]
pub struct PyPrediction {
    /// Zero-based start coordinates.
    start: Py<PyTuple>,
    /// Exclusive end coordinates.
    end: Py<PyTuple>,
    /// Z-Hunt scores.
    zhscore: Py<PyTuple>,
    /// Window sequences in uppercase ASCII.
    sequence: Py<PyTuple>,
    /// Anti/syn conformations for each window.
    conformation: Py<PyTuple>,
}

#[pymethods]
impl PyPrediction {
    #[new]
    fn new(
        py: Python<'_>,
        start: Vec<usize>,
        end: Vec<usize>,
        zhscore: Vec<f64>,
        sequence: Vec<String>,
        conformation: Vec<String>,
    ) -> PyResult<Self> {
        let count = start.len();
        if end.len() != count
            || zhscore.len() != count
            || sequence.len() != count
            || conformation.len() != count
        {
            return Err(PyValueError::new_err(
                "Prediction columns must have equal lengths",
            ));
        }
        Ok(Self {
            start: PyTuple::new(py, start)?.unbind(),
            end: PyTuple::new(py, end)?.unbind(),
            zhscore: PyTuple::new(py, zhscore)?.unbind(),
            sequence: PyTuple::new(py, sequence)?.unbind(),
            conformation: PyTuple::new(py, conformation)?.unbind(),
        })
    }

    fn __len__(&self, py: Python<'_>) -> usize {
        self.start.bind(py).len()
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!("Prediction(windows={})", self.__len__(py))
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> PyResult<bool> {
        Ok(self.start.bind(py).eq(other.start.bind(py))?
            && self.end.bind(py).eq(other.end.bind(py))?
            && self.zhscore.bind(py).eq(other.zhscore.bind(py))?
            && self.sequence.bind(py).eq(other.sequence.bind(py))?
            && self.conformation.bind(py).eq(other.conformation.bind(py))?)
    }

    #[allow(clippy::type_complexity)]
    fn __reduce__<'py>(
        &self,
        py: Python<'py>,
    ) -> (
        Bound<'py, PyType>,
        (
            Py<PyTuple>,
            Py<PyTuple>,
            Py<PyTuple>,
            Py<PyTuple>,
            Py<PyTuple>,
        ),
    ) {
        (
            py.get_type::<Self>(),
            (
                self.start.clone_ref(py),
                self.end.clone_ref(py),
                self.zhscore.clone_ref(py),
                self.sequence.clone_ref(py),
                self.conformation.clone_ref(py),
            ),
        )
    }
}

impl PyPrediction {
    fn from_prediction(py: Python<'_>, prediction: zhuntrs::Prediction) -> PyResult<Self> {
        Ok(Self {
            start: PyTuple::new(py, prediction.start)?.unbind(),
            end: PyTuple::new(py, prediction.end)?.unbind(),
            zhscore: PyTuple::new(py, prediction.zhscore)?.unbind(),
            sequence: PyTuple::new(
                py,
                prediction.sequence.into_iter().map(|window| {
                    window
                        .into_iter()
                        .map(|n| n.to_string())
                        .collect::<String>()
                }),
            )?
            .unbind(),
            conformation: PyTuple::new(
                py,
                prediction.conformation.into_iter().map(|window| {
                    window
                        .into_iter()
                        .map(|dn| dn.to_string())
                        .collect::<String>()
                }),
            )?
            .unbind(),
        })
    }
}

/// Validated, immutable Z-Hunt parameters.
///
/// Window lengths are in dinucleotides. Only scores at least `threshold` are
/// returned. `wrap` enables circular sequences. Instances can be shared across
/// threads and pickled; each prediction call or iterator uses separate buffers.
#[pyclass(name = "Engine", module = "zhuntrs", frozen, eq)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyEngine {
    engine: zhuntrs::Engine,
}

#[pymethods]
impl PyEngine {
    #[new]
    #[pyo3(signature = (mindn, maxdn, threshold = 0.0, wrap = false))]
    fn new(mindn: NonZeroUsize, maxdn: NonZeroUsize, threshold: f64, wrap: bool) -> PyResult<Self> {
        Ok(Self {
            engine: zhuntrs::Engine::new(mindn, maxdn, threshold, wrap)?,
        })
    }

    #[getter]
    fn mindn(&self) -> usize {
        self.engine.mindn().get()
    }

    #[getter]
    fn maxdn(&self) -> usize {
        self.engine.maxdn().get()
    }

    #[getter]
    fn threshold(&self) -> f64 {
        self.engine.threshold()
    }

    #[getter]
    fn wrap(&self) -> bool {
        self.engine.wrap()
    }

    /// Return a Prediction with starts, ends, scores, sequences, and conformations.
    ///
    /// Input must be bytes containing ASCII A, C, G, T, or N (case insensitive).
    /// Coordinates are zero based, with exclusive ends. Empty input yields an
    /// empty Prediction. Computation runs detached from the Python interpreter.
    fn predict(&self, py: Python<'_>, sequence: &[u8]) -> PyResult<PyPrediction> {
        let prediction = py.detach(|| self.engine.predict(sequence))?;
        PyPrediction::from_prediction(py, prediction)
    }

    /// Stream predictions for up to `chunk` start positions per chunk.
    #[pyo3(signature = (sequence, chunk = 1024))]
    fn stream(
        slf: Py<Self>,
        py: Python<'_>,
        sequence: &[u8],
        chunk: usize,
    ) -> PyResult<PyPredictionIterator> {
        let chunk = NonZeroUsize::new(chunk).ok_or(Error::ZeroStreamSize)?;
        PyPredictionIterator::new(py, slf, sequence, chunk)
    }

    fn __repr__(&self) -> String {
        format!(
            "Engine(mindn={}, maxdn={}, threshold={:?}, wrap={})",
            self.mindn(),
            self.maxdn(),
            self.threshold(),
            if self.wrap() { "True" } else { "False" },
        )
    }

    // Reconstruct through the public validating constructor, without mutable
    // state restoration. Working memory is never part of the pickle.
    fn __reduce__<'py>(&self, py: Python<'py>) -> (Bound<'py, PyType>, (usize, usize, f64, bool)) {
        (
            py.get_type::<Self>(),
            (self.mindn(), self.maxdn(), self.threshold(), self.wrap()),
        )
    }
}

/// An iterator of prediction chunks, safe for concurrent consumers.
#[pyclass(name = "PredictionIterator", module = "zhuntrs", frozen)]
pub struct PyPredictionIterator {
    state: Mutex<zhuntrs::IteratorState>,
    engine: Py<PyEngine>,
    chunk: NonZeroUsize,
}

impl PyPredictionIterator {
    fn new(
        py: Python<'_>,
        engine: Py<PyEngine>,
        sequence: &[u8],
        chunk: NonZeroUsize,
    ) -> PyResult<Self> {
        let state = py
            .detach(|| zhuntrs::IteratorState::new(&engine.get().engine, sequence))
            .map_err(PyErr::from)?;
        Ok(Self {
            state: Mutex::new(state),
            engine,
            chunk,
        })
    }

    fn next_output(&self, py: Python<'_>, chunk: NonZeroUsize) -> PyResult<Option<PyPrediction>> {
        let prediction = py.detach(|| -> crate::Result<_> {
            let engine = &self.engine.get().engine;
            let mut state = self.state.lock().map_err(|_| Error::StreamPoisoned)?;
            Ok(state.next(engine, chunk))
        })?;
        prediction
            .map(|prediction| PyPrediction::from_prediction(py, prediction))
            .transpose()
    }
}

#[pymethods]
impl PyPredictionIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self, py: Python<'_>) -> PyResult<Option<PyPrediction>> {
        self.next_output(py, self.chunk)
    }

    /// Return all remaining predictions, or raise StopIteration if exhausted.
    fn consume(&self, py: Python<'_>) -> PyResult<PyPrediction> {
        self.next_output(py, NonZeroUsize::MAX)?
            .ok_or(Error::StreamExhausted)
            .map_err(PyErr::from)
    }
}

#[pymodule(gil_used = false)]
#[pyo3(name = "zhuntrs")]
pub fn _zhuntrs(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPrediction>()?;
    module.add_class::<PyEngine>()?;
    module.add_class::<PyPredictionIterator>()?;
    Ok(())
}
