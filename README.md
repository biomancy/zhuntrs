## Z-Hunt[rs]

This repository contains a Rust port of the [Z-Hunt algorithm](https://pubmed.ncbi.nlm.nih.gov/1601856/), which is
designed to predict the propensity of DNA to adopt a left-handed Z-DNA conformation. Unlike the original implementation,
this version is intended to be used as a standalone Python library.

This Rust port is fully compatible with the original C implementation. Please refer to
the [Compatibility](#compatibility) section for more details.

### Installation

The library can be installed from PyPI using `pip`:

```shell
pip install zhuntrs
```

Wheels are built and tested for CPython 3.10–3.15 on Linux, Windows, and macOS (x86-64),
including the supported free-threaded builds of Python 3.14 and 3.15.

### Usage

Here is a simple example of how to use the library in Python:

```python
import zhuntrs

sequence = "nATGCGCGCGGCATGC".encode("ASCII")  # The target sequence (only ATGCN are accepted, case insensitive)
mindn, maxdn = 3, 6  # Minimum and maximum length of evaluated DNA windows in dinucleotides
threshold = 1_000  # Windows with a ZH-score below this threshold will not be reported

# Validate and freeze the parameters, then predict any number of sequences
engine = zhuntrs.Engine(mindn, maxdn, threshold)
prediction: zhuntrs.Prediction = engine.predict(sequence)

assert len(prediction) == 3  # Number of reported windows

# Zero-based starts and exclusive ends of the locally optimal windows
assert prediction.start == (0, 1, 2)
assert prediction.end == (10, 11, 10)

# Scores, sequences, and anti/syn conformations of those windows
assert [int(score) for score in prediction.zhscore] == [2220, 1057, 3428]
assert prediction.sequence == ('NATGCGCGCG', 'ATGCGCGCGG', 'TGCGCGCG')
assert prediction.conformation == ('ASASASASAS', 'SASASASASA', 'ASASASAS')
```

`Prediction` is a Python extension class with read-only `start`, `end`, `zhscore`,
`sequence`, and `conformation` columns. Each column is an immutable tuple constructed once;
repeated attribute access returns the same tuple without copying its contents.

The original implementation treats sequences as circular. Enable this with
`zhuntrs.Engine(mindn, maxdn, threshold, wrap=True)`.

The engine supports streaming predictions for large sequences. It reads and validates the
entire sequence, then scores it in chunks so results can be processed incrementally:

```python
slice_size = 512  # Size of each slice, default is 1024

for batch in engine.stream(sequence, chunk=slice_size):
    assert isinstance(batch, zhuntrs.Prediction)
    # Filtering can produce an empty prediction chunk
    assert 0 <= len(batch) <= slice_size
    assert len(batch.start) == len(batch.end) == len(batch.zhscore)
    assert len(batch.sequence) == len(batch.conformation) == len(batch)
    ...
```

Concatenating the streamed results is equivalent to calling `engine.predict` on the entire sequence.

An engine can be shared between threads, including on free-threaded Python. Multiple threads can also 
consume one iterator; its calls are serialized so each chunk is returned once.

Engines and predictions support pickling:

```python
import pickle

restored = pickle.loads(pickle.dumps(engine))
assert restored == engine
assert restored.predict(sequence) == engine.predict(sequence)
assert pickle.loads(pickle.dumps(prediction)) == prediction
```

The wheel includes type stubs for engines, predictions, and prediction iterators.

### Compatibility

This version has been tested against the original C implementation compiled with all compile optimizations turned off.
The output for human mtDNA, a random 1kb sequence, and subsequences of human chrY or HSV-1 genome matches the output of
the C implementation (see tests; 12 digits after the decimal point are identical).

Reference outputs were generated using the `tests/resources/make-references.sh` script. The reference implementation used is also provided in the same folder (`zhunt-reference.c`). Note that the license
of the reference implementation (MIT).

### Development

Use [Pixi](https://pixi.sh) to install the pinned development environment and run checks:

```shell
pixi run --locked rs-ci     # Rust formatting, compilation, Clippy, and tests
pixi run --locked py-ci     # Build the extension and check Python bindings and typing
pixi run --locked ci        # The full CI suite
```

The default environment uses free-threaded Python 3.15. Named environments `py-3-10`
through `py-3-15`, plus `py-3-14t` and `py-3-15t`, support the wheel matrix on Linux,
macOS, and Windows.

Publishing to PyPI is done automatically by the CI pipeline after checks and tests. To publish a new version, update the version in `Cargo.toml`, then create and push the matching tag (for example, `v0.0.7`).
