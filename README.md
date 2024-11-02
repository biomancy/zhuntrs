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

### Usage

Here is a simple example of how to use the library in Python:

```python
import zhuntrs

sequence = "nATGCGCGCGGCATGC".encode("ASCII")  # The target sequence (only ATGCN are accepted, case insensitive)
mindn, maxdn = 3, 6  # Minimum and maximum length of evaluated DNA windows in dinucleotides
threshold = 1_000  # Windows with a ZH-score below this threshold will not be reported

# Make a prediction
start, end, score, window, conformation = zhuntrs.predict(sequence, mindn, maxdn, threshold)

# start: start coordinates of evaluated DNA windows
assert start == [0, 1, 2]

# end: end coordinates of DNA windows starting at the above coordinates and 
#      having the highest ZH-score among [mindn; maxdn] windows
assert end == [10, 11, 10]

# score: ZH-score for each window
assert [int(zh) for zh in score] == [2220, 1057, 3428]

# window: DNA sequence for each window:
#         window[i] = sequence[start[i]:end[i]]
assert window == ['NATGCGCGCG', 'ATGCGCGCGG', 'TGCGCGCG']

# conformation: DNA conformation (anti/syn) for each window
assert conformation == ['ASASASASAS', 'SASASASASA', 'ASASASAS']
```

Note: the default behaviour of the original implementation is to "wrap around" windows at the end of the sequence.
It can be enabled by passing `wrap=True` to the `predict` function.

The library also supports stream predictions, wherein the sequence is evaluated slice by slice to minimize the memory
footprint and allow interruption when required. This feature is especially useful when processing large sequences:

```python
slice_size = 512  # Size of each slice, default is 1024

for start, end, score, window, conformation in zhuntrs.stream(sequence, mindn, maxdn, threshold, ssize=slice_size):
    length = len(start)
    # The last slice may be shorter than slice_size but not empty
    assert 0 < length <= slice_size
    # All arrays have the same length and semantics as before
    assert len(start) == len(end) == len(score) == len(window) == len(conformation)
    ...
```

This mode is equivalent to calling `zhuntrs.predict` on the entire sequence and processing the output slice by slice,
where each slice has a size equal to `ssize`, except, perhaps, for the last one. Note that the 'wrapping around' of the
sequence in this mode occurs as usual, i.e., only once, and all coordinates are provided relative to the original
sequence.

### Compatibility

This version has been tested against the original C implementation compiled with all compile optimizations turned off.
The output for human mtDNA, a random 1kb sequence, and subsequences of human chrY or HSV-1 genome matches the output of
the C implementation (see tests; 12 digits after the decimal point are identical).

Reference outputs were generated using the `tests/resources/make-references.sh` script. The reference implementation
used is also provided in the same folder (`zhunt3-reference.c`).

### Development

Publishing to PyPI is done automatically by the CI pipeline. To publish a new version, simply create a new tag and push
it to the repository. The CI pipeline will take care of the rest.
