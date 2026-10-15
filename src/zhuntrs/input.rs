use std::io::{BufReader, Read};
use std::num::NonZeroUsize;
use std::ops::Range;

use crate::{Error, Result};

use super::nucleotide::Nucleotide;

/// Owns a fully validated sequence with all required window context.
pub(crate) struct Input {
    buffer: Vec<Nucleotide>,
    window_size: NonZeroUsize,
    total_windows: usize,
}

/// Contiguous context for consecutive candidate starts in the original sequence.
pub(crate) struct Windows<'a> {
    pub(crate) bases: &'a [Nucleotide],
    pub(crate) starts: Range<usize>,
}

impl Input {
    /// Read, validate, and append all circular context.
    pub(crate) fn new<R: Read>(reader: R, window_size: NonZeroUsize, wrap: bool) -> Result<Self> {
        let mut sequence = BufReader::new(reader)
            .bytes()
            .map(|byte| Nucleotide::try_from(byte?))
            .collect::<Result<Vec<_>>>()?;

        let length = sequence.len();
        let context = window_size.get() - 1;
        let total_windows = if wrap {
            length
        } else {
            length.saturating_sub(context)
        };

        if wrap && length > 0 {
            length
                .checked_add(context)
                .filter(|&length| length <= isize::MAX as usize)
                .ok_or(Error::InputTooLarge)?;
            sequence.reserve_exact(context);
            for index in 0..context {
                let base = sequence[index % length];
                sequence.push(base);
            }
        }

        Ok(Self {
            buffer: sequence,
            window_size,
            total_windows,
        })
    }

    pub fn total_windows(&self) -> usize {
        self.total_windows
    }

    /// Borrow at most `upto` consecutive windows starting at `start`. A request beyond the end returns `None`.
    pub(crate) fn windows(&self, start: usize, upto: NonZeroUsize) -> Option<Windows<'_>> {
        let limit = self.total_windows();
        if start >= limit {
            return None;
        }

        // Clamp before adding so even usize::MAX is a valid requested count.
        let end = start + upto.get().min(limit - start);
        // Construction guarantees all eligible starts have this much context.
        let input_end = end + (self.window_size.get() - 1);

        Some(Windows {
            bases: &self.buffer[start..input_end],
            starts: start..end,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nz(value: usize) -> NonZeroUsize {
        NonZeroUsize::new(value).unwrap()
    }

    #[test]
    fn windows_match_explicit_examples() -> Result<()> {
        use Nucleotide::{A, C, G, N, T};

        #[derive(Debug)]
        struct Case {
            sequence: &'static [u8],
            width: usize,
            wrap: bool,
            total: usize,
            request: (usize, usize),
            expected: Option<(Range<usize>, &'static [Nucleotide])>,
        }

        let cases = [
            Case {
                sequence: b"acgtN",
                width: 4,
                wrap: false,
                total: 2,
                request: (1, 1),
                expected: Some((1..2, &[C, G, T, N])),
            },
            Case {
                sequence: b"ACGT",
                width: 4,
                wrap: true,
                total: 4,
                request: (2, 2),
                expected: Some((2..4, &[G, T, A, C, G])),
            },
            Case {
                sequence: b"AC",
                width: 6,
                wrap: true,
                total: 2,
                request: (1, 1),
                expected: Some((1..2, &[C, A, C, A, C, A])),
            },
            Case {
                sequence: b"AC",
                width: 4,
                wrap: false,
                total: 0,
                request: (0, 1),
                expected: None,
            },
        ];

        for case in cases {
            let input = Input::new(case.sequence, nz(case.width), case.wrap)?;
            assert_eq!(input.total_windows(), case.total, "{case:?}");
            assert_eq!(
                input
                    .windows(case.request.0, nz(case.request.1))
                    .map(|view| (view.starts, view.bases)),
                case.expected,
                "{case:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn window_views_can_be_borrowed_simultaneously() -> Result<()> {
        use Nucleotide::{A, C, G, T};

        let input = Input::new(b"ACGT".as_slice(), nz(4), true)?;
        let later = input.windows(2, nz(2)).unwrap();
        let earlier = input.windows(0, nz(2)).unwrap();
        assert_eq!(earlier.starts, 0..2);
        assert_eq!(earlier.bases, &[A, C, G, T, A]);
        assert_eq!(later.starts, 2..4);
        assert_eq!(later.bases, &[G, T, A, C, G]);
        Ok(())
    }

    #[test]
    fn empty_input_never_needs_context_even_with_the_largest_width() -> Result<()> {
        for wrap in [false, true] {
            let input = Input::new(b"".as_slice(), NonZeroUsize::MAX, wrap)?;
            assert_eq!(input.total_windows(), 0);
            assert!(input.windows(0, NonZeroUsize::MAX).is_none());
        }
        Ok(())
    }

    #[test]
    fn unaddressable_context_is_rejected_during_construction() {
        for width in [isize::MAX as usize, usize::MAX] {
            assert!(matches!(
                Input::new(b"AC".as_slice(), nz(width), true),
                Err(Error::InputTooLarge)
            ));
        }
    }
}
