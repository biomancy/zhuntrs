use std::num::NonZeroUsize;

use zhuntrs::{Engine, Prediction, Result};

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
fn streams_match_predictions_across_boundaries() -> Result<()> {
    let sequence = b"nATGCGCGCGGCATGCAT";
    for length in 0..=sequence.len() {
        let sequence = &sequence[..length];
        for maxdn in 1..=4 {
            for wrap in [false, true] {
                for threshold in [0.0, 1000.0, f64::INFINITY] {
                    let engine = Engine::new(nz(1), nz(maxdn), threshold, wrap)?;
                    let expected = engine.predict(sequence)?;

                    for chunk in [1, 2, 3, 5, 8, 32, usize::MAX] {
                        let chunks = engine
                            .stream(sequence, nz(chunk))?
                            .collect::<Result<Vec<_>>>()?;
                        assert!(chunks.iter().all(|result| result.start.len() <= chunk));
                        assert_eq!(
                            concatenate(chunks),
                            expected,
                            "length={length}, maxdn={maxdn}, wrap={wrap}, threshold={threshold}, chunk={chunk}"
                        );
                    }
                }
            }
        }
    }
    Ok(())
}
