use std::fs;
use std::path::PathBuf;

use itertools::Itertools;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn compare_implementations(seqpath: &str, expected: &str) -> Result<()> {
    let file = PathBuf::from(file!());
    let directory = file.parent().ok_or("Unable to get parent directory")?;

    let seqcontent = fs::read_to_string(directory.join(seqpath))
        .map_err(|err| format!("Unable to read test sequence {seqpath}: {err}"))?;
    let sequence = seqcontent.trim();

    let expected = fs::read_to_string(directory.join(expected))
        .map_err(|err| format!("Unable to read expected result {expected}: {err}"))?;
    let mut lines = expected.lines();

    // Parse algorithm parameters
    let (_, length, minw, maxw) = lines
        .next()
        .ok_or("Unable to get Z-Hunt parameters: target file is empty")?
        .split_whitespace()
        .collect_tuple()
        .ok_or("Unable to parse Z-Hunt parameters: expected 4 values")?;

    let length: usize = length
        .parse()
        .map_err(|err| format!("Unable to parse Z-Hunt parameter length: {err}"))?;
    assert_eq!(
        length,
        sequence.len(),
        "Sequence length mismatch: expected {}, got {}",
        length,
        sequence.len()
    );

    let mindn = minw
        .parse()
        .map_err(|err| format!("Unable to parse Z-Hunt parameter minw: {err}"))?;
    let maxdn = maxw
        .parse()
        .map_err(|err| format!("Unable to parse Z-Hunt parameter maxw: {err}"))?;

    let predicted = zhuntrs::Engine::new(mindn, maxdn, 0.0, true)?.predict(sequence.as_bytes())?;
    for count in [
        predicted.start.len(),
        predicted.end.len(),
        predicted.zhscore.len(),
        predicted.sequence.len(),
        predicted.conformation.len(),
    ] {
        assert_eq!(count, length, "Prediction count mismatch");
    }

    let mut errors = vec![];
    for (ind, line) in lines.enumerate() {
        let (_, _, exp_zhscore, exp_conformation) = line
            .split_whitespace()
            .collect_tuple()
            .ok_or_else(|| format!("Unable to parse expected result in the line: {line}"))?;

        let exp_start = ind;
        let exp_end = exp_start + exp_conformation.len();

        let exp_zhscore: f64 = exp_zhscore
            .parse()
            .map_err(|err| format!("Unable to parse expected ZH-score {exp_zhscore}: {err}"))?;

        let pred_start = predicted.start[ind];
        let pred_end = predicted.end[ind];
        let pred_zhscore = predicted.zhscore[ind];
        let pred_sequence: String = predicted.sequence[ind]
            .iter()
            .map(ToString::to_string)
            .collect();
        let pred_conformation: String = predicted.conformation[ind]
            .iter()
            .map(ToString::to_string)
            .collect();

        let diff = (exp_zhscore - pred_zhscore).abs();
        if exp_start != pred_start
            || exp_end != pred_end
            || exp_conformation != pred_conformation
            || !pred_zhscore.is_finite()
            || diff >= 1e-12
        {
            errors.push(format!(
                "Incompatible lines ({pred_sequence}: {diff}):\n\
                \tExpected: {exp_start} {exp_end} {exp_zhscore} {exp_conformation}\n\
                \tGot: {pred_start} {pred_end} {pred_zhscore} {pred_conformation}"
            ));
        }
    }
    if errors.is_empty() {
        return Ok(());
    }

    let total = errors.len();
    let freq = total as f64 / predicted.start.len() as f64 * 100.0;
    let errors = errors.join("\n");
    Err(format!("Total: {total} ({freq:.2}%)\n{errors}").into())
}

#[test]
fn test_chr_y() -> Result<()> {
    compare_implementations(
        "resources/chrY.NC_000024.10-crop.txt",
        "resources/chrY.NC_000024.10-crop.txt.Z-SCORE",
    )?;
    Ok(())
}

#[test]
fn test_mt_dna() -> Result<()> {
    compare_implementations(
        "resources/mtDNA.NC_012920.1.txt",
        "resources/mtDNA.NC_012920.1.txt.Z-SCORE",
    )?;
    Ok(())
}

#[test]
fn test_hsv1_genome() -> Result<()> {
    compare_implementations(
        "resources/HSV-1.NC_001806.2-crop.txt",
        "resources/HSV-1.NC_001806.2-crop.txt.Z-SCORE",
    )?;
    Ok(())
}

#[test]
fn test_random_sequence() -> Result<()> {
    compare_implementations("resources/random.txt", "resources/random.txt.Z-SCORE")?;
    Ok(())
}

#[test]
fn test_short_sequence() -> Result<()> {
    compare_implementations("resources/short.txt", "resources/short.txt.Z-SCORE")?;
    Ok(())
}

#[test]
fn test_test_sequence() -> Result<()> {
    compare_implementations("resources/test.txt", "resources/test.txt.Z-SCORE")?;
    Ok(())
}
