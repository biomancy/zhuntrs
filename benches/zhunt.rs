#![feature(test)]
extern crate test;

use std::fs;
use std::path::PathBuf;
use test::{black_box, Bencher};

use eyre::{ContextCompat, Result, WrapErr};

use zhuntrs;

fn load_sequence(seqpath: &str) -> Result<Vec<u8>> {
    let file = PathBuf::from(file!());
    let directory = file.parent().wrap_err("Unable to get parent directory")?;

    let seqcontent = fs::read_to_string(directory.join(seqpath))
        .wrap_err(format!("Unable to read test sequence: {seqpath}"))?;
    let sequence = seqcontent.trim();
    Ok(sequence.bytes().collect())
}

#[bench]
fn bench_mtdna_1_6(b: &mut Bencher) -> Result<()> {
    let sequence = load_sequence("resources/mtDNA.NC_012920.1.txt")
        .wrap_err("Unable to load mtDNA sequence")?;

    b.iter(|| black_box(zhuntrs::py::predict(&sequence, 1, 6, 0.0, false)));
    Ok(())
}

#[bench]
fn bench_random_5_12(b: &mut Bencher) -> Result<()> {
    let sequence =
        load_sequence("resources/random.txt").wrap_err("Unable to load random sequence")?;

    b.iter(|| black_box(zhuntrs::py::predict(&sequence, 5, 12, 0.0, false)));
    Ok(())
}
