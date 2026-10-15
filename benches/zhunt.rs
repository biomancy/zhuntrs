#![feature(test)]
extern crate test;

use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::{fs, io};
use test::{Bencher, black_box};

fn load_sequence(seqpath: &str) -> io::Result<Vec<u8>> {
    let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests");
    let seqcontent = fs::read_to_string(directory.join(seqpath))?;
    let sequence = seqcontent.trim();
    Ok(sequence.bytes().collect())
}

#[bench]
fn bench_mtdna_1_6(b: &mut Bencher) -> io::Result<()> {
    let sequence = load_sequence("resources/mtDNA.NC_012920.1.txt")?;

    let engine = zhuntrs::Engine::new(
        NonZeroUsize::new(1).unwrap(),
        NonZeroUsize::new(6).unwrap(),
        0.0,
        false,
    )
    .expect("valid benchmark parameters");
    b.iter(|| black_box(engine.predict(sequence.as_slice())));
    Ok(())
}

#[bench]
fn bench_random_5_12(b: &mut Bencher) -> io::Result<()> {
    let sequence = load_sequence("resources/random.txt")?;

    let engine = zhuntrs::Engine::new(
        NonZeroUsize::new(5).unwrap(),
        NonZeroUsize::new(12).unwrap(),
        0.0,
        false,
    )
    .expect("valid benchmark parameters");
    b.iter(|| black_box(engine.predict(sequence.as_slice())));
    Ok(())
}
