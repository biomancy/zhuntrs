use std::fmt::{Display, Formatter};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[repr(u8)]
pub enum DiNucleotide {
    AS = 0,
    SA = 1,
}

impl Display for DiNucleotide {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            DiNucleotide::AS => "AS".to_string(),
            DiNucleotide::SA => "SA".to_string(),
        };
        write!(f, "{}", str)
    }
}

#[rustfmt::skip]
pub const DBZED: [[f64; 17]; 4] = [
    // AS-AS
    [4.40, 6.20, 3.40, 5.20, 2.50, 4.40, 1.40, 3.30, 3.30, 5.20, 2.40, 4.20, 1.40, 3.40, 0.66, 2.40, 4.26, ],
    // SA-SA
    [4.40, 2.50, 3.30, 1.40, 6.20, 4.40, 5.20, 3.40, 3.40, 1.40, 2.40, 0.66, 5.20, 3.30, 4.20, 2.40, 4.26, ],
    // AS-SA
    [6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00, 4.00, 4.26, ],
    // SA-AS
    [6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00, 4.00, 4.26, ],
];

pub fn dbzed(prevconf: Option<DiNucleotide>, curconf: DiNucleotide) -> usize {
    match prevconf {
        None => match curconf {
            DiNucleotide::AS => 0,
            DiNucleotide::SA => 1,
        },
        Some(DiNucleotide::AS) => match curconf {
            DiNucleotide::AS => 0,
            DiNucleotide::SA => 3,
        },
        Some(DiNucleotide::SA) => match curconf {
            DiNucleotide::SA => 1,
            DiNucleotide::AS => 2,
        },
    }
}
