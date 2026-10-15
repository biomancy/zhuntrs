use std::sync::OnceLock;

use super::dbzed::{DBZED, DiNucleotide, dbzed};
use super::input::Windows;
use super::nucleotide::Nucleotide;
use super::{Engine, Prediction};

#[derive(Debug, Clone, Default)]
struct State {
    pub esum: f32,
    pub bzenergy: Vec<f64>,
    pub conformation: Vec<DiNucleotide>,
}

pub(crate) struct Predictor {
    best: Vec<State>,
    cache: State,

    logcoef: Vec<f64>,
    exponent: Vec<f64>,

    // Indices used for bzenergy calculation
    bzindex: Vec<usize>,

    // Constants computed only once
    bztwist: Vec<f64>,
    expdbzed: &'static [[f64; 17]; 4],
}

impl Predictor {
    pub(crate) fn new() -> Self {
        static EXPDBZED: OnceLock<[[f64; 17]; 4]> = OnceLock::new();
        let expdbzed = EXPDBZED.get_or_init(|| {
            let mut values = [[0.0; 17]; 4];
            let rt = 0.59004; // 0.00198*298
            for i in 0..4 {
                for j in 0..17 {
                    values[i][j] = (-DBZED[i][j] / rt).exp();
                }
            }
            values
        });

        Self {
            best: Vec::new(),
            cache: State::default(),
            logcoef: Vec::new(),
            exponent: Vec::new(),
            bzindex: Vec::new(),
            bztwist: Vec::new(),
            expdbzed,
        }
    }

    fn reserve(&mut self, maxdn: usize) {
        self.bztwist.clear();
        self.bztwist.reserve(maxdn);
        let mut ab = 0.4 + 0.4; // b + b with b = 0.4;
        for _ in 0..maxdn {
            ab += 0.357; // a with a = 0.357
            self.bztwist.push(ab);
        }

        // Clear and resize all arrays
        self.logcoef.resize(maxdn, 0.0);
        self.exponent.resize(maxdn, 0.0);

        self.cache.esum = 0.0;
        self.cache.bzenergy.resize(maxdn, 0.0);
        self.cache.conformation.resize(maxdn, DiNucleotide::AS);

        self.bzindex.clear();

        self.best.resize(maxdn + 1, State::default());
        for (dn, state) in self.best.iter_mut().enumerate() {
            state.esum = 0.0;
            state.conformation.resize(dn, DiNucleotide::AS);
            state.bzenergy.resize(dn, 0.0);
        }
    }

    fn estimate_zh_score(&self, dl: f64) -> f64 {
        /* calculate the probability of the value 'dl' in a Gaussian distribution */
        /* from "Data Reduction and Error Analysis for the Physical Science" */
        /* Philip R. Bevington, 1969, McGraw-Hill, Inc */
        const AVERAGE: f64 = 29.6537135;
        const STDV: f64 = 2.71997;
        const SQRT_2: f64 = std::f64::consts::FRAC_1_SQRT_2; // 1/sqrt(2)
        const SQRT_PI: f64 = 0.564189583546; // 1/sqrt(pi)

        let mut z = (dl - AVERAGE).abs() / STDV;
        let mut x = z * SQRT_2;
        let y = SQRT_PI * (-x * x).exp();

        z *= z;
        let mut k = 1.0;
        let mut sum = 0.0;

        loop {
            sum += x;
            k += 2.0;

            x *= z / k;
            if sum + x <= sum {
                //
                break;
            }
        }

        // Z-scores
        // z - area under the curve (from -inf to z || from z to +inf)
        // dl > average => 1 - cdf (just right tail)
        // dl <= average => 1 / cdf (1 / left tail)
        z = 0.5 - (y * sum); // probability of each tail
        if dl > AVERAGE { z } else { 1.0 / z }
    }

    fn find_root<F>(&mut self, left: f64, right: f64, tole: f64, func: F) -> f64
    where
        F: Fn(&mut Self, f64) -> f64,
    {
        let fleft = func(self, left);

        // If the function has the same sign at both ends of the interval, then the root is not bracketed
        if fleft * func(self, right) >= 0.0 {
            return right;
        }

        // Move from the negative end to the positive end
        let (mut x, mut dx) = if fleft < 0.0 {
            (left, right - left)
        } else {
            (right, left - right)
        };

        // Halve the interval until the root is bracketed
        loop {
            dx *= 0.5;

            let xmid = x + dx;
            let fmid = func(self, xmid);

            if fmid <= 0.0 {
                x = xmid;
            }

            if dx.abs() <= tole {
                break;
            }
        }
        x
    }

    fn find_optimal_conformation(&mut self, dn: usize, maxdn: usize) {
        // Update the optimum
        if self.cache.esum < self.best[dn].esum {
            self.best[dn].esum = self.cache.esum;
            self.best[dn]
                .conformation
                .copy_from_slice(&self.cache.conformation[..dn]);
            self.best[dn]
                .bzenergy
                .copy_from_slice(&self.cache.bzenergy[..dn]);
        }
        // End of the recursion
        if dn == maxdn {
            return;
        }

        let baseline = self.cache.esum;

        // Explore AS conformation
        self.cache.conformation[dn] = DiNucleotide::AS;
        let i: usize = if dn == 0 {
            dbzed(None, DiNucleotide::AS)
        } else {
            dbzed(Some(self.cache.conformation[dn - 1]), DiNucleotide::AS)
        };
        let energy = DBZED[i][self.bzindex[dn]] as f32;

        self.cache.esum += energy;
        self.cache.bzenergy[dn] = self.expdbzed[i][self.bzindex[dn]];
        self.find_optimal_conformation(dn + 1, maxdn);
        self.cache.esum -= energy;

        // Explore SA conformation
        self.cache.conformation[dn] = DiNucleotide::SA;
        let i = if dn == 0 {
            dbzed(None, DiNucleotide::SA)
        } else {
            dbzed(Some(self.cache.conformation[dn - 1]), DiNucleotide::SA)
        };

        // This will NOT work correctly because of the way floats are handled in C.
        // let energy = DBZED[i][self.bzindex[dn]] as f32;
        // self.cache.esum += energy;

        // Equivalent to the original code behaviour. Adding f64 to f32 in C results in:
        // 1. f32 is promoted to f64
        // 2. f64 is added to f64
        // 3. f64 is demoted to f32
        self.cache.esum = (self.cache.esum as f64 + DBZED[i][self.bzindex[dn]]) as f32;

        self.cache.bzenergy[dn] = self.expdbzed[i][self.bzindex[dn]];
        self.find_optimal_conformation(dn + 1, maxdn);

        self.cache.esum = baseline;
    }

    fn find_optimal_delta_linking(&mut self, mindn: usize, maxdn: usize) -> (f64, usize) {
        let mut best = (50.0, 0);

        for dn in mindn..=maxdn {
            self.cache.bzenergy.clear();
            self.cache.bzenergy.resize(dn, 1.0);

            for i in 0..dn {
                let mut sum = 0.0;
                for j in 0..dn - i {
                    self.cache.bzenergy[j] *= self.best[dn].bzenergy[i + j];
                    sum += self.cache.bzenergy[j];
                }
                self.logcoef[i] = sum.ln();
            }

            let dl = self.find_root(10.0, 50.0, 0.001, |s, dl| s.delta_linking(dl, dn));
            if dl < best.0 {
                best = (dl, dn);
            }
        }
        best
    }

    fn delta_linking(&mut self, dl: f64, dinucleotides: usize) -> f64 {
        const K_RT: f64 = -0.2521201; // -1100/4363
        const SIGMA: f64 = 16.94800353; // 10/RT
        const EXPLIMIT: f64 = -600.0;

        let mut expmini = 0.0;
        for i in 0..dinucleotides {
            let z = self.logcoef[i] + K_RT * (dl - self.bztwist[i]).powi(2);
            self.exponent[i] = z;
            if z < expmini {
                expmini = z;
            }
        }
        expmini = if expmini < EXPLIMIT {
            EXPLIMIT - expmini
        } else {
            0.0
        };

        let mut sump = 0.0;
        let mut sumq = 0.0;
        for i in 0..dinucleotides {
            let z = (self.exponent[i] + expmini).exp();
            sumq += z;
            sump += self.bztwist[i] * z;
        }
        sumq += (K_RT * dl * dl + SIGMA + expmini).exp();

        let deltatwist = 0.357 / 2.0 * dinucleotides as f64;
        deltatwist - sump / sumq
    }

    fn precompute_bzenergy_index(&mut self, sequence: &[Nucleotide]) {
        let (pairs, remainder) = sequence.as_chunks::<2>();
        debug_assert!(
            remainder.is_empty(),
            "Prediction windows must contain an even number of nucleotides"
        );
        self.bzindex.clear();

        for &[c1, c2] in pairs {
            let idx = match c1 {
                Nucleotide::A => match c2 {
                    Nucleotide::A => 0,
                    Nucleotide::T => 1,
                    Nucleotide::G => 2,
                    Nucleotide::C => 3,
                    Nucleotide::N => 16,
                },
                Nucleotide::T => match c2 {
                    Nucleotide::A => 4,
                    Nucleotide::T => 5,
                    Nucleotide::G => 6,
                    Nucleotide::C => 7,
                    Nucleotide::N => 16,
                },
                Nucleotide::G => match c2 {
                    Nucleotide::A => 8,
                    Nucleotide::T => 9,
                    Nucleotide::G => 10,
                    Nucleotide::C => 11,
                    Nucleotide::N => 16,
                },
                Nucleotide::C => match c2 {
                    Nucleotide::A => 12,
                    Nucleotide::T => 13,
                    Nucleotide::G => 14,
                    Nucleotide::C => 15,
                    Nucleotide::N => 16,
                },
                Nucleotide::N => 16,
            };
            self.bzindex.push(idx);
        }
    }

    pub(crate) fn predict(&mut self, engine: &Engine, windows: Windows<'_>) -> Prediction {
        let Windows { bases, starts } = windows;
        let mut output = Prediction::default();
        let (mindn, maxdn, threshold) = (
            engine.mindn().get(),
            engine.maxdn().get(),
            engine.threshold(),
        );
        self.reserve(maxdn);

        let window_size = engine.window_size().get();
        let initial_esum = 10.0 * maxdn as f32;
        for (i, start) in starts.enumerate() {
            self.precompute_bzenergy_index(&bases[i..i + window_size]);

            for best in self.best.iter_mut() {
                best.esum = initial_esum;
            }
            self.cache.esum = 0.0;
            self.find_optimal_conformation(0, maxdn);

            let (dl, dn) = self.find_optimal_delta_linking(mindn, maxdn);
            let length = dn * 2;

            let probability = self.estimate_zh_score(dl);
            if probability < threshold {
                continue;
            }

            output.start.push(start);
            output.end.push(start + length);
            output.zhscore.push(probability);
            output.sequence.push(bases[i..i + length].to_vec());
            output.conformation.push(self.best[dn].conformation.clone());
        }
        output
    }
}
