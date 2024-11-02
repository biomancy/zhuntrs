use std::io::{BufRead, BufReader, Read};

use eyre::{eyre, Result};

use super::dbzed::{dbzed, DiNucleotide, DBZED};
use super::nucleotide::Nucleotide;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TailExtension<'a> {
    WrapAround(usize),
    WrapAroundFrom(&'a [u8], usize),
    None,
}

#[derive(Debug, Clone, Default)]
struct State {
    pub esum: f32,
    pub bzenergy: Vec<f64>,
    pub conformation: Vec<DiNucleotide>,
}

#[derive(Debug, Clone, Default)]
pub struct Prediction {
    pub start: Vec<usize>,
    pub end: Vec<usize>,
    pub zhscore: Vec<f64>,
    pub sequence: Vec<Vec<Nucleotide>>,
    pub conformation: Vec<Vec<DiNucleotide>>,
}

pub struct Engine {
    best: Vec<State>,
    cache: State,

    logcoef: Vec<f64>,
    exponent: Vec<f64>,

    // Parsed nucleotides + index of nucleotides to use for bzenergy calculation
    sequence: Vec<Nucleotide>,
    bzindex: Vec<usize>,

    // Constants computed only once
    bztwist: Vec<f64>,
    expdbzed: [[f64; 17]; 4],
}

impl Engine {
    pub fn new(reserve: Option<usize>) -> Engine {
        let mut expdbzed = [[0.0; 17]; 4];
        let rt = 0.59004; // 0.00198*298
        for i in 0..4 {
            for j in 0..17 {
                expdbzed[i][j] = (-DBZED[i][j] / rt).exp();
            }
        }

        let reserve = reserve.unwrap_or(32);
        Engine {
            best: Vec::with_capacity(reserve),
            cache: State::default(),
            bztwist: Vec::with_capacity(reserve),
            logcoef: Vec::with_capacity(reserve),
            exponent: Vec::with_capacity(reserve),
            sequence: Vec::with_capacity(reserve),
            bzindex: Vec::with_capacity(reserve),
            expdbzed,
        }
    }

    fn reserve(&mut self, maxdn: usize, buffer: &mut Prediction) {
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

        self.sequence.clear();
        self.bzindex.clear();

        self.best.resize(maxdn + 1, State::default());
        for (dn, state) in self.best.iter_mut().enumerate() {
            state.esum = 0.0;
            state.conformation.resize(dn, DiNucleotide::AS);
            state.bzenergy.resize(dn, 0.0);
        }

        buffer.start.clear();
        buffer.end.clear();
        buffer.zhscore.clear();
        buffer.sequence.clear();
        buffer.conformation.clear();
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
        if dl > AVERAGE {
            z
        } else {
            1.0 / z
        }
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
    fn precompute_bzenergy_index(&mut self, maxdn: usize, start: usize) -> Result<()> {
        self.bzindex.clear();

        let mut i = start;
        while i < (start + maxdn * 2) {
            let (c1, c2) = (self.sequence[i], self.sequence[i + 1]);
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
            i += 2;
        }
        Ok(())
    }
    fn input_sequence<R: Read>(&mut self, reader: R, tail: TailExtension) -> Result<()> {
        let mut reader = BufReader::new(reader);
        let mut buffer = String::new();
        while reader.read_line(&mut buffer)? > 0 {
            for c in buffer.bytes() {
                self.sequence.push(c.try_into()?);
            }
            buffer.clear();
        }

        match tail {
            TailExtension::WrapAround(size) => {
                // Wrap the sequence around itself
                let length = self.sequence.len();
                for i in 0..size {
                    let i = i % length;
                    self.sequence.push(self.sequence[i]);
                }
            }
            TailExtension::WrapAroundFrom(sequence, size) => {
                let sequence = sequence[..size]
                    .iter()
                    .map(|&c| c.try_into())
                    .collect::<Result<Vec<Nucleotide>>>()?;

                let length = sequence.len();
                for i in 0..size {
                    let i = i % length;
                    self.sequence.push(sequence[i]);
                }
            }
            TailExtension::None => {}
        }
        Ok(())
    }
    pub fn predict<R: Read>(
        &mut self,
        reader: R,
        mindn: usize,
        maxdn: usize,
        threshold: f64,
        wrap: TailExtension,
        buffer: &mut Prediction,
    ) -> Result<usize> {
        if mindn > maxdn {
            return Err(eyre!(
                "Minimum length of the DNA window must be less than or equal to the maximum length"
            ));
        } else if threshold < 0.0 {
            return Err(eyre!("Threshold must be >= 0"));
        } else if mindn == 0 || maxdn == 0 {
            return Err(eyre!(
                "Minimum/Maximum length of the DNA window must be > 0"
            ));
        }

        self.reserve(maxdn, buffer);
        self.input_sequence(reader, wrap)?;

        let initial_esum = 10.0 * maxdn as f32;
        for i in 0..self.sequence.len() - 2 * maxdn {
            self.precompute_bzenergy_index(maxdn, i)?;

            for best in self.best.iter_mut() {
                best.esum = initial_esum;
            }
            // self.anti_syn_energy_stack(maxdn);
            self.cache.esum = 0.0;
            self.find_optimal_conformation(0, maxdn);

            let (dl, dn) = self.find_optimal_delta_linking(mindn, maxdn);
            let length = dn * 2;

            let probability = self.estimate_zh_score(dl);
            if probability < threshold {
                continue;
            }

            buffer.start.push(i);
            buffer.end.push(i + length);
            buffer.zhscore.push(probability);
            buffer.sequence.push(self.sequence[i..i + length].to_vec());
            buffer.conformation.push(self.best[dn].conformation.clone());
        }
        Ok(buffer.start.len())
    }

    // pub fn stream<'a>(
    //     &'a mut self, sequence: &'a [u8], ssize: Option<usize>,
    //     mindn: usize, maxdn: usize, threshold: f64, wrap: TailExtension<'a>,
    // ) -> PredictionsStream {
    //     let ssize = ssize.unwrap_or(1024);
    //     PredictionsStream::new(self, sequence, ssize, mindn, maxdn, threshold, wrap)
    // }
}

// pub struct PredictionsStream<'a, 'b>
// {
//     engine: &'a mut Engine,
//     sequence: &'a [u8],
//     mindn: usize,
//     maxdn: usize,
//     threshold: f64,
//     wrap: TailExtension<'b>,
//
//     ssize: usize,
//     index: usize,
// }
//
// impl<'a, 'b> PredictionsStream<'a, 'b> {
//     fn new(
//         engine: &'a mut Engine, sequence: &'a [u8], ssize: usize,
//         mindn: usize, maxdn: usize, threshold: f64, wrap: TailExtension<'b>,
//     ) -> PredictionsStream<'a, 'b> {
//         PredictionsStream { engine, sequence, mindn, maxdn, threshold, wrap, ssize, index: 0 }
//     }
//
//     pub fn next(&mut self, buffer: &mut Prediction) -> Result<usize> {
//         if self.index == self.sequence.len() {
//             return Ok(0);
//         }
//
//         let start = self.index;
//         let end = (self.index + self.ssize).min(self.sequence.len());
//         debug_assert!(start < end);
//         self.index = end;
//
//         // Handle tail extension strategy
//         let strategy = if end == self.sequence.len() {
//             match self.wrap {
//                 TailExtension::WrapAround(size) => {
//                     if size < self.sequence.len() {
//                         return Err(eyre!("Wrap around size must be >= sequence length"));
//                     }
//                     TailExtension::WrapAroundFrom(&self.sequence[..size], size)
//                 }
//                 TailExtension::WrapAroundFrom(sequence, size) => TailExtension::WrapAroundFrom(sequence, size),
//                 TailExtension::None => TailExtension::None,
//             }
//         } else {
//             TailExtension::None
//         };
//
//         // Run & return the total number of predictions
//         self.engine.predict(
//             &self.sequence[start..end], self.mindn, self.maxdn, self.threshold, strategy, buffer,
//         )
//     }
//
//     pub fn consume(&mut self, buffer: &mut Prediction) -> Result<usize> {
//         self.ssize = self.sequence.len() - self.index;
//         self.next(buffer)
//     }
// }

// fn anti_syn_energy_dynamic_programming(&mut self, maxdn: usize) {
//     const CONFORMATIONS: [DiNucleotide; 2] = [DiNucleotide::AS, DiNucleotide::SA];
//     self.anti_syn_energy_stack(maxdn);
//
//     // Optimum energy for structures ending with particular conformation
//     let mut prev_energies = CONFORMATIONS.iter().map(
//         |conformation| DBZED[dbzed(None, *conformation)][self.bzindex[0]]
//     ).collect::<Vec<f64>>();
//     let mut next_energies = vec![0.0; CONFORMATIONS.len()];
//
//     let mut prev_score = CONFORMATIONS.iter().map(
//         |conformation| self.expdbzed[dbzed(None, *conformation)][self.bzindex[0]]
//     ).collect::<Vec<f64>>();
//     let mut next_score = vec![0.0; CONFORMATIONS.len()];
//
//     // Backtrace for the best conformation ( => which conformation was used for the previous dinucleotide)
//     let mut backtrace = (0..CONFORMATIONS.len()).into_iter().map(
//         |ind| {
//             let mut cache = Vec::with_capacity(maxdn);
//             cache.push(ind);
//             cache
//         }
//     ).collect::<Vec<Vec<usize>>>();
//
//     let mut history = Vec::new();
//
//     // Calculate cache for all dinucleotide
//     for dinuc in 1..maxdn {
//         // Next dinucleotide can be either AS or SA
//         for (next_ind, next_conf) in CONFORMATIONS.iter().enumerate() {
//             let mut best = (f64::MAX, 0, f64::MAX); // Placeholder
//
//             // Previous dinucleotide can be either AS or SA
//             for (prev_ind, prev_conf) in CONFORMATIONS.iter().enumerate() {
//                 let i = dbzed(Some(*prev_conf), *next_conf);
//                 let energy = prev_energies[prev_ind] + DBZED[i][self.bzindex[dinuc]];
//                 let score = energy +  self.expdbzed[i][self.bzindex[dinuc]];
//
//                 if score < best.2 {
//                     best = (energy, prev_ind, score);
//                 }
//             }
//
//             next_energies[next_ind] = best.0;
//             backtrace[next_ind].push(best.1);
//             next_score[next_ind] = best.2;
//         }
//         history.push(prev_energies.clone());
//         std::mem::swap(&mut prev_energies, &mut next_energies);
//         std::mem::swap(&mut prev_score, &mut next_score);
//     }
//     history.push(prev_energies.clone());
//     let energies = prev_energies;
//
//
//     // Select the best conformation
//     let best = energies
//         .iter()
//         .enumerate()
//         .min_by(|(_, a), (_, b)| a.partial_cmp(&b).unwrap())
//         .expect("Failed to find the best conformation").0;
//
//     let best_esum = energies[best];
//
//     // Reconstruct the best conformation
//     let mut current = best;
//     for dn in (0..maxdn).rev() {
//         self.conformation[dn] = CONFORMATIONS[current];
//
//         // Move to the previous dinucleotide
//         current = backtrace[current][dn];
//     }
//
//     // Fill the bzenergy array
//     let ind = dbzed(None, self.conformation[0]);
//     self.bzenergy[0] = self.expdbzed[ind][self.bzindex[0]];
//
//     for dn in 1..maxdn {
//         let i = dbzed(Some(self.conformation[dn - 1]), self.conformation[dn]);
//         self.bzenergy[dn] = self.expdbzed[i][self.bzindex[dn]];
//     }
//
//     // Fill pseudo bzenergy array
//     let ind = dbzed(None, self.best_conformation[0]);
//     let mut reconstructed_esum= DBZED[ind][self.bzindex[0]];
//     let mut reconstructed_bzenergy = vec![self.expdbzed[ind][self.bzindex[0]]];
//     for dn in 1..maxdn {
//         let i = dbzed(Some(self.best_conformation[dn - 1]), self.best_conformation[dn]);
//         reconstructed_bzenergy.push(self.expdbzed[i][self.bzindex[dn]]);
//         reconstructed_esum += DBZED[i][self.bzindex[dn]];
//     }
//
//     for i in 0..maxdn {
//         let a = self.bzenergy[i];
//         let b = self.best_bzenergy[i];
//         // println!("{} {} {}", self.conformation[i], self.bzenergy[i], self.best_bzenergy[i]);
//         assert!(true);
//     }
//
//     let best_bz = self.best_bzenergy.iter().sum::<f64>();
//     let bz = self.bzenergy.iter().sum::<f64>();
//
//
//     self.update_optimum(maxdn, best_esum, self.conformation.clone());
// }

// fn anti_syn_energy_stack(&mut self, maxdn: usize) {
//     // It's impossible to rewrite this function using linear programming, because the original implementation
//     // favored a path with as many AS states as possible.
//     // However, there are more than one combination of AS and SA states that can lead to the same overall esum.
//     // But having identical esum != having identical ZH score. That's where the problem lies.
//
//     let mut stack: Vec<(usize, f64, DiNucleotide)> = Vec::new();
//     // push the initial state to the stack (dinucleotide, esum, next_state)
//     stack.push((0, 0.0, DiNucleotide::AS));
//
//     while let Some((dn, mut esum, state)) = stack.pop() {
//         if esum < self.best[dn].esum {
//             self.best[dn].esum = esum;
//             self.best[dn].conformation.copy_from_slice(&self.cache.conformation[..dn]);
//             self.best[dn].bzenergy.copy_from_slice(&self.cache.bzenergy[..dn]);
//         }
//         if dn == maxdn {
//             continue;
//         }
//
//         match state {
//             DiNucleotide::AS => {
//                 // Explore alternative path on return
//                 stack.push((dn, esum, DiNucleotide::SA));
//
//                 // Explore current path first
//                 self.cache.conformation[dn] = DiNucleotide::AS;
//                 let i = dbzed(
//                     dn.checked_sub(1).and_then(|idx| Some(self.cache.conformation[idx])),
//                     DiNucleotide::AS,
//                 );
//                 esum += DBZED[i][self.bzindex[dn]];
//                 self.cache.bzenergy[dn] = self.expdbzed[i][self.bzindex[dn]];
//                 stack.push((dn + 1, esum, DiNucleotide::AS));
//             }
//             DiNucleotide::SA => {
//                 // Explore SA conformation, AS conformation was already explored here
//                 self.cache.conformation[dn] = DiNucleotide::SA;
//                 let i = dbzed(
//                     dn.checked_sub(1).and_then(|idx| Some(self.cache.conformation[idx])),
//                     DiNucleotide::SA,
//                 );
//                 esum += DBZED[i][self.bzindex[dn]];
//                 self.cache.bzenergy[dn] = self.expdbzed[i][self.bzindex[dn]];
//
//                 // Move to the next dinucleotide and try AS first
//                 stack.push((dn + 1, esum, DiNucleotide::AS));
//             }
//         }
//     }
// }
