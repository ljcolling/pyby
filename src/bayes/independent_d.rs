// use rayon::prelude::*;
// use itertools::Itertools;

// use crate::make_random;
// use crate::welfords;
// use super::super::stats::summary::between_d;

/*
pub fn bf_sim_independent(
    effsize: f64,
    n_min: i32,
    n_max: i32,
    step_size: i32,
    reps: i32,
    seed: Option<u64>,
    // prior: Prior,
    // alternative: &str,
) -> Vec<(usize, i32, f64, f64, f64, f64, f64, f64, f64)> {
    let ran1: Vec<(usize, i32, f64, f64, f64)> = make_random(effsize, n_max * reps, seed)
        .into_par_iter()
        .chunks(n_max as usize)
        .enumerate()
        .map(|(i, group)| {
            group
                .iter()
                .scan((i, 0f64, 0f64, 0f64), welfords)
                .collect_vec()
        })
        .flatten()
        .filter(|(_, n, _, _, _)| (n >= &n_min) && (*n % step_size) == 0)
        .collect();

    let ran2: Vec<(usize, i32, f64, f64, f64)> = make_random(0.0, n_max * reps, seed)
        .into_par_iter()
        .chunks(n_max as usize)
        .enumerate()
        .map(|(i, group)| {
            group
                .iter()
                .scan((i, 0f64, 0f64, 0f64), welfords)
                .collect_vec()
        })
        .flatten()
        .filter(|(_, n, _, _, _)| (n >= &n_min) && (*n % step_size) == 0)
        .collect();

    ran1.into_par_iter()
        .zip(ran2)
        .map(|((group, n, value1, m1, s1), (_, _, value2, m2, s2))| {
            // calculate the bf value and structure the output here
            let d = between_d(m1, m2, s1, s2, n as f64);
            (group, n, value1, value2, m1, m2, s1, s2, d)
        })
        .collect::<Vec<(usize, i32, f64, f64, f64, f64, f64, f64, f64)>>()
} */
