use pyo3::prelude::*;

mod bayes;
mod components;
mod randomisers;
mod stats;
mod types;
use crate::components::priors::Prior;
// use crate::randomisers::random::make_random;

// use crate::stats::algo::welfords;

use crate::bayes::paired_d::bf_sim_paired;
use crate::types::SamplingRule;

/// bf_sim(effsize, comp, prior, min_n, max_n, step_size, alternative, reps, seed)
/// ---
/// Parameters:
///     effsize: The effect size
///     comp: Type of comparison. Either "t.paired" or "between_t"
///     prior:
#[pyfunction]
#[pyo3(text_signature = "(effsize, comp, prior, sampling_rule, alternative, reps, seed)")]
fn bf_sim(
    effsize: f64,
    comp: &str,
    prior: Prior,
    sampling_rule: SamplingRule,
    alternative: &str,
    reps: i32,
    seed: Option<u64>,
) -> Vec<(usize, f64, i32, Option<f64>, f64, f64, f64)> {
    let d = match comp {
        "t.paired" => Some(bf_sim_paired(
            effsize,
            sampling_rule,
            reps,
            seed,
            prior,
            alternative,
        )),
        _ => unimplemented!(),
    };

    d.unwrap().table()
}

#[pymodule]
fn pyby(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bf_sim, m)?)?;
    Ok(())
}

