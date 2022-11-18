use crate::SamplingRule;
use crate::Prior;
use crate::types::SimTable;
use crate::types::SimRow;
use crate::randomisers::random::make_random;
use crate::stats::algo::welfords;
use crate::stats::teststats::one_sample_t;
use crate::components::predictive::create_predictive;
use crate::components::likelihoods::noncentral_d_likelihood;

use crate::components::priors::Function;

use itertools::Itertools;
use gkquad::single::Integrator;

use std::f64::{INFINITY, NEG_INFINITY};

pub fn bf_sim_paired(
    effsize: f64,
    sampling_rule: SamplingRule,
    reps: i32,
    seed: Option<u64>,
    prior: Prior,
    alternative: &str,
) -> SimTable {
    let (n_min, n_max, step_size) = sampling_rule.unpack();
    // make_random1(effsize, n_max , reps, seed)
    // let d: Vec<(usize, f64, i32, Option<f64>, f64, f64, f64)> =
    let d: Vec<SimRow> = make_random(effsize, n_max * reps, seed)
        .chunks(n_max as usize)
        // .iter()
        .enumerate()
        .map(|(i, group)| {
            group
                .iter()
                .scan((i, 0f64, 0f64, 0f64), welfords)
                .collect_vec()
        })
        .flatten()
        .filter(|(_, n, _, _, _)| (n >= &n_min) && (*n % step_size) == 0)
        .map(
            |(id, n, _, mean, sd)| {
                let bf = if n < 2 {
                    Some(0f64)
                } else {
                    bf_onesample_default_t((mean / sd) as f64, n, &prior, alternative)
                };
                let emp_es = mean / sd;
                let se = sd / (n as f64).sqrt();
                let t = mean / se;
                let pvalue = one_sample_t(t, n, 0);
                SimRow {
                    id,
                    effsize,
                    n,
                    bf,
                    emp_es,
                    t,
                    pvalue,
                }
                // (id, effsize, n, bf, emp_es, t, pvalue)
            },
        )
        .collect();

    SimTable {
        rows: d.into_iter().collect(),
    }
}


pub fn bf_onesample_default_t(d: f64, n: i32, prior: &Prior, alternative: &str) -> Option<f64> {
    // let d = mean / sd;
    let likelihood = noncentral_d_likelihood(d, n as f64);

    let h1_prior = prior.function(alternative);
    let null = likelihood(0.0);
    let h1 = create_predictive(likelihood, h1_prior);
    let alt = Integrator::new(h1)
        // .algorithm(QAGS::new())
        .run(NEG_INFINITY..INFINITY)
        .estimate()
        .unwrap();
    Some((alt / null).ln())
}

#[cfg(test)]
mod tests {

    use float_cmp::approx_eq;
    use crate::bayes::paired_d::bf_onesample_default_t;
    use crate::components::priors::Prior;
    const TOL: f64 = 0.00001;

    #[test]
    fn test1() {
        let alternative = "two.sided";
        let prior = Prior {
            family: "Cauchy".to_string(),
            params: vec![0., (2f64).sqrt() / 2.0],
        };
        let got = bf_onesample_default_t(
            0.508114037354055 / 1.06850369859537,
            140,
            &prior,
            alternative,
        )
        .unwrap();
        let want = (113731.66890204486844595522 as f64).ln();
        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );
    }

    #[test]
    fn test2() {
        let prior = Prior {
            family: "Cauchy".to_string(),
            params: vec![0., (2f64).sqrt() / 2.0],
        };
        let alternative = "two.sided";
        let want = (304475.15582726168213412166 as f64).ln();
        let got = bf_onesample_default_t(
            0.508114037354055 / 1.06850369859537,
            150,
            &prior,
            alternative,
        )
        .unwrap();
        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );
    }


}

