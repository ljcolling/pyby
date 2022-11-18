
use dict_derive::{FromPyObject, IntoPyObject};
use gkquad::single::algorithm::QAGS;
use gkquad::single::Integrator;
use itertools::Itertools;
use pyo3::prelude::*;
use r_mathlib::cauchy_cdf as pcauchy;
use r_mathlib::cauchy_pdf as dcauchy;
use r_mathlib::non_central_t_pdf;
use r_mathlib::students_t_cdf;
use rand::distributions::Normal as NormalRand;
use rand::prelude::*;
use rand::rngs::SmallRng;
// use rand::thread_rng;
use itertools::izip;
use rand::seq::SliceRandom;
use rayon::prelude::*;
use std::f64::{INFINITY, NEG_INFINITY};
use itertools::enumerate;

mod components;

use crate::components::priors::Prior;
use crate::components::priors::Function;
use crate::components::priors::cauchy_auc;

#[derive(FromPyObject, IntoPyObject, Debug)]
struct SamplingRule {
    n_min: i32,
    n_max: i32,
    step_size: i32,
}

impl SamplingRule  {
    fn unpack(&self) -> (i32, i32, i32) {
        (self.n_min, self.n_max, self.step_size)
    }
}



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
        _ => None,
    };

    d.unwrap().table()
}

impl SimTable {
    fn table (self) -> Vec<(usize, f64, i32, Option<f64>, f64, f64, f64)> {

    let mut id = Vec::new();
    let mut effsize = Vec::new();
    let mut n = Vec::new();
    let mut bf = Vec::new();
    let mut emp_es = Vec::new();
    let mut t = Vec::new();
    let mut pvalue = Vec::new();

    self.rows.into_iter().for_each(|x| {
        id.push(x.id);
        effsize.push(x.effsize);
        n.push(x.n);
        bf.push(x.bf);
        emp_es.push(x.emp_es);
        t.push(x.t);
        pvalue.push(x.pvalue);

    });

    izip!(id, effsize, n, bf, emp_es, t, pvalue).collect_vec()

    }
}



fn bf_sim_paired(
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
    let d: Vec<SimRow> =
        make_random(effsize, n_max * reps, seed)
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
                #[inline(always)]
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
                    SimRow{id, effsize, n, bf, emp_es, t, pvalue}
                    // (id, effsize, n, bf, emp_es, t, pvalue)
                },
            )
            .collect();

    SimTable{rows: d.into_iter().collect()}
}


#[derive(Debug)]
struct SimRow {
    id: usize,
    effsize: f64,
    n: i32,
    bf: Option<f64>,
    emp_es: f64,
    t: f64,
    pvalue: f64,
}

struct SimTable {
    rows: Vec<SimRow>
}

#[derive(Debug)]
struct PairedSim {
    id: Vec<usize>,
    effsize: Vec<f64>,
    n: Vec<i32>,
    bf: Vec<Option<f64>>,
    emp_es: Vec<f64>,
    t: Vec<f64>,
    pvalue: Vec<f64>,
}
impl PairedSim {
    fn new() -> PairedSim {
        PairedSim {
            id: Vec::new(),
            effsize: Vec::new(),
            n: Vec::new(),
            bf: Vec::new(),
            emp_es: Vec::new(),
            t: Vec::new(),
            pvalue: Vec::new(),
        }
    }
}
impl FromIterator<(usize, f64, i32, Option<f64>, f64, f64, f64)> for PairedSim {
    fn from_iter<I: IntoIterator<Item = (usize, f64, i32, Option<f64>, f64, f64, f64)>>(
        iter: I,
    ) -> Self {
        let mut c = PairedSim::new();

        for item in iter {
            c.id.push(item.0);
            c.effsize.push(item.1);
            c.n.push(item.2);
            c.bf.push(item.3);
            c.emp_es.push(item.4);
            c.t.push(item.5);
            c.pvalue.push(item.6);
        }

        c
    }
}


pub fn one_sample_t(t: f64, n: i32, _tail: i16) -> f64 {
    2f64 * students_t_cdf(t, (n - 1) as f64, false, false)
}

#[pyfunction]
#[pyo3(text_signature = "(effsize, n_min, n_max, step_size, reps, seed)")]
fn bf_sim_independent(
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
}


fn welfords(
    (i, count, mean, squared_distances): &mut (usize, f64, f64, f64),
    new_value: &f64,
) -> Option<(usize, i32, f64, f64, f64)> {
    // Use Welford's method to work out the running mean and sd
    *count += 1f64;
    let delta = new_value - *mean;
    *mean += delta / *count;
    let delta2 = new_value - *mean;
    *squared_distances += delta * delta2;
    let sd = (*squared_distances / (*count - 1f64)).sqrt();

    Some((
        *i,
        *count as i32,     // item
        *new_value as f64, // current value
        *mean,             // running mean
        sd,                // running sd
                           // bf,
    ))
}


#[inline(always)]
fn bf_onesample_default_t(d: f64, n: i32, prior: &Prior, alternative: &str) -> Option<f64> {
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


#[pymodule]
fn pyby(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bf_sim, m)?)?;
    m.add_function(wrap_pyfunction!(bf_sim_independent, m)?)?;
    m.add_function(wrap_pyfunction!(make_random2, m)?)?;
    Ok(())
}

#[pyfunction]
fn make_random2() -> Vec<f64> {
    let seed = Some(123);
    let mut rng = match seed {
        None => SmallRng::from_entropy(),
        Some(seed) => SmallRng::seed_from_u64(seed),
    };
    // StandardNormal.sample_iter(&mut rng).take(n).collect()
    let sizes = vec![1., 2., 3.];
    let mut sizes_ret = sizes.repeat(1000 / 3);
    sizes_ret.shuffle(&mut rng);
    sizes_ret
        .iter()
        .map(|s: &f64| {
            NormalRand::new(*s, 1.0)
                .sample_iter(&mut rng)
                .take(100 as usize)
                .collect_vec()
        })
        .flatten()
        .collect_vec()

}

// Generate a bunch of random numbers
#[inline(always)]
fn make_random(effsize: f64, n: i32, seed: Option<u64>) -> Vec<f64> {
    // If there's no seed then seed from entropy
    let mut rng = match seed {
        None => SmallRng::from_entropy(),
        Some(seed) => SmallRng::seed_from_u64(seed),
    };
    // StandardNormal.sample_iter(&mut rng).take(n).collect()
    NormalRand::new(effsize, 1.0)
        .sample_iter(&mut rng)
        .take(n as usize)
        .collect()
}

fn make_random1(effsize: f64, n: i32, reps: i32, seed: Option<u64>) -> Vec<f64> {
    // If there's no seed then seed from entropy
    let mut rng = match seed {
        None => SmallRng::from_entropy(),
        Some(seed) => SmallRng::seed_from_u64(seed),
    };
    // StandardNormal.sample_iter(&mut rng).take(n).collect()
    let sizes = vec![effsize].repeat(reps as usize);
    sizes
        .iter()
        .map(|s: &f64| {
            NormalRand::new(*s, 1.0)
                .sample_iter(&mut rng)
                .take(n as usize)
                .collect_vec()
        })
        .flatten()
        .collect_vec()
}

fn between_d(m1: f64, m2: f64, s1: f64, s2: f64, n: f64) -> f64 {
    let md_diff = m1 - m2;
    let sd_pooled = sd_pooled(s1, s2, n);
    let d = md_diff / sd_pooled;
    return d;
}

fn sd_pooled(s1: f64, s2: f64, n: f64) -> f64 {
    ((((n - 1.0) * s1.powi(2)) + ((n - 1.0) * s2.powi(2))) / (n + n - 2.0)).sqrt()
}

fn between_t(m1: f64, m2: f64, s1: f64, s2: f64, n: f64) -> (f64, f64, f64) {
    let md_diff = m1 - m2;
    let sd_pooled = sd_pooled(s1, s2, n);

    let var_pooled = sd_pooled.powi(2);
    let t = md_diff / ((var_pooled / n) + (var_pooled / n)).sqrt();
    let df = n + n - 2.0;
    return (t, df, 0.0);
}

/// here is the stuff that will get moved to a library

#[inline(always)]
fn noncentral_d_likelihood(d: f64, n: f64) -> impl Fn(f64) -> f64 {
    move |x: f64| dt(d * n.sqrt(), n - 1., n.sqrt() * x)
}

#[inline(always)]



#[inline(always)]
pub fn dt(x: f64, df: f64, ncp: f64) -> f64 {
    non_central_t_pdf(x, df, ncp, false)
}

#[inline(always)]
fn create_predictive(
    likelihood: impl Fn(f64) -> f64,
    prior: impl Fn(f64) -> f64,
) -> impl Fn(f64) -> f64 {
    move |x: f64| likelihood(x) * prior(x)
}

#[cfg(test)]
mod tests {

    use super::*;
    use float_cmp::approx_eq;

    const TOL: f64 = 0.00001;

    #[test]
    fn test0() {
        let got = 1.0;
        let want = 1.0;
        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );
    }

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

    #[test]
    fn test3() {
        let got = dt(0.1, 10., 2.);
        let want = 0.06434707;

        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );
    }

    #[test]
    fn cauchy_area() {
        let want = 1.;
        let got = cauchy_auc(0., 1., None, None);

        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );

        let want = 0.5;
        let got = cauchy_auc(0., 1., Some(0.), None);
        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );

        let want = 0.3788811;
        let got = cauchy_auc(0., 1., Some(0.4), None);
        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );

        let want = 0.6211189;
        let got = cauchy_auc(0., 1., Some(-0.4), None);
        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );

        let want = 0.3361116;
        let got = cauchy_auc(0., 0.707, None, Some(-0.4));
        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );

        let want = 0.3277769;
        let got = cauchy_auc(0.0, 0.707, Some(-0.4), Some(0.4));
        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );
        let want = 0.07613601;
        let got = cauchy_auc(0.0, 0.707, Some(0.2), Some(0.4));
        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );
    }
}

