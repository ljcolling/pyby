// #[global_allocator]
// static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;
// use gkquad::single::integral_with_config;
use gkquad::single::Integrator;
use gkquad::single::algorithm::QAGS;
// use gkquad::single::algorithm::QAGP;
// use gkquad::single::IntegrationConfig;
// use gkquad::Tolerance::AbsAndRel;
use itertools::Itertools;
use pyo3::prelude::*;
use rand::distributions::Normal as NormalRand;
use rand::prelude::*;
use rand::rngs::SmallRng;
// use std::f32::consts::SQRT_2;
use std::f64::{INFINITY, NEG_INFINITY};
// use mathru::special::gamma::gamma;
// use statrs::distribution::Continuous;
// use statrs::distribution::Normal;
// use statrs::distribution::StudentsT;
// use statrs::distribution::Univariate;
use r_mathlib::students_t_cdf;
use r_mathlib::non_central_t_pdf;
// use r_mathlib::students_t_pdf;
use r_mathlib::cauchy_pdf;
// use r_mathlib::students_t_quantile;
// use rand::FromEntropy;
// use itertools::collect_tuple;
use rayon::prelude::*;
// use pyo3::wrap_pyfunction;
// extern crate rgsl;

// use std::fs::File;
// use std::io::{self, BufRead};
// use std::path::Path;
/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: f64, b: f64) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// Generate data
#[pyfunction]
fn gen_data(
    effsize: f64,
    n: usize,
    max_n: usize,
    seed: Option<u64>,
) -> Vec<(usize, f64, i32, Option<f64>, f64, f64, f64)> {
    // generate the random numbers
    
    // let x = make_random(effsize, n, seed);

    make_random(effsize, n, seed)
    // let mut list0 = Vec::new();
    // 
    // if let Ok(lines) = read_lines("/Users/lc663/GitHub/pyby/1_file.dat") {
    //     for line in lines {
    //         if let Ok(item) = line {
    //             list0.push(item.parse::<f64>().unwrap());
    //         }
    //     }
    // }

    // if let Ok(lines) = read_lines("/Users/lc663/GitHub/pyby/2_file.dat") {
    //     for line in lines {
    //         if let Ok(item) = line {
    //             list0.push(item.parse::<f64>().unwrap());
    //         }
    //     }
    // }
    // let x = list0;
    // // group into runs
    // // let max_n: usize = 200;
    // // let mut data_grouped = Vec::new();
    // // for chunk in &x.into_iter().chunks(max_n) {
    // //     data_grouped.push(chunk.collect_vec())
    // //     // data_grouped.push(chunk.scan(0f64, |state, x| {
    // //     //     *state = *state + x;
    // //     //     Some(*state)
    // //     //     }).collect_vec());

    // // }
    // let data_grouped: Vec<&[f64]> = x.chunks(max_n).collect();
    .chunks(max_n)
    .collect::<Vec<&[f64]>>()
    // let data_grouped: Vec<Vec<f64>> = x.chunks(max_n).map(|x| x.to_vec()).collect();
    // process each run
    // let runs = n / max_n;
    // let mut iter = Vec::new();
    // for i in 0..runs {
    // iter.push(data_grouped[i].iter().scan(0f64, |state, x| {
    //     *state = *state + x;
    //     Some(*state)
    //     }).collect_vec());
    // };
    //

    // let results: Vec<(usize, i32, f64, f64, f64)> = data_grouped
    // data_grouped
        .into_par_iter()
        .enumerate()
        .map(|(i, group)| {
            // let mut k = 0f64;
            // let mut S = 0f64;
            group
                .iter()
                .scan(
                    (0f64, 0f64, 0f64),
            #[inline(always)]
                    |(count, mean, squared_distances), new_value| {
                        // Use Welford's method to work out the running mean and sd
                        // see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
                        *count += 1f64;
                        let delta = new_value - *mean;
                        *mean = *mean + (delta / *count);
                        let delta2 = new_value - *mean;
                        *squared_distances = *squared_distances + (delta * delta2);
                        let sd = (*squared_distances / (*count - 1f64)).sqrt();
                        // let bf = if *count <  2.0 {
                        //     Some(0f64)
                        // } else {
                        //     bayesfactor(*mean, sd, *count as i32)
                        // };

                        Some((
                            i + 1usize,        // group
                            *count as i32,     // item
                            *new_value as f64, // current value
                            *mean,             // running mean
                            sd,                // running sd
                                               // bf,
                        ))
                    },
                )
                .collect_vec()
        })
        .flatten()
        .filter(|(_, n, _, _, _)| n >= &2) // filter out initial result
        // .collect()
    // return results
    // results
        // .filter(|(_, n, _, _, _)| n % &10 == 0 && n != &10) // filter out initial result
        // .filter(|(_, n, _, _, _)| n != &1) // filter out initial result
        // .into_par_iter()
        .map(
            #[inline(always)]
            |(id, n, _, mean, sd)| {
            let bf = if n < 2 {
                Some(0f64)
            } else {
             bayesfactor(mean, sd, n)
            };

            let emp_es = mean / sd;
            let se = sd / (n as f64).sqrt();
            let t = mean / se;
            let pvalue = 2f64 * students_t_cdf(t, (n as f64) - 1., false, false);
            (id, effsize, n, bf, emp_es, t, pvalue)
        })
        .collect()
}

// fn increment_block(i: usize, block : &mut u32, max_n: usize) -> u32 {
//     if i % max_n == 0 {
//         *block = *block + 1u32;
//     };
//     return *block
// }
#[inline(always)]
fn bayesfactor(mean: f64, sd: f64, n: i32) -> Option<f64> {
    let (lb, ub) = tail_to_bounds(0.0);
    let d = mean / sd;
    // let t = mean / (sd / (n as f64).sqrt());
    // if t > 10f64 {
    //     return Some(1000.0);
    // }
    let likelihood = noncentral_d_likelihood(d, n as f64);

    let null = likelihood(0.0);
    // let name = String::from("student_t");
    let params = vec![0.0, 2f64.sqrt() / 2f64, 1.0, lb, ub];
    // let prior = t_prior_one_tailed(location, scale, df, lb, ub);
    let prior = make_prior("student_t", params);
    let h1 = create_predictive(likelihood, prior);
     // .max_iters(100)
        // .points(&[0.])
    // let mut config = IntegrationConfig::default();
    // config.limit = 1000000;
    // config.tolerance = AbsAndRel(1e-9, 1e-9);
    // let alt = integral_with_config(|x: f64| h1(x), lb..ub, config).estimate();
    let alt = Integrator::new(|x: f64| h1(x)).algorithm(QAGS::new())
            // .max_iters(100)
            .points(&[0.])
            .run(NEG_INFINITY..INFINITY)
            .estimate().unwrap();
            
    // match alt {
    //     Err(_) => return None,
    //     _ => return Some(alt.unwrap() / null ),
    // };
    Some((alt  / null).ln())
}

// /// A Python module implemented in Rust.
#[pymodule]
fn pyby(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(gen_data, m)?)?;

    Ok(())
}

// Generate a bunch of random numbers

#[inline(always)]
fn make_random(effsize: f64, n: usize, seed: Option<u64>) -> Vec<f64> {
    // If there's no seed then seed from entropy
    let mut rng = match seed {
        None => SmallRng::from_entropy(),
        Some(seed) => SmallRng::seed_from_u64(seed),
    };
    // StandardNormal.sample_iter(&mut rng).take(n).collect()
    NormalRand::new(effsize, 1.0)
        .sample_iter(&mut rng)
        .take(n)
        .collect()
}

/// here is the stuff that will get moved to a library

#[inline(always)]
fn noncentral_d_likelihood(d: f64, n: f64) -> impl Fn(f64) -> f64 {
    move |x: f64| dt(d * n.sqrt(), n - 1., n.sqrt() * x)
}

#[inline(always)]
fn t_prior(location: f64, scale: f64, df: f64, lb: f64, ub: f64) -> Box<dyn Fn(f64) -> f64> {
    let k = if lb == NEG_INFINITY && ub == INFINITY {
        1f64
    } else {
        range_area_studentt(location, scale, df, lb, ub)
    };

    // return Box::new(move |x: f64| students_t_pdf((x - location) / scale, df, false) * k);
    // let dist = StudentsT::new(location, scale, df).unwrap();
    // return Box::new(move |x: f64| dist.pdf(x) * k);
    return Box::new(move |x: f64 | cauchy_pdf(x, location, scale, false) * k)
}

  // stats::dt((x - mean) / sd, df, ncp = ncp, log = FALSE) / sd
#[inline(always)]
fn range_area_studentt(location: f64, scale: f64, df: f64, ll: f64, ul: f64) -> f64 {
    if location >= 0f64 {
        if ll == 0.0 && ul == INFINITY {
            return 1.0 / pt(location, scale, df, ll);
        } else if ll == NEG_INFINITY && ul == 0.0 {
            return 1.0 / (1.0 - pt(location, scale, df, 0.0));
        } else {
            return 1.0;
        }
    } else {
        if ll == 0.0 && ul == INFINITY {
            return 1.0 / (1.0 - pt(location, scale, df, ll));
        } else if ll == NEG_INFINITY && ul == 0.0 {
            return 1.0 / (pt(location, scale, df, 0.0));
        } else {
            return 1.0;
        }
    }
}

#[inline(always)]
fn pt(_location: f64, _scale: f64, df: f64, q: f64) -> f64 {
    // let dist = StudentsT::new(location, scale, df).unwrap();
    // if q >= location {
    //     return dist.cdf(q);
    // } else {
    //     return 1.0 - dist.cdf(q);
    // }
    return students_t_cdf(q, df, true , false) 

}

#[inline(always)]
fn tail_to_bounds(tail: f64) -> (f64, f64) {
    let (lb, ub) = if tail == 1.0 {
        (0f64, INFINITY)
    } else if tail == -1.0 {
        (NEG_INFINITY, 0f64)
    } else {
        (NEG_INFINITY, INFINITY)
    };

    return (lb, ub);
}

#[inline(always)]
pub fn dt(x: f64, df: f64, ncp: f64) -> f64 {
    non_central_t_pdf(x, df, ncp, false)
}



#[inline(always)]
fn make_prior(_name: &str, params: Vec<f64>) -> Box<impl Fn(f64) -> f64> {
    // if name == "normal" {
        // Box::new(normal_prior(params[0], params[1]))
    // } else {
        Box::new(t_prior(
            params[0], params[1], params[2], params[3], params[4],
        ))
    // }
}

#[inline(always)]
fn create_predictive(
    likelihood: impl Fn(f64) -> f64,
    prior: impl Fn(f64) -> f64,
) -> impl Fn(f64) -> f64 {
    move |x: f64| likelihood(x) * prior(x)
}

// #[inline(always)]
// fn normal_prior(mean: f64, sd: f64) -> Box<dyn Fn(f64) -> f64> {
//     let dist = Normal::new(mean, sd).unwrap();
//     return Box::new(move |x: f64| dist.pdf(x));
// }

//         // let mut w = rgsl::IntegrationWorkspace::new(10000000).expect("IntegrationWorkspace::new failed");
//         // let f = |x: f64| predictive(x);

//         // let (_, result2, _) = w.qagi(f,1e-10, 1e-10,10000000);
//         // let result2 = qagi(f, f64::NEG_INFINITY, f64::INFINITY, 0, 1e-3);

// fn bayesfactor2(mean: f64, sd: f64, n: i32) -> Option<f64> {
//     let (lb, ub) = tail_to_bounds(0.0);
//     let d = mean / sd;
//     // let t = mean / (sd / (n as f64).sqrt());
//     // if t > 10f64 {
//     //     return Some(1000.0);
//     // }
//     let likelihood = noncentral_d_likelihood(d, n as f64);

//     let null = likelihood(0.0);
//     // let name = String::from("student_t");
//     let params = vec![0.0, 2f64.sqrt() / 2f64, 1.0, lb, ub];
//     // let prior = t_prior_one_tailed(location, scale, df, lb, ub);
//     let prior = make_prior("student_t", params);
//     let h1 = create_predictive(likelihood, prior);

//     let mut w = rgsl::IntegrationWorkspace::new(1000000).expect("IntegrationWorkspace::new failed");
//     let f = |x: f64| h1(x);

//     let (_, alt, _) = w.qagi(f, 1e-9, 1e-9, 100000);

//     return Some(alt / null);
// }

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
        let got = bayesfactor( 0.508114037354055, 1.06850369859537 , 140 ).unwrap();
        let want = 113731.66890204486844595522;
        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );
    }


    #[test]
    fn test2() {
        let want = 304475.15582726168213412166;
        let got = bayesfactor(0.508114037354055, 1.06850369859537, 150).unwrap();
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
        let want  = 0.06434707;

        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );
    }
}


// fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
// where
//     P: AsRef<Path>,
// {
//     let file = File::open(filename)?;
//     Ok(io::BufReader::new(file).lines())
// }
