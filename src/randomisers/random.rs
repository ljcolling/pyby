use rand::prelude::*;
use rand::rngs::SmallRng;
use rand::distributions::Normal as NormalRand;
// use itertools::Itertools;




// Generate a bunch of random numbers
pub fn make_random(effsize: f64, n: i32, seed: Option<u64>) -> Vec<f64> {
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

/*
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
*/
