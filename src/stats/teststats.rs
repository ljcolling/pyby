use r_mathlib::students_t_cdf;

/*
fn between_t(m1: f64, m2: f64, s1: f64, s2: f64, n: f64) -> (f64, f64, f64) {
    let md_diff = m1 - m2;
    let sd_pooled = sd_pooled(s1, s2, n);

    let var_pooled = sd_pooled.powi(2);
    let t = md_diff / ((var_pooled / n) + (var_pooled / n)).sqrt();
    let df = n + n - 2.0;
    return (t, df, 0.0);
}

pub fn sd_pooled(s1: f64, s2: f64, n: f64) -> f64 {
    ((((n - 1.0) * s1.powi(2)) + ((n - 1.0) * s2.powi(2))) / (n + n - 2.0)).sqrt()
}  */

pub fn one_sample_t(t: f64, n: i32, _tail: i16) -> f64 {
    2f64 * students_t_cdf(t, (n - 1) as f64, false, false)
}


