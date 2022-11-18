use dict_derive::{FromPyObject, IntoPyObject};
use r_mathlib::non_central_t_pdf;

#[derive(FromPyObject, IntoPyObject, Debug)]
struct Likelihood {
    family: String,
    params: Vec<f64>,
}


// FIXME: Still haven't finalised API Design

pub fn noncentral_d_likelihood(d: f64, n: f64) -> impl Fn(f64) -> f64 {
    move |x: f64| dt(d * n.sqrt(), n - 1., n.sqrt() * x)
}

fn dt(x: f64, df: f64, ncp: f64) -> f64 {
    non_central_t_pdf(x, df, ncp, false)
}


#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::approx_eq;
    const TOL: f64 = 0.00001;

    #[test]
    fn dt_test() {
        let got = dt(0.1, 10., 2.);
        let want = 0.06434707;

        assert!(
            approx_eq!(f64, got, want, epsilon = TOL),
            "got = {}; want = {}",
            got,
            want
        );
    }
}
