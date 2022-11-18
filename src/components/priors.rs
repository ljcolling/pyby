use dict_derive::{FromPyObject, IntoPyObject};
use r_mathlib::cauchy_cdf as pcauchy;
use r_mathlib::cauchy_pdf as dcauchy;

// Definitions for priors
#[derive(FromPyObject, IntoPyObject, Debug, Clone)]
pub struct Prior {
    pub family: String,
    pub params: Vec<f64>,
}

fn alt_limits(alternative: &str) -> (Option<f64>, Option<f64>) {
    let (ll, ul) = match alternative {
        "two.sided" => (None, None),
        "greater" => (Some(0.), None),
        "less" => (None, Some(0.)),
        _ => (Some(0.), Some(0.)),
    };
    (ll, ul)
}

pub trait Function {
    fn function(&self, alternative: &str) -> Box<dyn Fn(f64) -> f64>;
}

impl Function for Prior {
    fn function(&self, alternative: &str) -> Box<dyn Fn(f64) -> f64> {
        new(self, alternative)
    }
}

// FIXME ME: This is just temporary while I finalise the API design
fn new(prior: &Prior, alternative: &str) -> Box<dyn Fn(f64) -> f64> {
    let location = prior.params[0];
    let scale = prior.params[1];
    let (ll, ul) = alt_limits(alternative);
    let k = 1.0 / cauchy_auc(location, scale, ll, ul);
    return Box::new(move |x: f64| dcauchy(x, location, scale, false) * k);
}

pub fn cauchy_auc(location: f64, scale: f64, ll: Option<f64>, ul: Option<f64>) -> f64 {
    match (ll, ul) {
        (None, None) => return 1.,
        (Some(_), None) => return pcauchy(ll.unwrap(), location, scale, false, false),
        (None, Some(_)) => return pcauchy(ul.unwrap(), location, scale, true, false),
        (Some(_), Some(_)) => {
            1.0 - ((1.0 - pcauchy(ll.unwrap(), location, scale, false, false))
                + (1.0 - pcauchy(ul.unwrap(), location, scale, true, false)))
        }
    }
}


#[cfg(test)]
mod tests {

    use super::*;
    use float_cmp::approx_eq;
    const TOL: f64 = 0.00001;

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
