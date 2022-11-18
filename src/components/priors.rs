
// use crate::Prior;
use r_mathlib::cauchy_cdf as pcauchy;
use r_mathlib::cauchy_pdf as dcauchy;
use dict_derive::{FromPyObject, IntoPyObject};

// Definitions for priors and likelihoods

#[derive(FromPyObject, IntoPyObject, Debug, Clone)]
pub struct Prior {
    family: String,
    params: Vec<f64>,
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
/*
        let f = match self.family.as_str() {
            "Cauchy" => CreatePrior::new::<CauchyPrior>(self,alternative),
                _ => CreatePrior::new::<CauchyPrior>(self, alternative),
        }; */

        new(self, alternative)
        
    }
}

/*
pub trait CreatePrior {
    // fn new(prior: &Prior, alternative: &str) -> Self;
    fn new(prior: &Prior, alternative: &str) -> Box<dyn Fn(f64) -> f64>;
}
*/
// pub fn new<T: CreatePrior>(prior: &Prior, alternative: &str) -> T {
//     CreatePrior::new(prior, alternative)
// }

// pub fn new<T: CreatePrior>(prior: &Prior, alternative: &str) -> Box<dyn Fn(f64) -> f64> {
//     CreatePrior::new(prior, alternative)
// }

// impl dyn CreatePrior  {
    fn new(prior: &Prior, alternative: &str) -> Box<dyn Fn(f64) -> f64> {
        let location = prior.params[0];
        let scale = prior.params[1];
        let (ll, ul) =  alt_limits(alternative);
        let k = 1.0 / cauchy_auc(location, scale, ll, ul);
        return Box::new(move |x: f64| dcauchy(x, location, scale, false) * k);
    }
// }
// TODO: So far I've only implemented the cauchy prior


/* struct CauchyPrior {
    location: f64,
    scale: f64,
    ll: Option<f64>,
    ul: Option<f64>,
}
*/
// pub trait PriorFunction {
//     fn prior_function(&'static self) -> Box<dyn Fn(f64) -> f64>;
// }

// impl PriorFunction for CauchyPrior {
//     fn prior_function(&self) -> Box<dyn Fn(f64) -> f64> {
//         let k = 1.0 / cauchy_auc(self.location, self.scale, self.ll, self.ul);
//         let location = self.location;
//         let scale = self.scale;
//         return Box::new(move |x: f64| dcauchy(x, location, scale, false) * k);
//     }
// }






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

// pub fn cauchy_prior(
//     location: f64,
//     scale: f64,
//     ll: Option<f64>,
//     ul: Option<f64>,
// ) -> Box<impl Fn(f64) -> f64> {
// }

/*
    pub fn function(&self, alternative: &str) -> Box<impl Fn(f64) -> f64> {
        let (ll, ul) = match alternative {
            "two.sided" => (None, None),
            "greater" => (Some(0.), None),
            "less" => (None, Some(0.)),
            _ => (Some(0.), Some(0.)),
        };

        match self.family.as_str() {
            "Cauchy" => cauchy_prior(self.params[0], self.params[1], ll, ul),
            "student_t" => cauchy_prior(0., 0., Some(0.), Some(0.)),
            _ => cauchy_prior(0.0, 0.0, Some(0.0), Some(0.0)),
        }
    } */
