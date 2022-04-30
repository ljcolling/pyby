use rand::Rng;
use rand::distributions::{Sample, IndependentSample};
use statistics::*;
use distribution::{Univariate, Discrete, Distribution, Binomial};
use Result;

/// Implements the [Bernoulli](https://en.wikipedia.org/wiki/Bernoulli_distribution)
/// distribution which is a special case of the [Binomial](https://en.wikipedia.org/wiki/Binomial_distribution)
/// distribution where `n = 1` (referenced [Here](./struct.Binomial.html))
///
/// # Examples
///
/// ```
/// use statrs::distribution::{Bernoulli, Discrete};
/// use statrs::statistics::Mean;
///
/// let n = Bernoulli::new(0.5).unwrap();
/// assert_eq!(n.mean(), 0.5);
/// assert_eq!(n.pmf(0), 0.5);
/// assert_eq!(n.pmf(1), 0.5);
/// ```
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Bernoulli {
    b: Binomial,
}

impl Bernoulli {
    /// Constructs a new bernoulli distribution with
    /// the given `p` probability of success.
    ///
    /// # Errors
    ///
    /// Returns an error if `p` is `NaN`, less than `0.0`
    /// or greater than `1.0`
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Bernoulli;
    ///
    /// let mut result = Bernoulli::new(0.5);
    /// assert!(result.is_ok());
    ///
    /// result = Bernoulli::new(-0.5);
    /// assert!(result.is_err());
    /// ```
    pub fn new(p: f64) -> Result<Bernoulli> {
        Binomial::new(p, 1).map(|b| Bernoulli { b: b })
    }

    /// Returns the probability of success `p` of the
    /// bernoulli distribution.
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Bernoulli;
    ///
    /// let n = Bernoulli::new(0.5).unwrap();
    /// assert_eq!(n.p(), 0.5);
    /// ```
    pub fn p(&self) -> f64 {
        self.b.p()
    }

    /// Returns the number of trials `n` of the
    /// bernoulli distribution. Will always be `1.0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Bernoulli;
    ///
    /// let n = Bernoulli::new(0.5).unwrap();
    /// assert_eq!(n.n(), 1);
    /// ```
    pub fn n(&self) -> u64 {
        1
    }
}

impl Sample<f64> for Bernoulli {
    /// Generate a random sample from a bernoulli
    /// distribution using `r` as the source of randomness.
    /// Refer [here](#method.sample-1) for implementation details
    fn sample<R: Rng>(&mut self, r: &mut R) -> f64 {
        super::Distribution::sample(self, r)
    }
}

impl IndependentSample<f64> for Bernoulli {
    /// Generate a random independent sample from a bernoulli
    /// distribution using `r` as the source of randomness.
    /// Refer [here](#method.sample-1) for implementation details
    fn ind_sample<R: Rng>(&self, r: &mut R) -> f64 {
        super::Distribution::sample(self, r)
    }
}

impl Distribution<f64> for Bernoulli {
    /// Generate a random sample from the
    /// bernoulli distribution using `r` as the source
    /// of randomness where the generated
    /// values are `1` with probability `p` and `0`
    /// with probability `1-p`.
    ///
    /// # Examples
    ///
    /// ```
    /// # extern crate rand;
    /// # extern crate statrs;
    /// use rand::StdRng;
    /// use statrs::distribution::{Bernoulli, Distribution};
    ///
    /// # fn main() {
    /// let mut r = rand::StdRng::new().unwrap();
    /// let n = Bernoulli::new(0.5).unwrap();
    /// print!("{}", n.sample::<StdRng>(&mut r));
    /// # }
    /// ```
    fn sample<R: Rng>(&self, r: &mut R) -> f64 {
        self.b.sample(r)
    }
}

impl Univariate<u64, f64> for Bernoulli {
    /// Calculates the cumulative distribution
    /// function for the bernoulli distribution at `x`.
    ///
    /// # Panics
    ///
    /// If `x < 0.0` or `x > 1.0`
    ///
    /// # Formula
    ///
    /// ```ignore
    /// if x < 0 { 0 }
    /// else if x >= 1 { 1 }
    /// else { 1 - p }
    /// ```
    fn cdf(&self, x: f64) -> f64 {
        self.b.cdf(x)
    }
}

impl Min<u64> for Bernoulli {
    /// Returns the minimum value in the domain of the
    /// bernoulli distribution representable by a 64-
    /// bit integer
    ///
    /// # Formula
    ///
    /// ```ignore
    /// 0
    /// ```
    fn min(&self) -> u64 {
        0
    }
}

impl Max<u64> for Bernoulli {
    /// Returns the maximum value in the domain of the
    /// bernoulli distribution representable by a 64-
    /// bit integer
    ///
    /// # Formula
    ///
    /// ```ignore
    /// 1
    /// ```
    fn max(&self) -> u64 {
        1
    }
}

impl Mean<f64> for Bernoulli {
    /// Returns the mean of the bernoulli
    /// distribution
    ///
    /// # Formula
    ///
    /// ```ignore
    /// p
    /// ```
    fn mean(&self) -> f64 {
        self.b.mean()
    }
}

impl Variance<f64> for Bernoulli {
    /// Returns the variance of the bernoulli
    /// distribution
    ///
    /// # Formula
    ///
    /// ```ignore
    /// p * (1 - p)
    /// ```
    fn variance(&self) -> f64 {
        self.b.variance()
    }

    /// Returns the standard deviation of the bernoulli
    /// distribution
    ///
    /// # Formula
    ///
    /// ```ignore
    /// sqrt(p * (1 - p))
    /// ```
    fn std_dev(&self) -> f64 {
        self.b.std_dev()
    }
}

impl Entropy<f64> for Bernoulli {
    /// Returns the entropy of the bernoulli
    /// distribution
    ///
    /// # Formula
    ///
    /// ```ignore
    /// q = (1 - p)
    /// -q * ln(q) - p * ln(p)
    /// ```
    fn entropy(&self) -> f64 {
        self.b.entropy()
    }
}

impl Skewness<f64> for Bernoulli {
    /// Returns the skewness of the bernoulli
    /// distribution
    ///
    /// # Formula
    ///
    /// ```ignore
    /// q = (1 - p)
    /// (1 - 2p) / sqrt(p * q)
    /// ```
    fn skewness(&self) -> f64 {
        self.b.skewness()
    }
}

impl Median<f64> for Bernoulli {
    /// Returns the median of the bernoulli
    /// distribution
    ///
    /// # Formula
    ///
    /// ```ignore
    /// if p < 0.5 { 0 }
    /// else if p > 0.5 { 1 }
    /// else { 0.5 }
    /// ```
    fn median(&self) -> f64 {
        self.b.median()
    }
}

impl Mode<u64> for Bernoulli {
    /// Returns the mode of the bernoulli distribution
    ///
    /// # Formula
    ///
    /// ```ignore
    /// if p < 0.5 { 0 }
    /// else { 1 }
    /// ```
    fn mode(&self) -> u64 {
        self.b.mode()
    }
}

impl Discrete<u64, f64> for Bernoulli {
    /// Calculates the probability mass function for the
    /// bernoulli distribution at `x`.
    ///
    /// # Panics
    ///
    /// If `x > 1`
    ///
    /// # Formula
    ///
    /// ```ignore
    /// if x == 0 { 1 - p }
    /// else { p }
    /// ```
    fn pmf(&self, x: u64) -> f64 {
        self.b.pmf(x)
    }

    /// Calculates the log probability mass function for the
    /// bernoulli distribution at `x`.
    ///
    /// # Panics
    ///
    /// If `x > 1`
    ///
    /// # Formula
    ///
    /// ```ignore
    /// else if x == 0 { ln(1 - p) }
    /// else { ln(p) }
    /// ```
    fn ln_pmf(&self, x: u64) -> f64 {
        self.b.ln_pmf(x)
    }
}
