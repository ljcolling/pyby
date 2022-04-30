//! Normal distribution
use crate::{
    algebra::abstr::Real,
    special::error,
    statistics::distrib::{Continuous, Distribution},
    special::gamma::Gamma,
    special::error::Error,
};
use rand;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::clone::Clone;

/// Normal distribution
///
/// Fore more information:
/// <https://en.wikipedia.org/wiki/Normal_distribution>
///
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Copy, Debug)]
pub struct Normal<T>
{
    mean: T,
    variance: T,
}

impl<T> Normal<T> where T: Real
{
    /// Creates a probability distribution
    ///
    /// # Arguments
    ///
    /// * `mean`: Expected value
    /// * `variance`: variance > 0.0
    ///
    /// # Panics
    ///
    /// if variance <= 0.0
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::Normal;
    ///
    /// let distrib: Normal<f64> = Normal::new(0.3, 0.2);
    /// ```
    pub fn new(mean: T, variance: T) -> Self
    {
        if variance <= T::zero()
        {
            panic!();
        }

        Normal { mean, variance }
    }

    /// It is assumed that data are normal distributed.
    ///
    /// data.len() >= 2
    pub fn from_data(data: &Vec<T>) -> Self
    {
        let n: usize = data.len();
        if n < 2
        {
            panic!()
        }

        let mean: T = Normal::calc_mean(data);
        let variance: T = Normal::calc_variance(data, mean);

        Normal::new(mean, variance)
    }

    fn calc_mean(data: &Vec<T>) -> T
    {
        let n: usize = data.len();
        let mut sum: T = T::zero();

        for x in data.iter()
        {
            sum += *x;
        }

        sum / T::from_u64(n as u64)
    }

    fn calc_variance(data: &Vec<T>, mean: T) -> T
    {
        let n: usize = data.len();
        let mut sum: T = T::zero();

        for x in data.iter()
        {
            sum += (*x - mean).pow(T::from_f64(2.0));
        }

        sum / T::from_u64((n - 1) as u64)
    }
}

impl<T> Continuous<T> for Normal<T>
    where T: Real + Gamma + Error
{
    /// Probability density function
    ///
    /// # Arguments
    ///
    /// * `x`:  x &isin; &#x2115
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, Normal};
    ///
    /// let distrib: Normal<f64> = Normal::new(0.3, 0.2);
    /// let x: f64 = 5.0;
    /// let p: f64 = distrib.pdf(x);
    /// ```
    fn pdf(&self, x: T) -> T
    {
        let z: T = T::from_f64(-0.5) * ((x - self.mean) / self.variance).pow(T::from_f64(2.0));
        let f: T = T::one() / (self.variance * T::from_f64(2.0) * T::pi()).sqrt();

        f * z.exp()
    }

    /// Cumulative distribution function
    ///
    /// # Arguments
    ///
    /// * `x`:
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Continuous, Normal};
    ///
    /// let distrib: Normal<f64> = Normal::new(0.3, 0.2);
    /// let x: f64 = 0.4;
    /// let p: f64 = distrib.cdf(x);
    /// ```
    fn cdf(&self, x: T) -> T
    {
        let k: T = (x - self.mean) / ((T::from_f64(2.0) * self.variance).sqrt());
        let prob: T = T::from_f64(0.5) * (T::one() + error::erf(k));
        prob
    }

    /// Quantile: function of inverse cdf
    ///
    /// The Percentage Points of the Normal Distribution
    ///  Author(s): Michael J. Wichura
    /// Year 1988
    /// Journal of the Royal Statistical Society
    /// 0.0 < p < 1.0
    ///
    /// # Panics
    ///
    /// if  p <= 0.0 || p >= 1.0
    fn quantile(&self, p: T) -> T
    {
        if p <= T::zero() || p >= T::one()
        {
            panic!();
        }

        let mut ppnd16: T;
        let mut r: T;
        let q: T = p - T::from_f64(0.5);
        if q.abs() <= T::from_f64(0.425)
        {
            let r: T = T::from_f64(0.180625) - q * q;
            ppnd16 = q * Normal::r1(r);
        }
        else
        {
            if q < T::zero()
            {
                r = p
            }
            else
            {
                r = T::one() - p
            }
            if r <= T::zero()
            {
                ppnd16 = T::zero();
                return ppnd16;
            }
            r = (-r.ln()).sqrt();
            if r <= T::from_f64(5.)
            {
                r -= T::from_f64(1.6);
                ppnd16 = Normal::r2(r);
            }
            else
            {
                r -= T::from_f64(5.0);
                ppnd16 = Normal::r3(r);
            }
            if q <= T::zero()
            {
                ppnd16 = -ppnd16;
            }
        }

        self.mean + self.variance.sqrt() * ppnd16
    }

    /// Expected value
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{
    ///     self,
    ///     statistics::distrib::{Continuous, Normal},
    /// };
    ///
    /// let distrib: Normal<f64> = Normal::new(0.0, 0.2);
    /// let mean: f64 = distrib.mean();
    /// ```
    fn mean(&self) -> T
    {
        self.mean
    }

    /// Variance
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{
    ///     self,
    ///     statistics::distrib::{Continuous, Normal},
    /// };
    ///
    /// let distrib: Normal<f64> = Normal::new(0.0, 0.2);
    /// let var: f64 = distrib.variance();
    /// ```
    fn variance(&self) -> T
    {
        self.variance
    }

    /// Skewness
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{
    ///     self,
    ///     statistics::distrib::{Continuous, Normal},
    /// };
    /// let mean: f64 = 1.0;
    /// let variance: f64 = 0.5;
    /// let distrib: Normal<f64> = Normal::new(mean, variance);
    /// assert_eq!(0.0, distrib.skewness());
    /// ```
    fn skewness(&self) -> T
    {
        T::zero()
    }

    /// Median
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{
    ///     self,
    ///     statistics::distrib::{Continuous, Normal},
    /// };
    ///
    /// let mean: f64 = 0.0;
    ///
    /// let distrib: Normal<f64> = Normal::new(mean, 0.2);
    /// let median: f64 = distrib.median();
    /// ```
    fn median(&self) -> T
    {
        self.mean
    }

    /// Entropy
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{
    ///     self,
    ///     statistics::distrib::{Continuous, Normal},
    /// };
    /// use std::f64::consts::{E, PI};
    ///
    /// let mean: f64 = 1.0;
    /// let variance: f64 = 0.5;
    /// let distrib: Normal<f64> = Normal::new(mean, variance);
    ///
    /// let entropy: f64 =  distrib.entropy();
    /// ```
    fn entropy(&self) -> T
    {
        T::from_f64(2.0) * T::pi() * T::e() * self.variance
    }
}

impl<T> Distribution<T> for Normal<T> where T: Real
{
    ///
    ///  See Knuth The Art of Computer Programming Vol 2 3.4.1 C Algorithm P
    fn random(&self) -> T
    {
        let mut s: T = T::one();
        let mut v1: T = T::one();
        let mut v2: T;

        while s >= T::one()
        {
            let u1: T = T::from_f64(rand::random::<f64>());
            let u2: T = T::from_f64(rand::random::<f64>());
            v1 = T::from_f64(2.0) * u1 - T::one();
            v2 = T::from_f64(2.0) * u2 - T::one();
            s = v1 * v1 + v2 * v2
        }
        let x1: T = v1 * (-T::from_f64(2.0) * s.ln() / s).sqrt();
        x1 * self.variance.sqrt() + self.mean
    }
}

impl<T> Normal<T> where T: Real
{
    fn r1(t: T) -> T
    {
        let r: f64 = t.to_f64();
        let value: f64 = (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r
                               + 67265.770927008700853)
                              * r
                              + 45921.953931549871457)
                             * r
                             + 13731.693765509461125)
                            * r
                            + 1971.5909503065514427)
                           * r
                           + 133.14166789178437745)
                          * r
                          + 3.387132872796366608)
                         / (((((((r * 5226.495278852854561 + 28729.085735721942674) * r
                                 + 39307.89580009271061)
                                * r
                                + 21213.794301586595867)
                               * r
                               + 5394.1960214247511077)
                              * r
                              + 687.1870074920579083)
                             * r
                             + 42.313330701600911252)
                            * r
                            + 1.);
        T::from_f64(value)
    }

    fn r2(t: T) -> T
    {
        let r: f64 = t.to_f64();
        let value: f64 = (((((((r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r
                               + 0.24178072517745061177)
                              * r
                              + 1.27045825245236838258)
                             * r
                             + 3.64784832476320460504)
                            * r
                            + 5.7694972214606914055)
                           * r
                           + 4.6303378461565452959)
                          * r
                          + 1.42343711074968357734)
                         / (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4)
                                 * r
                                 + 0.0151986665636164571966)
                                * r
                                + 0.14810397642748007459)
                               * r
                               + 0.68976733498510000455)
                              * r
                              + 1.6763848301838038494)
                             * r
                             + 2.05319162663775882187)
                            * r
                            + 1.);

        T::from_f64(value)
    }

    fn r3(t: T) -> T
    {
        let r: f64 = t.to_f64();
        let value: f64 =
            (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r
                  + 0.0012426609473880784386)
                 * r
                 + 0.026532189526576123093)
                * r
                + 0.29656057182850489123)
               * r
               + 1.7848265399172913358)
              * r
              + 5.4637849111641143699)
             * r
             + 6.6579046435011037772)
            / (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r
                    + 1.846_318_317_510_054_681_8e-5)
                   * r
                   + 7.868_691_311_456_132_591e-4)
                  * r
                  + 0.014_875_361_290_850_614_8525)
                 * r
                 + 0.136_929_880_922_735_805_31)
                * r
                + 0.599_832_206_555_887_937_69)
               * r
               + 1.);

        T::from_f64(value)
    }
}
