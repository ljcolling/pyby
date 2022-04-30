//! Solves an ordinary differential equation using the 4th order Tsitouras algorithm
use crate::algebra::linear::{Vector};
use crate::algebra::abstr::Real;
use super::{explicit_method::ExplicitEmbeddedMethod};
use super::ExplicitODE;
use crate::analysis::differential_equation::ordinary::butcher::ButcherAdaptiveStepSize;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};


/// Solves an ordinary differential equation using the 4th order Tsitouras algorithm
///
///<http://users.uoa.gr/~tsitourasc/RK54_new_v2.pdf>
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Tsitouras<T>
{
    butcher: ButcherAdaptiveStepSize<T>,
}

impl<T> Default for Tsitouras<T> where T: Real
{
    fn default() -> Tsitouras<T>
    {
        let a: Vec<T> = vec![T::from_f64(0.161),
                             T::from_f64(-0.008480655492356992), T::from_f64(0.3354806554923570),
                             T::from_f64(2.8971530571054944), T::from_f64(-6.359448489975075), T::from_f64(4.362295432869581),
                             T::from_f64(5.3258648284392596), T::from_f64(-11.74888356406283), T::from_f64(7.495539342889836), T::from_f64(-0.09249506636175525),
                             T::from_f64(5.8614554429464203), T::from_f64(-12.92096931784711), T::from_f64(8.159367898576159), T::from_f64(-0.07158497328140100), T::from_f64(-0.02826905039406838),
                             T::from_f64(0.09646076681806523), T::from_f64(0.01), T::from_f64(0.4798896504144996), T::from_f64(1.379008574103742), T::from_f64(-3.290069515436081), T::from_f64(2.324710524099774)];

        let b: Vec<T> = vec![T::from_f64(0.09646076681806523), T::from_f64(0.01), T::from_f64(0.4798896504144996), T::from_f64(1.379008574103742), T::from_f64(-3.290069515436081), T::from_f64(2.324710524099774), T::zero()];
        let b_s: Vec<T> = vec![T::from_f64(0.001780011052226), T::from_f64(0.000816434459657), T::from_f64(-0.007880878010262), T::from_f64(0.144711007173263), T::from_f64(-0.582357165452555), T::from_f64(0.458082105929187), T::from_f64(1.0 / 66.0)];
        let c: Vec<T> = vec![T::from_f64(0.161), T::from_f64(0.327), T::from_f64(0.9), T::from_f64(0.9800255409045097), T::one(), T::one()];

        Tsitouras {
            butcher: ButcherAdaptiveStepSize::new(a, b, b_s, c),
        }
    }
}

impl<T> ExplicitEmbeddedMethod<T> for Tsitouras<T>
    where T: Real
{

    fn do_step<F>(&self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> (Vector<T>, Vector<T>)
        where F: ExplicitODE<T>
    {
        self.butcher.do_step(prob, t_n, x_n, h)
    }

    fn order(&self) -> (u8, u8)
    {
        (5, 4)
    }
}