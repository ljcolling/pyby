//! Solves an ODE using the 4th order Runge-Kutta algorithm.
use super::{explicit_method::ExplicitMethod, ExplicitODE};
use crate::{
    algebra::{abstr::Real, linear::Vector},
    analysis::differential_equation::ordinary::ButcherFixedStepSize,
};
use std::clone::Clone;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};


/// Solves an ODE using the 4th order Runge-Kutta algorithm.
///
///<https://en.wikipedia.org/wiki/Rung-Kutta_methods>
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct RungeKutta4<T>
{
    butcher: ButcherFixedStepSize<T>
}

impl<T> Default for RungeKutta4<T> where T: Real
{
    /// Creates a Kutta3 instance
    fn default() -> RungeKutta4<T>
    {
        let a: Vec<T> = vec![T::from_f64(0.5),
                             T::zero(), T::from_f64(0.5),
                             T::zero(), T::zero(), T::one()];
        let b: Vec<T> = vec![T::from_f64(1.0 / 6.0), T::from_f64(1.0 / 3.0), T::from_f64(1.0 / 3.0), T::from_f64(1.0 / 6.0)];
        let c: Vec<T> = vec![T::from_f64(0.5), T::from_f64(0.5), T::one()];

        RungeKutta4 {
            butcher: ButcherFixedStepSize::new(a, b, c)
        }
    }
}

impl<T> ExplicitMethod<T> for RungeKutta4<T> where T: Real
{
    fn do_step<F>(&self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ExplicitODE<T>
    {
        self.butcher.do_step(prob, t_n, x_n, h)
    }

    // Runge-Kutta 4 is a fourth order method
    fn order(&self) -> u8
    {
        4
    }
}
