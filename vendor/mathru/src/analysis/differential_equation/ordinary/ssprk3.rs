//! Solves an ODE using Ssprk3's method.
use super::{explicit_method::ExplicitMethod, ExplicitODE};
use crate::{
    algebra::{abstr::Real, linear::Vector},
    analysis::differential_equation::ordinary::ButcherFixedStepSize,
};
use std::clone::Clone;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Solves an ODE using Third-order Strong Stability Preserving Runge-Kutta(SSPRK3).
///
/// <https://en.wikipedia.org/wiki/List_of_Runge-Kutta_methods>
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Ssprk3<T>
{
    butcher: ButcherFixedStepSize<T>,
}

impl<T> Default for Ssprk3<T> where T: Real
{
    /// Creates a Ssprk3 instance
    fn default() -> Ssprk3<T>
    {
        let a: Vec<T> = vec![T::one(),
                             T::from_f64(1.0 / 4.0), T::from_f64(1.0 / 4.0)];
        let b: Vec<T> = vec![T::from_f64(1.0 / 6.0), T::from_f64(1.0 / 6.0), T::from_f64(2.0 / 3.0)];
        let c: Vec<T> = vec![T::one(), T::from_f64(0.5)];

        Ssprk3 {
            butcher: ButcherFixedStepSize::new(a, b, c)
        }
    }
}

impl<T> ExplicitMethod<T> for Ssprk3<T> where T: Real
{
    fn do_step<F>(&self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ExplicitODE<T>
    {
        self.butcher.do_step(prob, t_n, x_n, h)
    }

    // Ssprk3's method is a 3rd order method
    fn order(&self) -> u8
    {
        3
    }
}
