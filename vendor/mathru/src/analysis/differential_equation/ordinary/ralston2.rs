//! Solves an ODE using Ralston2's method.
use super::{explicit_method::ExplicitMethod, ExplicitODE};
use crate::{
    algebra::{abstr::Real, linear::Vector},
    analysis::differential_equation::ordinary::ButcherFixedStepSize
};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::clone::Clone;

/// Solves an ODE using Ralston's 2nd order method.
///
/// <https://en.wikipedia.org/wiki/List_of_Runge-Kutta_methods>
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Ralston2<T>
{
    butcher: ButcherFixedStepSize<T>
}

impl<T> Default for Ralston2<T> where T: Real
{
    /// Creates a Ralston2 instance
    fn default() -> Ralston2<T>
    {
        let a: Vec<T> = vec![T::from_f64(2.0/3.0)];
        let b: Vec<T> = vec![T::from_f64(1.0/4.0), T::from_f64(3.0/4.0)];
        let c: Vec<T> = vec![T::from_f64(2.0/3.0)];

        Ralston2 {
            butcher: ButcherFixedStepSize::new(a, b, c)
        }
    }
}

impl<T> ExplicitMethod<T> for Ralston2<T> where T: Real
{
    fn do_step<F>(&self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ExplicitODE<T>
    {
        self.butcher.do_step(prob, t_n, x_n, h)
    }

    // Ralston's method is a 2rd order method
    fn order(&self) -> u8
    {
        2
    }
}
