//! Solves an ODE using Ralston's 4th order method.
use super::{explicit_method::ExplicitMethod, ExplicitODE};
use crate::{
    algebra::{abstr::Real, linear::Vector},
    analysis::differential_equation::ordinary::ButcherFixedStepSize,
};
use std::clone::Clone;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Solves an ODE using Ralston's 4th order method.
///
/// <https://en.wikipedia.org/wiki/List_of_Runge-Kutta_methods>
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Ralston4<T>
{
    butcher: ButcherFixedStepSize<T>
}

impl<T> Default for Ralston4<T> where T: Real
{
    /// Creates a Ralston4 instance
    fn default() -> Ralston4<T>
    {
        const SQRT5: f64 = 2.236067977499789694091736687312762;
        let a: Vec<T> = vec![T::from_f64(0.4),
                             T::from_f64((-2889.0 + 1428.0 * SQRT5) / 1024.0), T::from_f64((3785.0 - 1620.0 * SQRT5) / 1024.0),
                             T::from_f64((-3365.0 + 2094.0 * SQRT5) / 6040.0), T::from_f64((-975.0 - 3046.0 * SQRT5) / 2552.0), T::from_f64((467040.0 + 203968.0 * SQRT5) / 240845.0)];
        let b: Vec<T> = vec![T::from_f64((263.0 + 24.0 * SQRT5) / 1812.0), T::from_f64((125.0 - 1000.0 * SQRT5) / 3828.0), T::from_f64(1024.0 * (3356.0 + 1623.0 * SQRT5) / 5924787.0), T::from_f64((30.0 - 4.0 * SQRT5) / 123.0)];
        let c: Vec<T> = vec![T::from_f64(0.4), T::from_f64((14.0 - 3.0 * SQRT5) / 16.0), T::one()];

        Ralston4 {
            butcher: ButcherFixedStepSize::new(a, b, c)
        }
    }
}

impl<T> ExplicitMethod<T> for Ralston4<T> where T: Real
{
    fn do_step<F>(&self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ExplicitODE<T>
    {
        self.butcher.do_step(prob, t_n, x_n, h)
    }

    // Ralston's method is a 4th order method
    fn order(&self) -> u8
    {
        4
    }
}
