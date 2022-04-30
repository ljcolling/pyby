//! Solves an ODE using Heun's method.
use super::{explicit_method::ExplicitMethod, ExplicitODE};
use crate::{
    algebra::{abstr::Real, linear::Vector},
    analysis::differential_equation::ordinary::ButcherFixedStepSize
};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::clone::Clone;


/// Solves an ODE using Heun's 3rd order method.
///
/// <https://en.wikipedia.org/wiki/Heun's_method>
///
/// # Example
///
/// For this example, we want to solve the following ordinary differential
/// equation:
/// ```math
/// \frac{dy}{dt} = ay = f(t, y)
/// ```
/// The initial condition is $y(0) = 0.5$ and we solve it in the interval
/// $\lbrack 0, 2\rbrack$ The following equation is the closed solution for
/// this ODE:
/// ```math
/// y(t) = C a e^{at}
/// ```
/// $C$ is a parameter and depends on the initial condition $y(t_{0})$
/// ```math
/// C = \frac{y(t_{0})}{ae^{at_{0}}}
/// ```
///
/// In this example, we set $a=2$
/// ```
/// # #[macro_use]
/// # extern crate mathru;
/// # fn main()
/// # {
/// use mathru::{
///     algebra::linear::Vector,
///     analysis::differential_equation::ordinary::{ExplicitODE, FixedStepper, Heun3},
/// };
///
/// pub struct ExplicitODE1
/// {
///     time_span: (f64, f64),
///     init_cond: Vector<f64>,
/// }
///
/// impl Default for ExplicitODE1
/// {
///     fn default() -> ExplicitODE1
///     {
///         ExplicitODE1 { time_span: (0.0, 2.0),
///                        init_cond: vector![0.5] }
///     }
/// }
///
/// impl ExplicitODE<f64> for ExplicitODE1
/// {
///     fn func(&self, _t: &f64, x: &Vector<f64>) -> Vector<f64>
///     {
///         return x * &2.0f64;
///     }
///
///     fn time_span(&self) -> (f64, f64)
///     {
///         return self.time_span;
///     }
///
///     fn init_cond(&self) -> Vector<f64>
///     {
///         return self.init_cond.clone();
///     }
/// }
///
/// // We instantiate Heun's algorithm with a step size of 0.001
/// let solver: FixedStepper<f64> = FixedStepper::new(0.001);
///
/// let problem: ExplicitODE1 = ExplicitODE1::default();
///
/// // Solve the ODE
/// let (t, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem, &Heun3::default()).unwrap();
///
/// # }
/// ```
///
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Heun3<T>
{
    butcher: ButcherFixedStepSize<T>
}

impl<T> Default for Heun3<T> where T: Real
{
    /// Creates a Heun3 instance
    fn default() -> Heun3<T>
    {
        let a: Vec<T> = vec![T::from_f64(1.0 / 3.0),
                             T::zero(), T::from_f64(2.0 / 3.0)];
        let b: Vec<T> = vec![T::from_f64(1.0 / 4.0), T::zero(), T::from_f64(3.0 / 4.0)];
        let c: Vec<T> = vec![T::from_f64(1.0 / 3.0), T::from_f64(2.0 / 3.0)];

        Heun3 {
            butcher: ButcherFixedStepSize::new(a, b, c)
        }
    }
}

impl<'a, T> ExplicitMethod<T> for Heun3<T> where T: Real
{
    fn do_step<F>(&self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ExplicitODE<T>
    {
        self.butcher.do_step(prob, t_n, x_n, h)
    }

    fn order(&self) -> u8
    {
        3
    }
}
