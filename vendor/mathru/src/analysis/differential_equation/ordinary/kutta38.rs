//! Solves an ODE using the 4th order 3/8-rule Runge-Kutta algorithm.
use super::{explicit_method::ExplicitMethod, ExplicitODE};
use crate::{
    algebra::{abstr::Real, linear::Vector},
    analysis::differential_equation::ordinary::ButcherFixedStepSize
};
use std::clone::Clone;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};


/// Solves an ODE using the 4th order 3/8-rule Runge-Kutta algorithm.
///
///<https://en.wikipedia.org/wiki/Rung-Kutta_methods>
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
///     analysis::differential_equation::ordinary::{ExplicitODE, FixedStepper, Kutta38},
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
///
/// let step_size: f64 = 0.001;
/// let solver: FixedStepper<f64> = FixedStepper::new(step_size);
///
/// let problem: ExplicitODE1 = ExplicitODE1::default();
///
/// // Solve the ODE
/// let (t, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem, &Kutta38::default()).unwrap();
///
/// # }
/// ```
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Kutta38<T>
{
    butcher: ButcherFixedStepSize<T>
}

impl<T> Default for Kutta38<T> where T: Real
{
    fn default() -> Kutta38<T>
    {
        let a: Vec<T> = vec![T::from_f64(1.0/3.0),
                             T::from_f64(-1.0/3.0) ,T::one(),
                             T::one(), -T::one(), T::one()];
        let b: Vec<T> = vec![T::from_f64(1.0/8.0), T::from_f64(3.0/8.0), T::from_f64(3.0/8.0), T::from_f64(1.0/8.0)];
        let c: Vec<T> = vec![T::from_f64(1.0/3.0), T::from_f64(2.0/3.0), T::one()];

        Kutta38 {
            butcher: ButcherFixedStepSize::new(a, b, c)
        }
    }
}

impl<T> ExplicitMethod<T> for Kutta38<T> where T: Real
{
    fn do_step<F>(&self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ExplicitODE<T>
    {
        self.butcher.do_step(prob, t_n, x_n, h)
    }

    /// Kutta's method is a 4th order method
    fn order(&self) -> u8
    {
        4
    }
}
