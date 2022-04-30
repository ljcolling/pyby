//! Solves an ODE using  Runge-Kutta-Fehlberg1(2) algorithm.
use super::{explicit_method::ExplicitEmbeddedMethod, ExplicitODE, };
use crate::algebra::{abstr::Real, linear::Vector};
use std::default::Default;
use std::clone::Clone;
use crate::analysis::differential_equation::ordinary::butcher::ButcherAdaptiveStepSize;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};


/// Solves an ODE using Runge-Kutta-Fehlberg1(2) algorithm.
///
///<https://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method>
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
///     analysis::differential_equation::ordinary::{ExplicitODE, ProportionalControl, Fehlberg12},
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
/// let h_0: f64 = 0.0001;
/// let fac: f64 = 0.9;
/// let fac_min: f64 = 0.01;
/// let fac_max: f64 = 2.0;
/// let n_max: u32 = 100;
/// let abs_tol: f64 = 10e-6;
/// let rel_tol: f64 = 10e-3;
///
/// let solver: ProportionalControl<f64> = ProportionalControl::new(n_max, h_0, fac, fac_min, fac_max, abs_tol, rel_tol);
///
/// let problem: ExplicitODE1 = ExplicitODE1::default();
///
/// // Solve the ODE
/// let (t, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem, &Fehlberg12::default()).unwrap();
///
/// # }
/// ```
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Fehlberg12<T>
{
    butcher: ButcherAdaptiveStepSize<T>
}


impl<T> Default for Fehlberg12<T> where T: Real
{
    fn default() -> Fehlberg12<T>
    {
        let a: Vec<T> = vec![T::from_f64(0.5),
                             T::from_f64(1.0 / 256.0), T::from_f64(255.0 / 256.0)];
        let b: Vec<T> = vec![T::from_f64(1.0 / 512.0), T::from_f64(255.0 / 256.0), T::from_f64(1.0 / 512.0)];
        let b_s: Vec<T> = vec![T::from_f64(1.0 / 256.0), T::from_f64(255.0 / 256.0), T::zero()];
        let c: Vec<T> = vec![T::from_f64(0.5), T::one()];

        Fehlberg12 {
            butcher: ButcherAdaptiveStepSize::new(a, b, b_s, c)
        }
    }
}

impl<T> ExplicitEmbeddedMethod<T> for Fehlberg12<T> where T: Real
{
    fn do_step<F>(&self,
                  prob: &F,
                  t_n: &T,
                  x_n: &Vector<T>,
                  h_n: &T)
                  -> (Vector<T>, Vector<T>)
        where F: ExplicitODE<T>
    {
        self.butcher.do_step(prob, t_n, x_n, h_n)
    }

    fn order(&self) -> (u8, u8)
    {
        (2, 1)
    }
}
