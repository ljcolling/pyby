//! Solves an ODE using the 4th order Runge-Kutta-Fehlberg algorithm.
use super::{explicit_method::ExplicitEmbeddedMethod, ExplicitODE,
};
use crate::algebra::{abstr::Real, linear::Vector};
use std::default::Default;
use std::clone::Clone;
use crate::analysis::differential_equation::ordinary::butcher::ButcherAdaptiveStepSize;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};


/// Solves an ODE using the 4th order Runge-Kutta-Fehlberg algorithm.
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
///     analysis::differential_equation::ordinary::{ExplicitODE, ProportionalControl, Fehlberg45},
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
/// let h_0: f64 = 0.1;
/// let fac: f64 = 0.9;
/// let fac_min: f64 = 0.01;
/// let fac_max: f64 = 2.0;
/// let n_max: u32 = 100;
/// let abs_tol: f64 = 10e-6;
/// let rel_tol: f64 = 10e-6;
///
/// let problem: ExplicitODE1 = ExplicitODE1::default();
/// let solver: ProportionalControl<f64> = ProportionalControl::new(n_max, h_0, fac, fac_min, fac_max, abs_tol, rel_tol);
///
/// // Solve the ODE
/// let (t, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem, &Fehlberg45::default()).unwrap();
///
/// # }
/// ```
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Fehlberg45<T>
{
    butcher: ButcherAdaptiveStepSize<T>,
}

impl<T> Default for Fehlberg45<T> where T: Real
{
    fn default() -> Fehlberg45<T>
    {
        let a: Vec<T> = vec![T::from_f64(1.0/4.0),
                             T::from_f64(3.0/32.0), T::from_f64(9.0/32.0),
                             T::from_f64(1932.0/2197.0), T::from_f64(-7200.0/2197.0),T::from_f64(7296.0/2197.0),
                             T::from_f64(439.0/216.0), T::from_f64(-8.0), T::from_f64(3680.0/513.0), T::from_f64(-845.0/4104.0),
                             T::from_f64(-8.0/27.0), T::from_f64(2.0), T::from_f64(-3544.0/2565.0), T::from_f64(1859.0/4104.0), T::from_f64(-11.0/40.0)];
        let b: Vec<T> = vec![T::from_f64(16.0/135.0), T::zero(), T::from_f64(6656.0/12825.0), T::from_f64(28561.0/56430.0), T::from_f64(-9.0/50.0), T::from_f64(2.0/55.0)];
        let b_s: Vec<T> = vec![T::from_f64(25.0/216.0), T::zero(), T::from_f64(1408.0/2565.0), T::from_f64(2197.0/4104.0), T::from_f64(-1.0/5.0), T::zero()];
        let c: Vec<T> = vec![T::from_f64(1.0 /4.0), T::from_f64(3.0/8.0), T::from_f64(12.0/13.0), T::one(), T::from_f64(0.5)];

        Fehlberg45 {
            butcher: ButcherAdaptiveStepSize::new(a, b, b_s, c),
        }
    }
}

impl<T> ExplicitEmbeddedMethod<T> for Fehlberg45<T> where T: Real
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
        (5, 4)
    }
}
