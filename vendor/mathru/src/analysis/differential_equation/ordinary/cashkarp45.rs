//! Solves an ODE using the 4th order Cash-Karp algorithm.
use std::clone::Clone;
use std::default::Default;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::algebra::{abstr::Real, linear::Vector};
use crate::analysis::differential_equation::ordinary::{
    explicit_method::ExplicitEmbeddedMethod,
    explicit_ode::ExplicitODE,
    butcher::ButcherAdaptiveStepSize,
};


/// Solves an ODE using the 4th order Cash-Karp algorithm.
///
/// <https://en.wikipedia.org/wiki/Cash-Karp_method>
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
///     analysis::differential_equation::ordinary::{CashKarp45, ProportionalControl, ExplicitODE},
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
/// // Solve the ODE with CashKarp45 algorithm
/// let (t, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem, &CashKarp45::default()).unwrap();
/// # }
/// ```
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct CashKarp45<T>
{
    butcher: ButcherAdaptiveStepSize<T>
}

impl<T> Default for CashKarp45<T> where T: Real
{
    fn default() -> CashKarp45<T>
    {
        let a: Vec<T> = vec![T::from_f64(1.0 / 5.0),
                             T::from_f64(3.0 / 40.0), T::from_f64(9.0 / 40.0),
                             T::from_f64(3.0 / 10.0), T::from_f64(-9.0 / 10.0), T::from_f64(6.0 / 5.0),
                             T::from_f64(-11.0 / 54.0), T::from_f64(5.0 / 2.0), T::from_f64(-70.0 / 27.0), T::from_f64(35.0 / 27.0),
                             T::from_f64(1631.0 / 55296.0), T::from_f64(175.0 / 512.0), T::from_f64(575.0 / 13824.0), T::from_f64(44275.0 / 110592.0), T::from_f64(253.0 / 4096.0)];
        let b_s: Vec<T> = vec![T::from_f64(37.0 / 378.0), T::zero(), T::from_f64(250.0 / 621.0), T::from_f64(125.0 / 594.0), T::zero(), T::from_f64(512.0 / 1771.0)];
        let b: Vec<T> = vec![T::from_f64(2825.0 / 27648.0), T::zero(), T::from_f64(18575.0 / 48384.0), T::from_f64(13525.0 / 55296.0), T::from_f64(277.0 / 14336.0), T::from_f64(1.0 / 4.0)];
        let c: Vec<T> = vec![T::from_f64(1.0 / 5.0), T::from_f64(3.0 / 10.0), T::from_f64(3.0 / 5.0), T::one(), T::from_f64(7.0 / 8.0)];

        CashKarp45 {
            butcher: ButcherAdaptiveStepSize::new(a, b, b_s, c)
        }
    }
}

impl<T> ExplicitEmbeddedMethod<T> for CashKarp45<T> where T: Real
{
    fn do_step<F>(&self,
                  prob: &F,
                  t_n: &T,
                  x_n: &Vector<T>,
                  h: &T)
                  -> (Vector<T>, Vector<T>)
        where F: ExplicitODE<T>
    {
        self.butcher.do_step(prob, t_n, x_n, h)
    }

    fn order(&self) -> (u8, u8)
    {
        (5, 4)
    }
}
