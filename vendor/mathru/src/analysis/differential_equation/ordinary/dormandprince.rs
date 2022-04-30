//! Solves an ODE using the 5th order Runge-Kutta-Dormand-Prince algorithm.
use super::{explicit_method::ExplicitEmbeddedMethod, ExplicitODE};
use crate::algebra::{
    abstr::{Real},
    linear::Vector,
};
use crate::analysis::differential_equation::ordinary::butcher::ButcherAdaptiveStepSize;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::clone::Clone;


/// Solves an ODE using the 4th order Runge-Kutta-Dormand-Prince algorithm.
///
///<https://en.wikipedia.org/wiki/Dormand-Prince_method>
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
///     analysis::differential_equation::ordinary::{DormandPrince, ProportionalControl, ExplicitODE},
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
/// let n_max: u32 = 500;
/// let abs_tol: f64 = 10e-8;
/// let rel_tol: f64 = 10e-6;
///
/// let solver: ProportionalControl<f64> = ProportionalControl::new(n_max, h_0, fac, fac_min, fac_max, abs_tol, rel_tol);
/// let problem: ExplicitODE1 = ExplicitODE1::default();
///
/// // Solve the ODE
/// let (t, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem, &DormandPrince::default()).unwrap();
///
/// # }
/// ```
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct DormandPrince<T>
{
    butcher: ButcherAdaptiveStepSize<T>,
}

impl<T> Default for DormandPrince<T> where T: Real
{
    fn default() -> DormandPrince<T>
    {
        let a: Vec<T> = vec![T::from_f64(1.0 / 5.0),
                            T::from_f64(3.0 / 40.0), T::from_f64(9.0 / 40.0),
                            T::from_f64(44.0 / 45.0), T::from_f64(-56.0 / 15.0), T::from_f64(32.0 / 9.0),
                            T::from_f64(19372.0 / 6561.0), T::from_f64(-25360.0 / 2187.0), T::from_f64(64448.0 / 6561.0), T::from_f64(-212.0 / 729.0),
                            T::from_f64(9017.0 / 3168.0), T::from_f64(-355.0 / 33.0), T::from_f64(46732.0 / 5247.0), T::from_f64(49.0 / 176.0), T::from_f64(-5103.0 / 18656.0),
                            T::from_f64(35.0 / 384.0), T::zero(), T::from_f64(500.0 / 1113.0), T::from_f64(125.0 / 192.0), T::from_f64(-2187.0 / 6784.0), T::from_f64(11.0 / 84.0)];
        let b: Vec<T> = vec![T::from_f64(35.0 / 384.0), T::zero(), T::from_f64(500.0 / 1113.0), T::from_f64(125.0 / 192.0), T::from_f64(-2187.0 / 6784.0), T::from_f64(11.0 / 84.0), T::zero()];
        let b_s: Vec<T> = vec![T::from_f64(5179.0 / 57600.0), T::zero(), T::from_f64(7571.0 / 16695.0), T::from_f64(393.0 / 640.0), T::from_f64(-92097.0 / 339200.0), T::from_f64(187.0 / 2100.0), T::from_f64(1.0 / 40.0)];
        let c: Vec<T> = vec![T::from_f64(1.0 / 5.0), T::from_f64(3.0 / 10.0), T::from_f64(4.0 / 5.0), T::from_f64(8.0 / 9.0), T::one(), T::one()];

        DormandPrince {
            butcher: ButcherAdaptiveStepSize::new(a, b, b_s, c),
        }
    }
}

impl<T> ExplicitEmbeddedMethod<T> for DormandPrince<T> where T: Real
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
