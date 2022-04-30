//! Solves an ODE using the Bogacki Shampine algorithm.
use super::{explicit_method::ExplicitEmbeddedMethod, ExplicitODE, };
use crate::algebra::{abstr::Real, linear::Vector};
use std::clone::Clone;
use crate::analysis::differential_equation::ordinary::butcher::ButcherAdaptiveStepSize;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};


/// Solves an ODE using the Bogacki Shampine algorithm of 2nd order.
///
///<https://en.wikipedia.org/wiki/Bogacki-Shampine_method>
///
/// # Example
///
/// ```
/// use mathru::{
///     algebra::linear::{Matrix, Vector},
///     analysis::differential_equation::ordinary::{BogackiShampine23, ProportionalControl, ExplicitODE},
///     *,
/// };
///
/// // Define ODE
/// // $y^{'} = ay = f(x, y) $
/// // $y = C a e^{at}$
/// // $'y(t_{s}) = C a e^{at_s} => C = \frac{y(t_s)}{ae^{at_s}}$
/// pub struct ExplicitODEProblem
/// {
///     time_span: (f64, f64),
///     init_cond: Vector<f64>,
/// }
///
/// impl Default for ExplicitODEProblem
/// {
///     fn default() -> ExplicitODEProblem
///     {
///         ExplicitODEProblem { time_span: (0.0, 2.0),
///                              init_cond: vector![0.5] }
///     }
/// }
///
/// impl ExplicitODE<f64> for ExplicitODEProblem
/// {
///     fn func(&self, t: &f64, x: &Vector<f64>) -> Vector<f64>
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
/// let problem: ExplicitODEProblem = ExplicitODEProblem::default();
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
/// // Solve the ODE
/// let (t, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem, &BogackiShampine23::default()).unwrap();
/// ```
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct BogackiShampine23<T>
{
    butcher: ButcherAdaptiveStepSize<T>
}

impl<T> Default for BogackiShampine23<T> where T: Real
{
    fn default() -> BogackiShampine23<T>
    {
        let a: Vec<T> = vec![T::from_f64(0.5),
                             T::zero(), T::from_f64(3.0/4.0),
                             T::from_f64(2.0/9.0), T::from_f64(1.0/3.0), T::from_f64(4.0/9.0)];
        let b: Vec<T> = vec![T::from_f64(2.0/9.0), T::from_f64(1.0/3.0), T::from_f64(4.0/9.0), T::zero()];
        let b_s: Vec<T> = vec![T::from_f64(7.0/24.0), T::from_f64(1.0/4.0), T::from_f64(1.0/3.0), T::from_f64(1.0/8.0)];
        let c: Vec<T> = vec![T::from_f64(0.5), T::from_f64(3.0/4.0), T::one()];

        BogackiShampine23 {
            butcher: ButcherAdaptiveStepSize::new(a, b, b_s, c)
        }
    }
}

impl<T> ExplicitEmbeddedMethod<T> for BogackiShampine23<T> where T: Real
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
        (3, 2)
    }
}
