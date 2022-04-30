//! Adaptive step size stepper

use super::{explicit_method::ExplicitEmbeddedMethod, ExplicitODE};
use crate::algebra::{abstr::Real, linear::Vector};
use std::default::Default;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::clone::Clone;



/// Proportional Control
///
/// Numerical Recipes
/// William H. Press, Saul Teukolsky, William T. Vetterling und Brian P. Flannery,
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Copy, Debug)]
pub struct ProportionalControl<T>
{
    /// Step size
    n_max: u32,
    h_0: T,
    fac: T,
    fac_min: T,
    fac_max: T,
    /// abs_tol: Absolute tolerance. This is the tolerance on local error
    /// estimates, not necessarily the global error. Defaults to 1e-6.
    abs_tol: T,
    /// reltol: Relative tolerance. This is the tolerance on local error
    /// estimates, not necessarily the global error. Defaults to 1e-3.
    rel_tol: T,
}

impl<T> Default for ProportionalControl<T> where T: Real
{
    fn default() -> ProportionalControl<T>
    {
        ProportionalControl::new(1000,
                                    T::from_f64(0.02),
                                    T::from_f64(0.8),
                                    T::from_f64(0.001),
                                    T::from_f64(3.0),
                                    T::from_f64(10e-6),
                                    T::from_f64(10e-3))
    }
}

impl<T> ProportionalControl<T> where T: Real
{
    /// Creates an instance with the given step size
    ///
    /// # Param
    ///
    /// *'fac_min':
    /// *'fac_max': 1.5 <= fac_max <= 5.0
    pub fn new(n_max: u32,
               h_0: T,
               fac: T,
               fac_min: T,
               fac_max: T,
               abs_tol: T,
               rel_tol: T)
               -> ProportionalControl<T>
    {
        ProportionalControl { n_max,
                                 h_0,
                                 fac,
                                 fac_min,
                                 fac_max,
                                 abs_tol,
                                 rel_tol }
    }

    /// Returns the absolute tolerance
    pub fn get_abs_tol(&self) -> &T
    {
        &self.abs_tol
    }

    /// Returns the relative tolerance
    pub fn get_rel_tol(&self) -> &T
    {
        &self.rel_tol
    }

    /// Sets the absolute tolerance
    ///
    /// # Parameters
    ///
    /// * 'abs_tol': abs_tol >= 0.0
    /// # Panics
    ///
    /// if 'abs_tol' < 0.0
    pub fn set_abs_tol(&mut self, abs_tol: T)
    {
        if abs_tol < T::zero()
        {
            panic!();
        }
        self.abs_tol = abs_tol;
    }

    /// Sets the relative tolerance
    ///
    /// # Parameters
    ///
    /// * 'rel_tol': rel_tol >= 0.0
    /// # Panics
    ///
    /// if 'rel_tol' < 0.0
    pub fn set_rel_tol(&mut self, rel_tol: T)
    {
        if rel_tol < T::zero()
        {
            panic!();
        }
        self.rel_tol = rel_tol;
    }

    /// Solves `func` using the 4th order Runge-Kutta-Fehlberg algorithm.
    ///
    /// # Arguments
    ///
    /// * 'func' is an explicit ordinary differential equation
    /// * 'init' is the initial value at the time 't_start'
    /// * 't_span' Time span t_span.0 = t_start, t_span.1 = t_stop
    ///
    /// # Return
    ///
    /// The solver returns a vector and a matrix, containing the times used in
    /// each step of the algorithm and the respectful values for that time.
    ///
    /// # Panic
    ///
    /// if t_span.0 > t_span.1
    pub fn solve<F, M>(&self,
                       prob: &F,
                       method: &M)
                       -> Result<(Vec<T>, Vec<Vector<T>>), &'static str>
        where F: ExplicitODE<T>,
              M: ExplicitEmbeddedMethod<T>
    {
        let t_span: (T, T) = prob.time_span();
        let t_start: T = t_span.0;
        let t_stop: T = t_span.1;
        if t_start > t_stop
        {
            panic!();
        }

        let order: (u8, u8) = method.order();
        let q: T = T::from_u8(order.0);
        let l: T = T::one() / (q + T::one());
        let mut x_n: Vector<T> = prob.init_cond();
        let mut t_n: T = t_start;
        let mut h: T = self.h_0;

        let mut t_vec: Vec<T> = vec![t_n];

        let mut res_vec: Vec<Vector<T>> = vec![x_n.clone()];

        let mut n: u32 = 0;

        while n < self.n_max && t_n < t_stop
        {
            h = h.min(t_stop - t_n);
            //
            let (x_n_new, x_ne): (Vector<T>, Vector<T>) = method.do_step(prob, &t_n, &x_n, &h);
            let err: T = self.calc_error(&x_n_new, &x_ne);

            if err <= T::one()
            {
                t_n += h;
                x_n = x_n_new;
                t_vec.push(t_n);
                res_vec.push(x_n.clone());
                n += 1;
            }

            if err != T::zero()
            {
                let mut s: T = self.fac * (T::one() / err).pow(l);

                if s < self.fac_min
                {
                    s = self.fac_min;
                }

                if s > self.fac_max
                {
                    s = self.fac_max;
                }

                h = s * h;
            }
        }
        if t_n < t_stop
        {
            return Err("Maximum number of iterations reached");
        }
        Ok((t_vec, res_vec))
    }

    fn calc_error(&self, y: &Vector<T>, y_hat: &Vector<T>) -> T
    {
        let (_m, n) = y.dim();

        let mut sum: T = T::zero();

        for i in 0..n
        {
            let y_i: T = y[i];
            let y_hat_i: T = y_hat[i];
            let sc_i: T = self.abs_tol + y_i.abs() * self.rel_tol;

            let k: T = (y_i - y_hat_i) / sc_i;
            sum += k * k;
        }

        (sum / T::from_f64(n as f64)).sqrt()
    }
}
