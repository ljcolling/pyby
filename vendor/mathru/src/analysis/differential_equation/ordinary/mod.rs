//! Ordinary differential equation
//!
//! Fore more information:<br>
//! <https://en.wikipedia.org/wiki/Ordinary_differential_equation>
//!
//! Because ODE higher order can always be reduced to a system of first order
//! ODEs,  the implemented algorithms only support to solve first order ODEs.
//!
//! ```math
//! \frac{dy}{dt}=f(t, y)
//! ```

mod explicit_euler;
mod heun2;
mod midpoint;
mod ralston2;
mod heun3;
mod ralston3;
mod ssprk3;
mod kutta3;
mod kutta38;
mod rungekutta4;
mod ralston4;

mod fehlberg12;
mod implicit_euler;
mod fehlberg45;
// mod tsitouras;
mod adamsbashforth;
mod adaptive_stepper;
pub mod explicit_method;
mod explicit_ode;
mod fixed_stepper;
mod implicit_method;
mod implicit_ode;
pub mod problem;
pub mod butcher;
mod bdf;
mod bogackishampine23;
mod cashkarp45;
mod dormandprince;

pub use explicit_euler::ExplicitEuler;
pub use bogackishampine23::BogackiShampine23;
pub use cashkarp45::CashKarp45;
pub use dormandprince::DormandPrince;
pub use heun2::Heun2;
pub use heun3::Heun3;
pub use kutta3::Kutta3;
pub use kutta38::Kutta38;
pub use midpoint::Midpoint;
pub use ralston2::Ralston2;
pub use ralston3::Ralston3;
pub use ssprk3::Ssprk3;
pub use rungekutta4::RungeKutta4;
pub use ralston4::Ralston4;
pub use fehlberg12::Fehlberg12;
pub use fehlberg45::Fehlberg45;
pub use bdf::BDF;
pub use implicit_euler::ImplicitEuler;
pub use adamsbashforth::AdamsBashforth;
// pub use tsitouras::Tsitouras;
pub use explicit_ode::ExplicitODE;
pub use implicit_ode::ImplicitODE;
pub use butcher::ButcherFixedStepSize;
pub use fixed_stepper::FixedStepper;
pub use adaptive_stepper::ProportionalControl;
