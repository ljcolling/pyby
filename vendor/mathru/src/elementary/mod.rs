//! Elementary functions
//!
//! Fore more information:
//! <https://en.wikipedia.org/wiki/Elementary_function>

mod power;
mod exponential;
mod trigonometry;
mod hyperbolic;

pub use self::{
    exponential::Exponential, hyperbolic::Hyperbolic, power::Power, trigonometry::Trigonometry,
};
