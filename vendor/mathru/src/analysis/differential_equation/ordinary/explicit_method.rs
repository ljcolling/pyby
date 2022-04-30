use super::ExplicitODE;
use crate::algebra::linear::Vector;

///
pub trait ExplicitMethod<T>
{
    ///
    fn do_step<F>(&self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ExplicitODE<T>;

    fn order(&self) -> u8;
}

pub trait ExplicitEmbeddedMethod<T>
{
    ///
    fn do_step<F>(&self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> (Vector<T>, Vector<T>)
        where F: ExplicitODE<T>;

    fn order(&self) -> (u8, u8);
}
