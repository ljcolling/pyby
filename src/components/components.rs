use dict_derive::{FromPyObject, IntoPyObject};

#[derive(FromPyObject, IntoPyObject, Debug)]
struct Likelihood {
    family: String,
    params: Vec<f64>,
}
