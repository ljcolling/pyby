
pub fn create_predictive(
    likelihood: impl Fn(f64) -> f64,
    prior: impl Fn(f64) -> f64,
) -> impl Fn(f64) -> f64 {
    move |x: f64| likelihood(x) * prior(x)
}
