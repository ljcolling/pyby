v0.6.0
 - `gamma::gamma_ur`, `gamma::gamma_ui`, `gamma::gamma_lr`, and `gamma::gamma_li` now follow strict gamma function domain, panicking if `a` or `x` are not in `(0, +inf)`
 - `beta::beta_reg` no longer allows `0.0` for `a` or `b` arguments
 - `InverseGamma` distribution no longer accepts `f64::INFINITY` as valid arguments for `shape` or `rate` as the value is nonsense
 - `Binomial::cdf` no longer accepts arguments outside the domain of `[0, n]`
 - `Bernoulli::cdf` no longer accepts arguments outside the domain of `[0, 1]`
 - `DiscreteUniform::cdf` no longer accepts arguments outside the domain of `[min, max]`
 - `Uniform::cdf` no longer accepts arguments outside the domain of `[min, max]`
 - `Triangular::cdf` no longer accepts arguments outside the domain of `[min, max]`
 - `FisherSnedecor` no longer accepts `f64::INFINITY` as a valid argument for `freedom_1` or `freedom_2`
 - `FisherSnedecor::cdf` no longer accepts arguments outside the domain of `[0, +inf)`
 - `Geometric::cdf` no longer accepts non-positive arguments
 - `Normal` now uses the Ziggurat method to generate random samples. This also affects all distributions depending on `Normal` for sampling
    including `Chi`, `LogNormal`, `Gamma`, and `StudentsT`
 - `Exponential` now uses the Ziggurat methd to generate random samples.
 - `Binomial` now implements `Univariate<u64, f64>` rather than `Univariate<i64, f64>`, meaning `Binomial::min` and `Binomial::max` now return `u64`
 - `Bernoulli` now implements `Univariate<u64, f64>` rather than `Univariate<i64, f64>`, meaning `Bernoulli::min` and `Bernoulli::min` now return `u64`
 - `Geometric` now implements `Univariate<u64, f64>` rather than `Univariate<i64, f64>`, meaning `Geometric::min` and `Geometric::min` now return `u64`
 - `Poisson` now implements `Univariate<u64, f64>` rather than `Univariate<i64, f64>`, meaning `Poisson::min` and `Poisson::min` now return `u64`
 - `Binomial` now implements `Mode<u64>` instead of `Mode<i64>`
 - `Bernoulli` now implements `Mode<u64>` instead of `Mode<i64>`
 - `Poisson` now implements `Mode<u64>` instead of `Mode<i64>`
 - `Geometric` now implements `Mode<u64>` instead of `Mode<i64>`
 - `Hypergeometric` now implements `Mode<u64>` instead of `Mode<i64>`
 - `Binomial` now implements `Discrete<u64, f64>` rather than `Discrete<i64, f64>`
 - `Bernoulli` now implements `Discrete<u64, f64>` rather than `Discrete<i64, f64>`
 - `Geometric` now implements `Discrete<u64, f64>` rather than `Discrete<i64, f64>`
 - `Hypergeometric` now implements `Discrete<u64, f64>` rather than `Discrete<i64, f64>`
 - `Poisson` now implements `Discrete<u64, f64>` rather than `Discrete<i64, f64>`

v0.5.1
 - Fixed critical bug in `normal::sample_unchecked` where it was returning `NaN`

v0.5.0
 - Implemented the `logistic::logistic` special function
 - Implemented the `logistic::logit` special function
 - Implemented the `factorial::multinomial` special function
 - Implemented the `harmonic::harmonic` special function
 - Implemented the `harmonic::gen_harmonic` special function
 - Implemented the `InverseGamma` distribution
 - Implemented the `Geometric` distribution
 - Implemented the `Hypergeometric ` distribution
 - `gamma::gamma_ur` now panics when `x > 0` or `a == f64::NEG_INFINITY`. In addition, it also returns `f64::NAN` when `a == f64::INFINITY` and `0.0` when `x == f64::INFINITY`
 - `Gamma::pdf` and `Gamma::ln_pdf` now return `f64::NAN` if any of `shape`, `rate`, or `x` are `f64::INFINITY`
 - `Binomial::pdf` and `Binomial::ln_pdf` now panic if `x > n` or `x < 0`
 - `Bernoulli::pdf` and `Bernoulli::ln_pdf` now panic if `x > 1` or `x < 0`

v0.4.0
- Implemented the `exponential::integral` special function
- Implemented the `Cauchy` (otherwise known as the `Lorenz`) distribution
- Implemented the `Dirichlet` distribution
- `Continuous` and `Discrete` traits no longer dependent on `Distribution` trait

v0.3.2
- Implemented the `FisherSnedecor` (F) distribution

v0.3.1
- Removed print statements from `ln_pdf` method in `Beta` distribution

v0.3.0
- Moved methods `min` and `max` out of trait `Univariate` into their own respective traits `Min` and `Max`
- Traits `Min`, `Max`, `Mean`, `Variance`, `Entropy`, `Skewness`, `Median`, and `Mode` moved from `distribution` module to `statistics` module
- `Mean`, `Variance`, `Entropy`, `Skewness`, `Median`, and `Mode` no longer depend on `Distribution` trait
- `Mean`, `Variance`, `Skewness`, and `Mode` are now generic over only one type, the return type, due to not depending on `Distribution` anymore
- `order_statistic`, `median`, `quantile`, `percentile`, `lower_quartile`, `upper_quartile`, `interquartile_range`, and `ranks` methods removed
    from `Statistics` trait. 
- `min`, `max`, `mean`, `variance`, and `std_dev` methods added to `Statistics` trait
- `Statistics` trait now implemented for all types implementing `IntoIterator` where `Item` implements `Borrow<f64>`. Slice now implicitly implements
    `Statistics` through this new implementation.
- Slice still implements `Min`, `Max`, `Mean`, and `Variance` but now through the `Statistics` implementation rather than its own implementation
- `InplaceStatistics` renamed to `OrderStatistics`, all methods in `InplaceStatistics` have `_inplace` trimmed from method name.
- Inverse DiGamma function implemented with signature `gamma::inv_digamma(x: f64) -> f64`

v0.2.0
- Created `statistics` module and `Statistics` trait
- `Statistics` trait implementation for `[f64]`
- Implemented `Beta` distribution
- Added `Modulus` trait and implementations for `f32`, `f64`, `i32`, `i64`, `u32`, and `u64` in `euclid` module
- Added periodic and sinusoidal vector generation functions in `generate` module
