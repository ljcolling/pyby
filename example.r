library(reticulate)
pyby <- reticulate::import("pyby")


list_vec_to_dataframe <- function(x, col_names) {
  n_cols <- length(x)
  x |>
    unlist() |>
    matrix(ncol = n_cols) |>
    t() |>
    as.data.frame() |>
    setNames(col_names)
}

prior <- list(family = "Cauchy", params = c(0, sqrt(2) / 2), alternative = "two.sided")
effsize <- 0.4
B <- 1000
sampling_rule <- list(n_min = 50L, n_max = 300L, step_size = 5L)


set.seed(123)
tictoc::tic()
random_values <- rnorm(B * sampling_rule$n_max, effsize, 1)
output <- pyby$sim_withrandom(random_values, prior, effsize, sampling_rule)


output <- lapply(output, unlist)
output_df <- list_vec_to_dataframe(output, c("id", "true.ES", "n", "logBF", "emp.ES", "statistic", "p.value"))


settings <- list(
  n.min = sampling_rule$n_min,
  n.max = sampling_rule$n_max,
  design = "sequential",
  prior = list(prior$famly, list(prior.location = prior$params[1], prior.scale = prior$params[2])),
  boundary = c(0, Inf),
  alternative = prior$alternative,
  type = "paired.t",
  options.sample = list(),
  extra = list(),
  packageVersion = "0.5.0"
)

sim.H1 <- list(sim = output_df, settings = settings)
class(sim.H1) <- "BFDA"
timing <- tictoc::toc()

rust <- timing$toc - timing$tic
(rust |> lubridate::seconds()) |> lubridate::as.duration()

BFDA::BFDA.analyze(sim.H1, boundary = 6) |> print()






