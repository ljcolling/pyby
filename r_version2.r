require(future)
require(furrr)
require(tidyverse)
plan("multicore")
require(BayesFactor)

effect_size <- 0.5
n_min <- 10
n_max <- 70
step_size <- 2
alternative <- "two.sided"
comparison_type <- "t.paired"
# prior = {"family": "Cauchy", "params": [0, math.sqrt(2)/2]}
samples <- (16 * 70) * 10

welfords <- function(x, new_value) {
  count <- x[1]
  mean <- x[2]
  squared_distances <- x[3]
  group <- x[4]
  count <- count + 1
  delta <- new_value - mean
  mean <- mean + (delta / count)
  delta2 <- new_value - mean
  squared_distances <- squared_distances + (delta * delta2)
  c(count, mean, squared_distances, group)
}

start_time <- Sys.time()
tibble(
  group = sort(rep(seq(1:samples), n_max)),
  value = rnorm(n_max * samples, mean = 0.5, sd = 1),
) |>
  group_by(group) |>
  group_split() |>
  future_map(function(x) {
    accumulate(x$value, welfords,
      .init = c(0, 0, 0, x$group[[1]])
    )
  }, .progress = TRUE) |>
  flatten() |>
  unlist() |>
  matrix(nrow = 4) |>
  t() |>
  as.data.frame() |>
  setNames(c("n", "mean", "sq", "group")) |>
  mutate(
    sd = sqrt(sq / (n - 1)),
    se = sd / sqrt(n),
    t = mean / se
  ) |>
  filter(n >= n_min, n %% step_size == 0) |>
  future_pmap_dfc(function(n, mean, sq, group, sd, se, t) {
    tibble(
      n = n, mean = mean, group = group, sd = sd, t = t,
      logBF = ttest.tstat(t, n, rscale = sqrt(2) / 2, simple = TRUE)
    )
  }, .progress = TRUE) |>
  mutate(p.value = 2 * pt(t, n - 1, FALSE, FALSE)) |>
  mutate(emp.ES = mean / sd, true.ES = effect_size) |>
  rename(id = group, statistic = t) |>
  select(id, true.ES, n, logBF, emp.ES, statistic, p.value) -> sim

end_time <- Sys.time()

functional_r_time <- end_time - start_time
