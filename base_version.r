### algo
welfords <- function(x, new_value) {
  count <- x[1]
  mean <- x[2]
  squared_distances <- x[3]
  count <- count + 1
  delta <- new_value - mean
  new_mean <- mean + (delta / count)
  delta2 <- new_value - new_mean
  squared_distances <- squared_distances + (delta * delta2)
  sd <- sqrt(squared_distances / (count - 1))
  se <- sd / sqrt(count)
  t <- mean / se
  d <- mean / sd
  c(count, new_mean, squared_distances, sd, se, d, t)
}

### reducer
welfords_reduce <- function(x) {
  d <- Reduce(welfords, x$d,
    init = c(0, 0, 0, 0, 0, 0, 0),
    accumulate = TRUE
  ) |>
    list_vec_to_dataframe(c("i", "mean", "sq_dist", "sd", "se", "d", "t"))
  d$group <- x$group[[1]]
  d
}

### helpers

list_vec_to_dataframe <- function(x, col_names) {
  n_cols <- length(x)
  x |>
    unlist() |>
    matrix(ncol = n_cols) |>
    t() |>
    as.data.frame() |>
    setNames(col_names)
}

dataframe_to_list_vec <- function(x) {
  x |>
    t() |>
    as.data.frame() |>
    as.list() |>
    unname()
}

### filter
filter_base <- function(x, n_min, step_size) {
  if ((x >= n_min) & (x %% step_size == 0)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

filter_func <- function(x) filter_base(x, 10, 1)

### Bayes factor function
bayes <- function(x) {
  n <- x[1]
  t <- x[7]
  bf <- BayesFactor::ttest.tstat(t, n, rscale = .707, simple = TRUE) |>
    unname() |>
    log()
  p <- 2 * pt(abs(t), n - 1, FALSE, FALSE)
  c(x, bf, p)
}

n_max <- 20
samples <- 10
start_time <- Sys.time()
data.frame(
  d = rnorm(n_max * samples, mean = 0.5, sd = 1),
  group = sort(rep(seq(1:samples), n_max))
) |>
  split(~group) |>
  parallel::mcMap(f = welfords_reduce) |>
  Reduce(f = function(x, y) rbind(x, y)) |>
  dataframe_to_list_vec() |>
  Filter(f = function(x) filter_func(x[1])) |>
  parallel::mcMap(f = bayes) |>
  list_vec_to_dataframe(c(
    "n", "mean", "sq_dist",
    "sd", "se", "emp.ES", "statistic",
    "id", "logBF", "p"
  )) -> sim
end_time <- Sys.time()
abs(start_time - end_time)

sim$true.ES <- 0.5
sim <- sim[, c("id", "true.ES", "n", "logBF", "emp.ES", "statistic", "p")]
