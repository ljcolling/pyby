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

bf_func <- function(n, mean, sq_dist, sd, se, d, t, group) {
  bf <- BayesFactor::ttest.tstat(t, n, rscale = sqrt(2) / 2, simple = TRUE)
  c(n, mean, sq_dist, sd, se, d, t, bf, group)
}

list_vec_to_dataframe <- function(x, col_names) {
  n_cols <- length(x)
  x |>
    unlist() |>
    matrix(ncol = n_cols) |>
    t() |>
    as.data.frame() |>
    setNames(col_names)
}

list_vec_to_tibble <- function(x, col_names) {
  n_cols <- length(x)
  x |>
    unlist() |>
    matrix(ncol = n_cols) |>
    t() |>
    tibble::as_tibble() |>
    setNames(col_names)
}


require(purrr)
require(dplyr)
require(furrr)
n_max <- 20
samples <- 10
n_min <- 10
step_size <- 1
# start_time <- Sys.time()

tibble(
  group = sort(rep(seq(1:samples), n_max)),
  value = rnorm(n_max * samples, mean = 0.5, sd = 1),
) |>
  group_by(group) |>
  group_split() |>
  future_map(function(x) {
    accumulate(x$value, welfords,
      .init = c(0, 0, 0, 0, 0, 0, 0, x$group[[1]])
    )
  }, .progress = TRUE) |>
  flatten() |>
  list_vec_to_tibble(c(
    "n", "mean", "sq_dist", "sd",
    "se", "d", "t", "group"
  )) |>
  filter(n >= n_min, n %% step_size == 0) |>
  furrr::future_pmap(\(n, mean, sq_dist, sd, se, d, t, group)
    bf_func(n, mean, sq_dist, sd, se, d, t, group)
  ) |>
  list_vec_to_tibble(c(
    "n", "mean", "sq_dist", "sd",
    "se", "d", "t", "bf", "group"
  ))


 # cars |> split(1:10) # can split dataframes!
welfords_reduce <- function(x) {
d<- Reduce(welfords, x$d, init = c(0, 0, 0, 0, 0, 0, 0, 0), accumulate = TRUE) |>
  list_vec_to_dataframe(c("i", "mean", "sq_dist", "sd", "se", "d", "t"))
  d$group <- x$group[[1]]
  d
}

select_items <- function(x) {
  if(x > 10) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

data.frame(
  d =  rnorm(n_max * samples, mean = 0.5, sd = 1),
  group = sort(rep(seq(1:samples), n_max))
) |>
  split(~group) |>
  parallel::mcMap(f = welfords_reduce) |>
  Reduce(f = function(x, y) rbind(x, y)) |>
  as.list() |>
  Filter(f = select_items,)
### <- now do the filtering
### <- and then do the bfs


  # (function(x,y) x$group = i, _, 1:10)
split(
  f = sort(rep(seq(1:samples), n_max)),
  drop = TRUE
) |>
  # Reduce(f = function(x,y) rbind(x,y))


|>
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

temp$logBF <- future_pmap_dbl(temp %>% select(t, n), function(t, n) {
  ttest.tstat(t, n, rscale = sqrt(2) / 2, simple = TRUE)
}, .progress = TRUE)

sim <- temp |>
  mutate(p.value = 2 * pt(t, n - 1, FALSE, FALSE)) |>
  mutate(emp.ES = mean / sd, true.ES = effect_size) |>
  rename(id = group, statistic = t) |>
  select(id, true.ES, n, logBF, emp.ES, statistic, p.value)

end_time <- Sys.time()

functional_r_time <- end_time - start_time
save.image("functional_r_version.Rdata")

d <- data.frame(i = c(1, 2, 3, 4, 5), n = rnorm(5))

# list_vec_to_dataframe
dataframe_to_list_vec <- function(x) { 
 x |> 
  t() |>
  as.data.frame() |>
  as.list() |> unname()
}

  |> 
  Filter(f = function(x) x[1] == 4)
