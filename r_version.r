
require(future)
plan("multicore")

start_time <- Sys.time()
B <- 1000
N <- 200
n <- rnorm(B * N, mean = 0.5, sd = 1)
ns <- split(n, ceiling(seq_along(n) / N))

d <- function(x, new_value) {
  count <- x[1]
  mean <- x[2]
  sqquared_distances <- x[3]
  group <- x[4]
  count <- count + 1
  delta <- new_value - mean
  mean <- mean + (delta / count)
  delta2 <- new_value - mean
  sqquared_distances <- sqquared_distances + (delta * delta2)
  # sd <- (sqquared_distances / sqrt(count - 1))
  c(count, mean, sqquared_distances, group)
}

z <- furrr::future_map2(1:1000, ns, function(x, y) {
  purrr::accumulate(y, d, .init = c(0, 0, 0, x))
},
.progress = TRUE
) |>
  purrr::flatten() |>
  furrr::future_map(function(x) {
    n <- x[1]
    mean <- x[2]
    sqquared_distances <- x[3]
    group <- x[4]
    if (n < 2) {
      return(NA)
    }
    if ((n %% 10) != 0) {
      return(NA)
    }
    sd <- sqrt(sqquared_distances / (n - 1))
    se <- sd / sqrt(n)
    t <- mean / se
    bf <- BayesFactor::ttest.tstat(t, n, rscale = sqrt(2) / 2, simple = TRUE)
    c(
      id = group, true.ES = 0.5, n = n, logBF = log(bf), emp.ES = mean / sd,
      statistic = mean / se, p.value = 0
    )
  }, .progress = TRUE)


is_not_na <- function(x) {
  !is.na(x[[1]])
}

z_mat <- Reduce(function(x, y) rbind(x, y), Filter(is_not_na, z))
z_mat |> as.data.frame(row.names = F) -> z_df

end_time <- Sys.time()

end_time - start_time

