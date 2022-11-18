require(tidyverse)
require(bayesplay)
# load the BFDA file


bfda <- "/Users/lc663/GitHub/BFDA/dev/BFDA/package/testdating"
load(bfda)
sim.H1$sim -> bfda_data

# load the rust data

rust <- "/Users/lc663/GitHub/pyby/testdata2.csv"
rust_data <- read.csv(rust)

rust_data %>%
  as_tibble() %>%
  dplyr::mutate(emp.ES = running_mean / running_sd) %>%
  dplyr::mutate(logBF = log(bf)) %>%
  dplyr::select(item, emp.ES, logBF) %>%
  dplyr::rename(n = item) -> rust_data

bfda_data %>% select(n, emp.ES, logBF) -> bfda_data

all(bfda_data$n == rust_data$n)




purrr::map2(bfda_data$emp.ES, rust_data$emp.ES, function(x, y) {
  testthat::expect_equal(x, y, tolerance = 0.001)
})

purrr::map2(bfda_data$logBF, rust_data$logBF, function(x, y) {
  testthat::expect_equal(x, y, tolerance = 0.0001)
})



