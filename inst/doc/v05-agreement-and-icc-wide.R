## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## -----------------------------------------------------------------------------
library(matrixCorr)

set.seed(40)
ref <- rnorm(50, mean = 100, sd = 10)
m1 <- ref + rnorm(50, sd = 2)
m2 <- ref + 1.2 + rnorm(50, sd = 3)

fit_ba <- ba(m1, m2)
fit_ccc <- ccc(data.frame(m1 = m1, m2 = m2), ci = TRUE)

print(fit_ba)
summary(fit_ccc)

## -----------------------------------------------------------------------------
wide_methods <- data.frame(
  J1 = ref + rnorm(50, sd = 1.5),
  J2 = ref + 4.0 + rnorm(50, sd = 1.8),
  J3 = ref - 3.0 + rnorm(50, sd = 2.0),
  J4 = ref + rnorm(50, sd = 1.6)
)

fit_icc_pair <- icc(
  wide_methods,
  model = "twoway_random",
  type = "agreement",
  unit = "single",
  scope = "pairwise"
)

fit_icc_overall <- icc(
  wide_methods,
  model = "twoway_random",
  type = "agreement",
  unit = "single",
  scope = "overall",
  ci = TRUE
)

print(fit_icc_pair, digits = 2)
summary(fit_icc_pair)
print(fit_icc_overall)

## -----------------------------------------------------------------------------
fit_icc_cons <- icc(
  wide_methods,
  model = "twoway_random",
  type = "consistency",
  unit = "single",
  scope = "overall",
  ci = FALSE
)

fit_icc_agr <- icc(
  wide_methods,
  model = "twoway_random",
  type = "agreement",
  unit = "single",
  scope = "overall",
  ci = FALSE
)

data.frame(
  type = c("consistency", "agreement"),
  selected_coefficient = c(
    attr(fit_icc_cons, "selected_coefficient"),
    attr(fit_icc_agr, "selected_coefficient")
  ),
  estimate = c(
    attr(fit_icc_cons, "selected_row")$estimate,
    attr(fit_icc_agr, "selected_row")$estimate
  )
)

