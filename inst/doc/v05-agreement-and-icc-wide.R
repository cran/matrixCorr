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
estimate(fit_ccc)
confint(fit_ccc)
ci(fit_ccc)
tidy(fit_ccc)

## -----------------------------------------------------------------------------
fit_ba_pairwise <- ba(data.frame(m1 = m1, m2 = m2, m3 = ref - 0.8 + rnorm(50, sd = 2.5)))
print(fit_ba_pairwise)
summary(fit_ba_pairwise)

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

## ----echo = FALSE-------------------------------------------------------------
agreement_summary <- data.frame(
  Function = c(
    "`ccc()`",
    "`ccc_rm_reml()`",
    "`ccc_rm_ustat()`",
    "`icc()`",
    "`icc_rm_reml()`",
    "`cia()`",
    "`cia_rm()`"
  ),
  `Index family` = c(
    "Lin CCC",
    "Repeated-measures CCC",
    "Repeated-measures CCC",
    "Classical ICC",
    "Repeated-measures ICC",
    "CIA",
    "Repeated-measures CIA"
  ),
  `Target question` = c(
    "Pairwise concordance combining precision and accuracy",
    "Pairwise repeated-measures concordance from fitted variance components",
    "Pairwise repeated-measures concordance from U-statistic distances",
    "Reliability or agreement under the selected ICC model, type, and unit",
    "Pairwise repeated-measures reliability/agreement from fitted variance components",
    "Individual agreement/interchangeability relative to within-method replicate disagreement",
    "Individual agreement/interchangeability across matched repeated conditions"
  ),
  `Data/design supported` = c(
    "Numeric wide data; rows are paired observational units",
    "Long repeated-measures data; subject plus optional method/time structure",
    "Long repeated-measures data; balanced method/time coverage per pair",
    "Numeric wide data; pairwise or overall all-column scope",
    "Long repeated-measures data; subject plus optional method/time structure",
    "Long replicated method-comparison data with replicate identifiers",
    "Long matched repeated-measures data with one observation per subject-method-condition cell"
  ),
  `Estimation approach` = c(
    "Moment CCC with Lin delta-method/Fisher-z CI when requested",
    "REML mixed-model variance components and fixed-effect bias term",
    "Nonparametric U-statistic estimator with Fisher-z CI when requested",
    "Classical ANOVA mean-square formulas",
    "REML mixed-model variance components and fixed-effect bias term",
    "Method-of-moments disagreement ratios; optional bounded variance-component variant",
    "Categorical-condition ANOVA estimator"
  ),
  `Inference reported` = c(
    "Estimate; optional confidence interval",
    "Estimate; optional confidence interval; variance-component diagnostics",
    "Estimate; optional confidence interval",
    "Pairwise: estimate and optional CI. Overall: coefficient table with F statistic, df, p-value, and optional CI",
    "Estimate; optional confidence interval; variance-component diagnostics",
    "Estimate; optional confidence interval",
    "Estimate; optional confidence interval; homogeneity F statistic and p-value"
  ),
  `Formal hypothesis test implemented? Yes/No/Requires verification` = c(
    "No",
    "No for the CCC parameter; variance-component selection tests are not CCC tests",
    "No",
    "Requires verification for overall p-values; no pairwise ICC test is documented",
    "No for the ICC parameter; variance-component selection tests are not ICC tests",
    "No",
    "Yes: homogeneity of agreement across conditions"
  ),
  check.names = FALSE
)

knitr::kable(agreement_summary, escape = FALSE)

