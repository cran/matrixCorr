## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## -----------------------------------------------------------------------------
library(matrixCorr)

set.seed(10)
z <- rnorm(80)
u <- rnorm(80)
X <- data.frame(
  x1 = z + rnorm(80, sd = 0.35),
  x2 = 0.85 * z + rnorm(80, sd = 0.45),
  x3 = 0.25 * z + 0.70 * u + rnorm(80, sd = 0.45),
  x4 = rnorm(80)
)

## -----------------------------------------------------------------------------
R_pear <- pearson_corr(X)
R_spr  <- spearman_rho(X)
R_ken  <- kendall_tau(X)
R_dcor <- dcor(X)

print(R_pear, digits = 2)
summary(R_spr)

## -----------------------------------------------------------------------------
plot(R_pear)

## -----------------------------------------------------------------------------
R_pear_ci <- pearson_corr(X, ci = TRUE)
summary(R_pear_ci)

## -----------------------------------------------------------------------------
set.seed(11)
x <- sort(rnorm(60))
y <- x^3 + rnorm(60, sd = 0.5)
dat_mon <- data.frame(x = x, y = y)

pearson_corr(dat_mon)
spearman_rho(dat_mon)
kendall_tau(dat_mon)

## -----------------------------------------------------------------------------
fit_spr_ci <- spearman_rho(X, ci = TRUE)
fit_ken_ci <- kendall_tau(X, ci = TRUE)

summary(fit_spr_ci)
summary(fit_ken_ci)

## -----------------------------------------------------------------------------
set.seed(12)
x <- runif(100, -2, 2)
y <- x^2 + rnorm(100, sd = 0.2)
dat_nonlin <- data.frame(x = x, y = y)

pearson_corr(dat_nonlin)
dcor(dat_nonlin)

## -----------------------------------------------------------------------------
fit_dcor_p <- dcor(dat_nonlin, p_value = TRUE)
summary(fit_dcor_p)

## -----------------------------------------------------------------------------
X_miss <- X
X_miss$x2[c(3, 7)] <- NA

try(pearson_corr(X_miss))
pearson_corr(X_miss, na_method = "pairwise")

