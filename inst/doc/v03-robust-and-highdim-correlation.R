## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## -----------------------------------------------------------------------------
library(matrixCorr)

set.seed(20)
Y <- data.frame(
  x1 = rnorm(60),
  x2 = rnorm(60),
  x3 = rnorm(60),
  x4 = rnorm(60)
)

idx <- sample.int(nrow(Y), 5)
Y$x1[idx] <- Y$x1[idx] + 8
Y$x2[idx] <- Y$x2[idx] - 6

R_pear <- pearson_corr(Y)
R_bicor <- bicor(Y)
R_pb <- pbcor(Y)
R_win <- wincor(Y)
R_skip <- skipped_corr(Y)

summary(R_pear)
summary(R_bicor)

## -----------------------------------------------------------------------------
summary(R_skip)

## -----------------------------------------------------------------------------
fit_bicor_ci <- bicor(Y, ci = TRUE)
summary(fit_bicor_ci)

fit_pb_inf <- pbcor(Y, ci = TRUE, p_value = TRUE, n_boot = 200, seed = 1)
summary(fit_pb_inf)

fit_win_inf <- wincor(Y, ci = TRUE, p_value = TRUE, n_boot = 200, seed = 1)
summary(fit_win_inf)

fit_skip_inf <- skipped_corr(Y, ci = TRUE, p_value = TRUE, n_boot = 200, seed = 1)
summary(fit_skip_inf)

## -----------------------------------------------------------------------------
set.seed(21)
n <- 300
x2 <- rnorm(n)
x1 <- 0.8 * x2 + rnorm(n, sd = 0.4)
x3 <- 0.8 * x2 + rnorm(n, sd = 0.4)
x4 <- 0.7 * x3 + rnorm(n, sd = 0.5)
x5 <- rnorm(n)
x6 <- rnorm(n)
X <- data.frame(x1, x2, x3, x4, x5, x6)

fit_pcor_sample <- pcorr(X, method = "sample")
fit_pcor_oas <- pcorr(X, method = "oas")
R_raw <- pearson_corr(X)

round(c(
  raw_x1_x3 = R_raw["x1", "x3"],
  partial_x1_x3 = fit_pcor_sample$pcor["x1", "x3"]
), 2)

print(fit_pcor_sample, digits = 2)
summary(fit_pcor_oas)

## -----------------------------------------------------------------------------
fit_pcor_inf <- pcorr(
  X,
  method = "sample",
  return_p_value = TRUE,
  ci = TRUE
)

summary(fit_pcor_inf)

## -----------------------------------------------------------------------------
set.seed(22)
n <- 40
p_block <- 30

make_block <- function(f) {
  sapply(seq_len(p_block), function(j) 0.8 * f + rnorm(n, sd = 0.8))
}

f1 <- rnorm(n)
f2 <- rnorm(n)
f3 <- rnorm(n)

Xd <- cbind(make_block(f1), make_block(f2), make_block(f3))
p <- ncol(Xd)
colnames(Xd) <- paste0("G", seq_len(p))

R_raw <- stats::cor(Xd)
fit_shr <- shrinkage_corr(Xd)
block <- rep(1:3, each = p_block)
within <- outer(block, block, "==") & upper.tri(R_raw)
between <- outer(block, block, "!=") & upper.tri(R_raw)

round(c(
  raw_within = mean(abs(R_raw[within])),
  raw_between = mean(abs(R_raw[between])),
  shrink_within = mean(abs(fit_shr[within])),
  shrink_between = mean(abs(fit_shr[between]))
), 3)

print(fit_shr, digits = 2, max_rows = 6, max_vars = 6)
summary(fit_shr)

