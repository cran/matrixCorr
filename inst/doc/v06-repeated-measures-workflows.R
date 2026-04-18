## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## -----------------------------------------------------------------------------
library(matrixCorr)

set.seed(50)
n_id <- 14
n_time <- 4

dat <- expand.grid(
  id = factor(seq_len(n_id)),
  time = factor(seq_len(n_time)),
  method = factor(c("A", "B"))
)

dat$time_index <- as.integer(dat$time)

subj <- rnorm(n_id, sd = 1.0)[dat$id]
subject_method <- rnorm(n_id * 2, sd = 0.25)
sm <- subject_method[(as.integer(dat$id) - 1L) * 2L + as.integer(dat$method)]
subject_time <- rnorm(n_id * n_time, sd = 0.75)
st <- subject_time[(as.integer(dat$id) - 1L) * n_time + as.integer(dat$time)]

dat$y <- subj + sm + st + 0.35 * (dat$method == "B") +
  rnorm(nrow(dat), sd = 0.35)

## -----------------------------------------------------------------------------
set.seed(51)
dat_rmcorr <- data.frame(
  id = rep(seq_len(n_id), each = n_time),
  x = rnorm(n_id * n_time),
  y = rnorm(n_id * n_time),
  z = rnorm(n_id * n_time)
)

dat_rmcorr$y <- 0.7 * dat_rmcorr$x +
  rnorm(n_id, sd = 1)[dat_rmcorr$id] +
  rnorm(nrow(dat_rmcorr), sd = 0.3)

fit_rmcorr <- rmcorr(dat_rmcorr, response = c("x", "y", "z"), subject = "id")
summary(fit_rmcorr)

## -----------------------------------------------------------------------------
fit_ba_rm <- ba_rm(
  dat,
  response = "y",
  subject = "id",
  method = "method",
  time = "time_index"
)

summary(fit_ba_rm)

## -----------------------------------------------------------------------------
fit_ccc_ustat <- ccc_rm_ustat(
  dat,
  response = "y",
  subject = "id",
  method = "method",
  time = "time_index"
)

fit_ccc_reml <- ccc_rm_reml(
  dat,
  response = "y",
  subject = "id",
  method = "method",
  time = "time_index",
  ci = FALSE
)

summary(fit_ccc_ustat)
summary(fit_ccc_reml)

## -----------------------------------------------------------------------------
fit_icc_cons <- icc_rm_reml(
  dat,
  response = "y",
  subject = "id",
  method = "method",
  time = "time_index",
  type = "consistency",
  ci = TRUE
)

fit_icc_agr <- icc_rm_reml(
  dat,
  response = "y",
  subject = "id",
  method = "method",
  time = "time_index",
  type = "agreement",
  ci = FALSE
)

summary(fit_icc_cons)
summary(fit_icc_agr)

## -----------------------------------------------------------------------------
data.frame(
  method = c("Repeated CCC (REML)", "Repeated ICC (agreement, REML)"),
  estimate = c(
    fit_ccc_reml[1, 2],
    fit_icc_agr[1, 2]
  )
)

