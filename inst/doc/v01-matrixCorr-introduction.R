## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## -----------------------------------------------------------------------------
library(matrixCorr)

set.seed(1)
X <- as.data.frame(matrix(rnorm(120), ncol = 4))
names(X) <- paste0("V", 1:4)

fit_pearson <- pearson_corr(X, ci = TRUE)
fit_spearman <- spearman_rho(X, ci = TRUE)

print(fit_pearson, digits = 2)
summary(fit_spearman)

## -----------------------------------------------------------------------------
set.seed(2)
ref <- rnorm(40, mean = 10, sd = 2)
alt <- ref + 0.3 + rnorm(40, sd = 0.8)

fit_ba <- ba(ref, alt)
print(fit_ba)

wide_methods <- data.frame(
  m1 = ref + rnorm(40, sd = 0.2),
  m2 = ref + 0.2 + rnorm(40, sd = 0.3),
  m3 = ref - 0.1 + rnorm(40, sd = 0.4)
)

fit_ccc <- ccc(wide_methods)
fit_icc <- icc(wide_methods, scope = "pairwise")

summary(fit_ccc)
summary(fit_icc)

## -----------------------------------------------------------------------------
fit_icc_overall <- icc(wide_methods, scope = "overall", ci = TRUE)
print(fit_icc_overall)
summary(fit_icc_overall)

## -----------------------------------------------------------------------------
set.seed(3)
n_subject <- 12
n_rep <- 3

subject <- rep(seq_len(n_subject), each = n_rep)
signal <- rnorm(n_subject * n_rep)
subject_x <- rnorm(n_subject, sd = 1.2)[subject]
subject_y <- rnorm(n_subject, sd = 1.0)[subject]

dat_rm <- data.frame(
  id = subject,
  x = subject_x + signal + rnorm(n_subject * n_rep, sd = 0.2),
  y = subject_y + 0.7 * signal + rnorm(n_subject * n_rep, sd = 0.3),
  z = subject_y - 0.4 * signal + rnorm(n_subject * n_rep, sd = 0.4)
)

fit_rmcorr <- rmcorr(dat_rm, response = c("x", "y", "z"), subject = "id")
print(fit_rmcorr, digits = 2)
summary(fit_rmcorr)

## -----------------------------------------------------------------------------
set.seed(4)
n_id <- 10
n_time <- 3

dat_agree <- expand.grid(
  id = factor(seq_len(n_id)),
  time = factor(seq_len(n_time)),
  method = factor(c("A", "B"))
)

subj <- rnorm(n_id, sd = 1.0)[dat_agree$id]
subj_method <- rnorm(n_id * 2, sd = 0.2)
sm <- subj_method[(as.integer(dat_agree$id) - 1L) * 2L + as.integer(dat_agree$method)]

dat_agree$y <- subj + sm + 0.25 * (dat_agree$method == "B") +
  rnorm(nrow(dat_agree), sd = 0.35)

fit_icc_rm <- icc_rm_reml(
  dat_agree,
  response = "y",
  subject = "id",
  method = "method",
  time = "time",
  type = "consistency"
)

summary(fit_icc_rm)

## ----eval = FALSE-------------------------------------------------------------
# options(
#   matrixCorr.print_max_rows = 20L,
#   matrixCorr.print_topn = 5L,
#   matrixCorr.print_max_vars = 10L,
#   matrixCorr.print_show_ci = "yes",
#   matrixCorr.summary_max_rows = 12L,
#   matrixCorr.summary_topn = 5L,
#   matrixCorr.summary_max_vars = 10L,
#   matrixCorr.summary_show_ci = "yes"
# )

