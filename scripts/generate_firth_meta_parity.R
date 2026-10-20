#!/usr/bin/env Rscript
# =============================================================================
# generate_firth_meta_parity.R
# =============================================================================
# Synthetic strata + R references for the inverse-variance meta of per-stratum
# Firth burden log-ORs (issue #198, step 1).
#
# Usage (from the repository root):
#   Rscript scripts/generate_firth_meta_parity.R
#
# Writes tests/fixtures/firth_meta_parity/{s1,s2,s3}.csv and r_reference.csv.
# Each stratum CSV holds y, covariates c1 (continuous) and c2 (binary), and ALT
# dosages G1_1..G1_6 (rare, with an effect) and G2_1..G2_3 (ultra-rare). s2 is a
# strongly imbalanced stratum (788 cases : 101 controls) in which every G2 carrier
# is a case (complete separation). References are computed from the CSVs as re-read
# from disk:
#   w = dbeta(stratum MAF, 1, 25); burden = G %*% w
#   logistf(y ~ burden + c1 + c2, control = logistf.control(xconv = 1e-10, gconv = 1e-10))
#   beta = coef["burden"]; se = Wald SE from the inverse Fisher information
#   solve(t(X) %*% diag(p * (1 - p)) %*% X) at logistf's estimate (Heinze & Schemper
#   2002; what variantcentrifuge's Firth reports). logistf 1.26's fit$var is instead
#   the inverse of the augmented-data information (weights p(1-p)(1+h)), ~1% smaller
#   here; it is written as quantity "se_logistf_var" for the record only. Neither is
#   the profile-likelihood CI logistf reports by default.
#   metafor::rma(yi = beta, sei = se, method = "FE") over the strata of each gene.
# Not run in CI; used by tests/unit/test_meta_firth_parity.py.
# =============================================================================

suppressPackageStartupMessages({ library(logistf); library(metafor) })
out <- file.path("tests", "fixtures", "firth_meta_parity")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(20261020)

af1 <- c(0.004, 0.008, 0.012, 0.02, 0.03, 0.05)
af2 <- c(0.002, 0.003, 0.004)
make_stratum <- function(id, n_case, n_ctrl, separate_g2) {
  n <- n_case + n_ctrl
  y <- c(rep(1, n_case), rep(0, n_ctrl))
  c1 <- round(rnorm(n), 4)
  c2 <- rbinom(n, 1, 0.5)
  g1 <- sapply(af1, function(p) rbinom(n, 2, p))
  # effect: cases carry G1 alleles more often
  extra <- which(y == 1)[rbinom(n_case, 1, 0.04) == 1]
  g1[cbind(extra, sample.int(ncol(g1), length(extra), replace = TRUE))] <- 1
  g2 <- sapply(af2, function(p) rbinom(n, 2, p))
  if (separate_g2) g2[y == 0, ] <- 0
  g2[sample(which(y == 1), 2), 1] <- 1  # at least two case carriers in every stratum
  df <- data.frame(y = y, c1 = c1, c2 = c2)
  df <- cbind(df, setNames(as.data.frame(g1), paste0("G1_", seq_along(af1))))
  df <- cbind(df, setNames(as.data.frame(g2), paste0("G2_", seq_along(af2))))
  write.csv(df, file.path(out, paste0(id, ".csv")), row.names = FALSE, quote = FALSE)
}
make_stratum("s1", 300, 300, FALSE)
make_stratum("s2", 788, 101, TRUE)
make_stratum("s3", 150, 450, FALSE)

ctl <- logistf.control(xconv = 1e-10, gconv = 1e-10, maxit = 200)
rows <- list()
for (gene in c("G1", "G2")) {
  yi <- c(); sei <- c()
  for (id in c("s1", "s2", "s3")) {
    d <- read.csv(file.path(out, paste0(id, ".csv")))
    z <- as.matrix(d[, grep(paste0("^", gene, "_"), names(d))])
    maf <- colMeans(z) / 2
    w <- dbeta(pmin(maf, 1 - maf), 1, 25)
    d$burden <- as.vector(z %*% w)
    fit <- logistf(y ~ burden + c1 + c2, data = d, control = ctl, plconf = 2)
    x <- model.matrix(~ burden + c1 + c2, d)
    p <- as.vector(plogis(x %*% coef(fit)))
    b <- unname(coef(fit)["burden"]); s <- sqrt(diag(solve(crossprod(x, x * (p * (1 - p))))))[2]
    yi <- c(yi, b); sei <- c(sei, s)
    rows[[length(rows) + 1]] <- data.frame(gene = gene, stratum = id,
                                           quantity = c("beta", "se", "se_logistf_var"),
                                           value = sprintf("%.17g", c(b, s, sqrt(diag(fit$var))[2])))
  }
  m <- rma(yi = yi, sei = sei, method = "FE")
  rows[[length(rows) + 1]] <- data.frame(gene = gene, stratum = "meta", quantity = c("beta", "se", "p"),
                                         value = sprintf("%.17g", c(m$beta[1], m$se, m$pval)))
}
ref <- file.path(out, "r_reference.csv")
writeLines(sprintf("# %s; logistf %s; metafor %s; scripts/generate_firth_meta_parity.R",
                   R.version.string, packageVersion("logistf"), packageVersion("metafor")), ref)
suppressWarnings(write.table(do.call(rbind, rows), ref, sep = ",", row.names = FALSE,
                             quote = FALSE, append = TRUE))
