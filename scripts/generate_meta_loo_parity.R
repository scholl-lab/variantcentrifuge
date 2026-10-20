#!/usr/bin/env Rscript
# =============================================================================
# generate_meta_loo_parity.R
# =============================================================================
# R MetaSKAT references for leave-one-stratum-out meta (issue #193), stored in
# tests/unit/test_r_parity_skato.py::R_META_LOO. Inputs are the synthetic COAST
# golden fixtures (tests/fixtures/coast_golden), split into three contiguous thirds
# of samples; covariate (i mod 5) - 2. drop=0 is the full 3-stratum meta; drop=d
# drops cohort d (= stratum s{d-1}).
#
# Usage (from the repository root; SKAT 2.2.5, MetaSKAT 0.90):
#   Rscript scripts/generate_meta_loo_parity.R [tests/fixtures/coast_golden]
#
# Output lines: "<scenario> drop=<d> burden=<p> skato=<p>" (%.10e).
# =============================================================================

suppressMessages(library(MetaSKAT))
args <- commandArgs(TRUE)
fx <- if (length(args)) args[1] else "tests/fixtures/coast_golden"
for (name in c("gckd_g1_all_categories", "gckd_g2_partial_bmv_ptv", "gckd_g4_all_bmv")) {
  Z <- as.matrix(read.csv(file.path(fx, paste0(name, "_geno.csv"))))
  y <- read.csv(file.path(fx, paste0(name, "_pheno.csv")))$pheno
  n <- length(y); i0 <- 0:(n - 1)
  b <- c(0, n %/% 3, (2 * n) %/% 3, n)
  coh <- findInterval(i0, b[2:3]) + 1  # 1,2,3 = thirds [0,n//3),[n//3,2n//3),[2n//3,n)
  cov <- (i0 %% 5) - 2
  for (d in c(0, 1, 2, 3)) {
    keep <- coh != d
    yy <- y[keep]; xx <- cov[keep]; cc <- as.integer(factor(coh[keep]))
    ns <- as.vector(table(cc))
    obj <- Meta_Null_Model(split(yy, cc), lapply(split(xx, cc), as.matrix), n.cohort = length(ns), out_type = "D")
    zz <- Z[keep, , drop = FALSE][order(cc), , drop = FALSE]
    p1 <- MetaSKAT_wZ(zz, obj, combined.weight = TRUE, weights.beta = c(1, 25), r.corr = 1)$p.value
    p3 <- MetaSKAT_wZ(zz, obj, combined.weight = TRUE, weights.beta = c(1, 25), method = "optimal")$p.value
    cat(sprintf("%s drop=%d burden=%.10e skato=%.10e\n", name, d, p1, p3))
  }
}
