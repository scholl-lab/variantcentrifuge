#!/usr/bin/env Rscript
# =============================================================================
# generate_skat_robust_parity.R
# =============================================================================
# Synthetic fixtures + R SKAT reference p-values for the imbalance-robust binary
# SKAT/SKAT-O kernel (association/spa.py; port of SKAT::SKATBinary_Robust, #223).
#
# Usage (from the repository root):
#   Rscript scripts/generate_skat_robust_parity.R
#
# Writes tests/fixtures/skat_robust_parity/{<cohort>.csv, r_reference.csv,
# r_single.csv, r_moments.csv}. Each cohort CSV holds y, x1, x2, the R null-model fitted mean mu
# (%.17g, so the kernel is compared on identical inputs) and genotype columns
# <gene>_<j>. Cohorts: "r19" (n = 2000, 10% cases) and "r81" (n = 900, 800 cases).
# Genes mix MAC <= 10 (exact/ER path) and MAC > 10 (saddlepoint path) variants;
# per-variant case enrichment pushes some past |z| > 2, where SPA/ER engage.
#   SKAT_Null_Model(y ~ x1 + x2, out_type = "D", Adjustment = FALSE)
#   SKATBinary_Robust(Z, obj, method = "SKAT" | "SKATO")   (default weights.beta = c(1, 25))
# r_reference.csv: p_value from SKATBinary_Robust; for SKAT also p_davies_1e12, the
# same statistic evaluated by R's Davies at acc = 1e-12.
# r_single.csv holds p.value_singlevariant (the SPA/ER per-variant p-values);
# r_moments.csv the calibrated weighted score and covariance (Score, Phi) that
# SKATBinary_spa passes to Met_SKAT_Get_Pvalue (captured by wrapping that function).
# Not run in CI; used by tests/unit/test_skat_robust_parity.py.
# =============================================================================

suppressPackageStartupMessages(library(SKAT))
out <- file.path("tests", "fixtures", "skat_robust_parity")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(20261021)

genes <- list(
  le10 = c(1, 2, 3, 5, 7, 9),
  gt10 = c(14, 22, 35, 60),
  mixed = c(2, 4, 8, 12, 25, 40),
  single = c(6)
)
refs <- list()
single <- list()
moments <- list()
# Capture the calibrated (Score, Phi) that SKATBinary_spa hands to Met_SKAT_Get_Pvalue
cap <- new.env()
orig_get_pvalue <- SKAT:::Met_SKAT_Get_Pvalue
assignInNamespace("Met_SKAT_Get_Pvalue", function(Score, Phi, ...) {
  cap$Score <- Score
  cap$Phi <- Phi
  orig_get_pvalue(Score, Phi, ...)
}, "SKAT")
for (cohort in c("r19", "r81")) {
  n <- if (cohort == "r19") 2000 else 900
  ncase <- if (cohort == "r19") 200 else 800
  x1 <- round(rnorm(n), 4)
  x2 <- rbinom(n, 1, 0.4)
  y <- sample(rep(c(1, 0), c(ncase, n - ncase)))
  df <- data.frame(y = y, x1 = x1, x2 = x2)
  for (gname in names(genes)) {
    for (j in seq_along(genes[[gname]])) {
      mac <- genes[[gname]][j]
      enrich <- c(1, 4, 12)[(j %% 3) + 1]
      # minority-class enrichment: carriers favour the rarer phenotype class
      minority <- if (ncase < n / 2) 1 else 0
      pr <- ifelse(y == minority, enrich, 1)
      g <- rep(0, n)
      hom <- if (mac >= 12) 1 else 0
      idx <- sample.int(n, mac - hom, prob = pr)
      g[idx] <- 1
      if (hom) g[idx[1]] <- 2
      df[[sprintf("%s_%d", gname, j)]] <- g
    }
  }
  obj <- SKAT_Null_Model(y ~ x1 + x2, data = df, out_type = "D", Adjustment = FALSE)
  df$mu <- sprintf("%.17g", obj$mu)
  write.csv(df, file.path(out, paste0(cohort, ".csv")), row.names = FALSE, quote = FALSE)
  for (gname in names(genes)) {
    z <- as.matrix(df[, grep(paste0("^", gname, "_"), names(df)), drop = FALSE])
    for (method in c("SKAT", "SKATO")) {
      re <- SKATBinary_Robust(z, obj, method = method)
      p12 <- NA
      if (method == "SKAT") {
        # Same (Score, Phi), Davies at acc = 1e-12: R's acc = 1e-6 is absolute, so its
        # small p-values are Liu fallbacks (Davies returns 0 there)
        lam <- eigen(as.matrix(cap$Phi) / 2, symmetric = TRUE, only.values = TRUE)$values
        lam <- lam[lam > mean(lam[lam > 0]) / 1e5]
        p12 <- if (length(lam) == 1) re$p.value else
          SKAT:::SKAT_davies(sum(cap$Score^2) / 2, lam, acc = 1e-12, lim = 1e6)$Qq
      }
      refs[[length(refs) + 1]] <- data.frame(cohort = cohort, gene = gname, method = method,
                                             p_value = sprintf("%.15g", re$p.value),
                                             p_davies_1e12 = sprintf("%.15g", p12))
      if (method == "SKAT") {
        single[[length(single) + 1]] <- data.frame(
          cohort = cohort, gene = gname, variant = seq_len(ncol(z)),
          p_value = sprintf("%.15g", re$p.value_singlevariant))
        phi <- as.matrix(cap$Phi)
        moments[[length(moments) + 1]] <- data.frame(
          cohort = cohort, gene = gname, row = seq_len(nrow(phi)),
          score = sprintf("%.17g", as.numeric(cap$Score)),
          phi = apply(phi, 1, function(r) paste(sprintf("%.17g", r), collapse = " ")))
      }
    }
  }
}
header <- sprintf("# %s; SKAT %s; SPAtest %s; scripts/generate_skat_robust_parity.R",
                  R.version.string, packageVersion("SKAT"), packageVersion("SPAtest"))
tabs <- list(r_reference.csv = refs, r_single.csv = single, r_moments.csv = moments)
for (f in names(tabs)) {
  writeLines(header, file.path(out, f))
  tab <- do.call(rbind, tabs[[f]])
  suppressWarnings(write.table(tab, file.path(out, f), sep = ",", row.names = FALSE,
                               quote = FALSE, append = TRUE))
}
cat("Wrote", length(refs), "references to", out, "\n")
