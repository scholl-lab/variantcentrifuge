#!/usr/bin/env Rscript
# =============================================================================
# generate_skat_null_parity.R
# =============================================================================
# Synthetic fixtures + R SKAT reference p-values for the Python SKAT backend's
# null model / score-test conventions (issues #179, #180, #182).
#
# Usage (from the repository root):
#   Rscript scripts/generate_skat_null_parity.R
#
# Writes tests/fixtures/skat_null_parity/{<spec>.csv, weights.csv, r_reference.csv}.
# Each spec CSV holds y, covariates x1..x3 (cov3 specs only) and genotypes g1..g12.
# References are computed from the CSVs as re-read from disk, so R and Python
# see identical inputs. Options are matched explicitly:
#   SKAT_Null_Model(y ~ 1 | y ~ x, out_type = "D" | "C", Adjustment = FALSE)
#   SKAT  : SKAT(g, obj, weights = w, method = "davies", r.corr = 0)
#   Burden: SKAT(g, obj, weights = w, method = "davies", r.corr = 1)
#   SKATO : SKAT(g, obj, weights = w, method = "optimal.adj")
#   LinearBurden (quantitative only): lm(y ~ b [+ x]) with b = g %*% w; rows
#   LinearBurden (p), LinearBurden_beta, LinearBurden_se for the burden term.
# Not run in CI; used by tests/unit/test_skat_null_parity.py.
# =============================================================================

suppressPackageStartupMessages(library(SKAT))
out <- file.path("tests", "fixtures", "skat_null_parity")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(20261020)

maf <- seq(0.03, 0.22, length.out = 12)
weights <- list()
refs <- list()
for (trait in c("binary", "quantitative")) for (n in c(60, 300)) for (nc in c(0, 3)) for (sig in c(0, 1)) {
  name <- sprintf("%s_n%d_cov%d_sig%d", trait, n, nc, sig)
  x <- matrix(round(rnorm(n * 3), 4), n, 3)
  g <- sapply(maf, function(p) rbinom(n, 2, p))
  eta <- 0.3 * x[, 1] - 0.25 * x[, 2] + sig * 0.65 * rowSums(g[, 1:4])
  y <- if (trait == "binary") rbinom(n, 1, plogis(eta)) else round(eta + 2 * rnorm(n), 4)
  df <- data.frame(y = y)
  if (nc > 0) df <- cbind(df, setNames(as.data.frame(x), paste0("x", 1:3)))
  df <- cbind(df, setNames(as.data.frame(g), paste0("g", 1:12)))
  write.csv(df, file.path(out, paste0(name, ".csv")), row.names = FALSE, quote = FALSE)

  d <- read.csv(file.path(out, paste0(name, ".csv")))
  gg <- as.matrix(d[, paste0("g", 1:12)])
  w <- dbeta(colMeans(gg) / 2, 1, 25)
  weights[[name]] <- data.frame(name = name, variant = 1:12, w = sprintf("%.17g", w))
  yy <- d$y
  if (nc > 0) {
    xx <- as.matrix(d[, paste0("x", 1:3)])
    obj <- SKAT_Null_Model(yy ~ xx, out_type = if (trait == "binary") "D" else "C", Adjustment = FALSE)
  } else {
    obj <- SKAT_Null_Model(yy ~ 1, out_type = if (trait == "binary") "D" else "C", Adjustment = FALSE)
  }
  p <- c(
    SKAT = SKAT(gg, obj, weights = w, method = "davies", r.corr = 0)$p.value,
    Burden = SKAT(gg, obj, weights = w, method = "davies", r.corr = 1)$p.value,
    SKATO = SKAT(gg, obj, weights = w, method = "optimal.adj")$p.value
  )
  if (trait == "quantitative") {
    b <- as.numeric(gg %*% w)
    fit <- if (nc > 0) lm(yy ~ b + xx) else lm(yy ~ b)
    co <- summary(fit)$coefficients["b", ]
    p <- c(p, LinearBurden = co[["Pr(>|t|)"]], LinearBurden_beta = co[["Estimate"]],
           LinearBurden_se = co[["Std. Error"]])
  }
  refs[[name]] <- data.frame(name = name, method = names(p), p_value = sprintf("%.15g", p))
}
write.csv(do.call(rbind, weights), file.path(out, "weights.csv"), row.names = FALSE, quote = FALSE)
ref_file <- file.path(out, "r_reference.csv")
writeLines(sprintf("# %s; SKAT %s; scripts/generate_skat_null_parity.R",
                   R.version.string, packageVersion("SKAT")), ref_file)
suppressWarnings(write.table(do.call(rbind, refs), ref_file, sep = ",", row.names = FALSE,
                             quote = FALSE, append = TRUE))
cat("Wrote", length(refs), "specs to", out, "\n")
