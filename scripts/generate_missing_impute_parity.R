#!/usr/bin/env Rscript
# =============================================================================
# generate_missing_impute_parity.R
# =============================================================================
# Synthetic genotypes WITH missing calls + R SKAT reference p-values for the
# genotype-matrix missing-data policies (issue #200).
#
# Usage (from the repository root):
#   Rscript scripts/generate_missing_impute_parity.R
#
# Writes tests/fixtures/missing_impute_parity/{<trait>.csv, weights.csv, r_reference.csv}.
# Each CSV holds y and genotypes g1..g10 with NA = missing call. Per-variant
# missing rates are fixed (none equal to the cutoff); g9 is 30% missing and is
# removed by the 0.15 site cutoff in both implementations. References are computed
# from the CSVs as re-read from disk. Options:
#   SKAT_Null_Model(y ~ 1, out_type = "D" | "C", Adjustment = FALSE)
#   default policy: SKAT(Z, obj, weights = w, method = "davies", r.corr = 0 | 1,
#                        missing_cutoff = 0.15,
#                        impute.method = "bestguess" (binary; SKATBinary default)
#                                      | "fixed" (quantitative; SKAT default))
#   missing_as_ref: Z[is.na(Z)] <- 0 first (bcftools +setGT -t . -n 0), same call.
# w = dbeta(observed MAF, 1, 25), written to weights.csv so both sides use it.
# Not run in CI; used by tests/unit/test_missing_policy_parity.py.
# =============================================================================

suppressPackageStartupMessages(library(SKAT))
out <- file.path("tests", "fixtures", "missing_impute_parity")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(20261021)

n <- 200
af <- c(0.02, 0.04, 0.06, 0.08, 0.10, 0.30, 0.35, 0.70, 0.05, 0.12)
miss <- c(0.00, 0.02, 0.05, 0.08, 0.10, 0.12, 0.14, 0.10, 0.30, 0.06)
weights <- list()
refs <- list()
for (trait in c("binary", "quantitative")) {
  g <- sapply(af, function(p) rbinom(n, 2, p))
  eta <- -0.5 + 0.6 * rowSums(g[, 1:4])
  y <- if (trait == "binary") rbinom(n, 1, plogis(eta)) else round(eta + rnorm(n), 4)
  for (j in seq_along(miss)) g[sample.int(n, round(miss[j] * n)), j] <- NA
  df <- cbind(data.frame(y = y), setNames(as.data.frame(g), paste0("g", seq_along(af))))
  write.csv(df, file.path(out, paste0(trait, ".csv")), row.names = FALSE, quote = FALSE, na = "NA")

  d <- read.csv(file.path(out, paste0(trait, ".csv")))
  z <- as.matrix(d[, paste0("g", seq_along(af))])
  maf <- colMeans(z, na.rm = TRUE) / 2
  w <- dbeta(pmin(maf, 1 - maf), 1, 25)
  weights[[trait]] <- data.frame(name = trait, variant = seq_along(af), w = sprintf("%.17g", w))
  obj <- SKAT_Null_Model(d$y ~ 1, out_type = if (trait == "binary") "D" else "C",
                         Adjustment = FALSE)
  imp <- if (trait == "binary") "bestguess" else "fixed"
  z0 <- z
  z0[is.na(z0)] <- 0
  run <- function(zz, rc) {
    suppressWarnings(SKAT(zz, obj, weights = w, method = "davies", r.corr = rc,
                          missing_cutoff = 0.15, impute.method = imp))$p.value
  }
  p <- c(SKAT_default = run(z, 0), Burden_default = run(z, 1),
         SKAT_missing_as_ref = run(z0, 0), Burden_missing_as_ref = run(z0, 1))
  refs[[trait]] <- data.frame(name = trait, method = names(p), p_value = sprintf("%.15g", p))
}
write.csv(do.call(rbind, weights), file.path(out, "weights.csv"), row.names = FALSE, quote = FALSE)
ref_file <- file.path(out, "r_reference.csv")
writeLines(sprintf("# %s; SKAT %s; scripts/generate_missing_impute_parity.R",
                   R.version.string, packageVersion("SKAT")), ref_file)
suppressWarnings(write.table(do.call(rbind, refs), ref_file, sep = ",", row.names = FALSE,
                             quote = FALSE, append = TRUE))
cat("Wrote", length(refs), "specs to", out, "\n")
