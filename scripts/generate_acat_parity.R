#!/usr/bin/env Rscript
# =============================================================================
# generate_acat_parity.R
# =============================================================================
# Reference values for the meta ACAT over masks (acat_masks_p, issue #196),
# stored in tests/unit/test_meta_multimask.py::R_ACAT, and for the weighted
# single-cohort cauchy_combination (issue #216), stored in
# tests/unit/test_association_acat.py::R_ACAT_WEIGHTED. The ACAT() function below
# is copied verbatim from the reference implementation
# https://github.com/yaowuliu/ACAT, R/ACAT.R at commit e737d7d (Liu et al. 2019),
# so no R package is needed.
#
# Usage (from the repository root):
#   Rscript scripts/generate_acat_parity.R
# =============================================================================

ACAT<-function(Pvals,weights=NULL,is.check=TRUE){
    Pvals<-as.matrix(Pvals)
    if (is.check){
        #### check if there is NA
        if (sum(is.na(Pvals))>0){
            stop("Cannot have NAs in the p-values!")
        }
        #### check if Pvals are between 0 and 1
        if ((sum(Pvals<0)+sum(Pvals>1))>0){
            stop("P-values must be between 0 and 1!")
        }
        #### check if there are pvals that are either exactly 0 or 1.
        is.zero<-(colSums(Pvals==0)>=1)
        is.one<-(colSums(Pvals==1)>=1)
        if (sum((is.zero+is.one)==2)>0){
            stop("Cannot have both 0 and 1 p-values in the same column!")
        }

        if (sum(is.zero)>0){
            warning("There are p-values that are exactly 0!")
        }
        if (sum(is.one)>0){
            warning("There are p-values that are exactly 1!")
        }

    }
    #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
    if (is.null(weights)){
        is.weights.null<-TRUE
    }else{
        is.weights.null<-FALSE
        weights<-as.matrix(weights)
        if (sum(dim(weights)!=dim(Pvals))>0){
            stop("The dimensions of weights and Pvals must be the same!")
        }else if (is.check & (sum(weights<0)>0)){
            stop("All the weights must be nonnegative!")
        }else{
            w.sum<-colSums(weights)
            if (sum(w.sum<=0)>0){
                stop("At least one weight should be positive in each column!")
            }else{
                for (j in 1:ncol(weights)){
                    weights[,j]<-weights[,j]/w.sum[j]
                }
            }
        }

    }

    #### check if there are very small non-zero p values and calcuate the cauchy statistics
    is.small<-(Pvals<1e-15)
    if (is.weights.null){
         Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
         Pvals[is.small]<-1/Pvals[is.small]/pi
         cct.stat<-colMeans(Pvals)
    }else{
         Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
         Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
         cct.stat<-colSums(Pvals)
    }
    #### return the ACAT p value(s).
    pval<-pcauchy(cct.stat,lower.tail = F)
    return(pval)
}

cases <- list(
  c(1e-16, 1.0),
  c(0.01, 1.0),
  c(0.01, 0.2),
  c(1e-20, 0.5),
  c(0.3, 0.5, 0.9),
  c(1e-300, 1.0),
  c(0.999, 1.0)
)
for (p in cases) {
  cat(sprintf("%s -> %.17g\n", paste(format(p, digits = 17), collapse = ","),
              suppressWarnings(ACAT(p))))
}

# Weighted cases (COAST omnibus weights 1 x 6, 6) for cauchy_combination (#216)
w7 <- c(1, 1, 1, 1, 1, 1, 6)
weighted <- list(
  list(c(0.01, 0.05, 0.001, 0.1, 0.005, 0.02, 1e-200), w7),
  list(c(0.01, 1.0, 0.3, 1.0, 0.5, 0.2, 0.04), w7),
  list(c(1e-6, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), w7),
  list(c(0.0, 0.5), c(1, 1)),
  list(c(1.0, 1.0), c(1, 1)),
  list(c(5e-16, 0.5), c(1, 3))
)
cat("# weighted\n")
for (pw in weighted) {
  cat(sprintf("%s | %s -> %.17g\n", paste(format(pw[[1]], digits = 17), collapse = ","),
              paste(pw[[2]], collapse = ","), suppressWarnings(ACAT(pw[[1]], pw[[2]]))))
}
