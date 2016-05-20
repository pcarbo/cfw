# SUMMARY
# -------
# This file contains functions that implement the polygenic model for
# estimating the proportion of variance expained by available genetic
# markers. Here is an overview of the functions defined in this file:
#
#   normalizelogweights(logw)
#   polygenic.model(X,y,h)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Takes as input an array of unnormalized log-importance weights and
# returns normalized importance weights such that the sum of the
# normalized importance weights is equal to one.
normalizelogweights <- function (logw) {

  # I guard against underflow or overflow by adjusting the log-importance
  # weights so that the largest importance weight is one.
  c <- max(logw)
  w <- exp(logw - c)

  # Normalize the importance weights.
  return(w / sum(w))
}

# ----------------------------------------------------------------------
# Estimates the proportion of variance explained under the polygenic
# model with a uniform prior on the proportion of variance
# explained. Technically, h is not the proportion of variance
# explained, but rather a ratio of expectations which acts as an
# unbiased estimate of the proportion of variance explained.
# Critically, this algorithm will only work correctly if y and X are
# centered so that vector y and each column of X has a mean of zero.
#
# There are two outputs: logw is the vector of log-importance weights
# for each setting of h, and sigma is the posterior estimate of the
# residual variance.
polygenic.model <- function (X, y, h) {

  # Get the number of samples (n), the number of SNPs genotyped (p),
  # and the number of hyperparameter settings (ns).
  n  <- nrow(X)
  p  <- ncol(X)
  ns <- length(h)

  # Get settings for the prior variance of the "background" polygenic
  # effects.
  cat("Acquiring prior settings.\n")
  sx <- sum(apply(X,2,sd)^2)
  sa <- p*h/(1-h)/sx
  
  # Compute the kinship matrix.
  cat("Computing kinship matrix.\n")
  K <- tcrossprod(X)/p

  # Initialize storage for the log-importance weights (logw), the
  # posterior estimates of the residual variance (sigma), and the
  # posterior estimates of the proportion of variance explained (hpe).
  logw  <- rep(-Inf,length(h))
  sigma <- rep(NA,length(h))
  hpe   <- rep(NA,length(h))
  
  # Compute the log-importance weights.
  cat("Computing importance weights for",ns,"hyperparameter settings.\n")
  for (i in 1:ns) {
    caterase(sprintf("(%d) h = %0.3f (sa/p = %0.4f)  ",i,h[i],sqrt(sa[i]/p)))

    # Compute the covariance of Y (divided by residual variance) and
    # its factors. If H is not positive definite, L = FALSE.
    H <- diag(n) + sa[i]*K
    L <- t(tryCatch(chol(H),error = function(e) FALSE))

    # Compute the log-importance weight and the posterior mean of the
    # residual variance.
    if (is.matrix(L)) {
      b        <- sum(y*solve(H,y))
      logw[i]  <- -determinant(b*H,logarithm = TRUE)$modulus/2
      sigma[i] <- b/(n - 2)
    }
  }
  cat("\n")

  return(list(logw = logw,sigma = sigma))
}
