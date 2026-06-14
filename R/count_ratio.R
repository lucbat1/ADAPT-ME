
#' @importFrom stats pchisq optim pnorm dnorm lm.fit sd

log_sum_exp <- function(x) {
  max_x <- max(x)
  if (!is.finite(max_x)) {
    return(max_x)
  }
  max_x + log(sum(exp(x - max_x)))
}

gauss_hermite_rule <- function(n = 15) {
  i <- seq_len(n - 1)
  jacobi <- matrix(0, n, n)
  offdiag <- sqrt(i / 2)
  jacobi[cbind(i, i + 1)] <- offdiag
  jacobi[cbind(i + 1, i)] <- offdiag
  eig <- eigen(jacobi, symmetric = TRUE)
  ord <- order(eig$values)
  list(
    nodes = eig$values[ord],
    weights = sqrt(pi) * eig$vectors[1, ord]^2
  )
}

tobit_mixed_loglik <- function(par, y, censor_point, censored, X, subject, ghq) {
  p <- ncol(X)
  beta <- par[seq_len(p)]
  sigma <- exp(par[p + 1])
  subject_sd <- exp(par[p + 2])
  split_index <- split(seq_along(y), subject)
  log_subject_lik <- vapply(split_index, function(idx) {
    node_lik <- vapply(seq_along(ghq$nodes), function(k) {
      random_intercept <- sqrt(2) * subject_sd * ghq$nodes[k]
      eta <- drop(X[idx, , drop = FALSE] %*% beta) + random_intercept
      exact <- !censored[idx]
      exact_ll <- if (any(exact)) {
        sum(dnorm(y[idx][exact], mean = eta[exact], sd = sigma, log = TRUE))
      } else {
        0
      }
      censored_ll <- if (any(!exact)) {
        sum(pnorm(censor_point[idx][!exact], mean = eta[!exact], sd = sigma, log.p = TRUE))
      } else {
        0
      }
      log(ghq$weights[k]) + exact_ll + censored_ll
    }, numeric(1))
    log_sum_exp(node_lik) - 0.5 * log(pi)
  }, numeric(1))
  sum(log_subject_lik)
}

fit_tobit_mixed <- function(y, censor_point, censored, X, subject, ghq_order = 15) {
  X <- as.matrix(X)
  subject <- droplevels(factor(subject))
  if (ncol(X) == 0 || nrow(X) != length(y) || length(y) != length(subject)) {
    stop("Invalid model inputs")
  }
  if (qr(X)$rank < ncol(X)) {
    stop("The model matrix is rank deficient")
  }

  observed_y <- y
  lm_fit <- lm.fit(x = X, y = observed_y)
  beta_start <- lm_fit$coefficients
  beta_start[is.na(beta_start)] <- 0
  beta_start <- beta_start[seq_len(ncol(X))]
  resid_start <- lm_fit$residuals
  sigma_start <- sd(resid_start)
  if (!is.finite(sigma_start) || sigma_start <= 0) {
    sigma_start <- sd(observed_y)
  }
  if (!is.finite(sigma_start) || sigma_start <= 0) {
    sigma_start <- 1
  }
  subject_means <- tapply(resid_start, subject, mean)
  subject_sd_start <- sd(subject_means)
  if (!is.finite(subject_sd_start) || subject_sd_start <= 0) {
    subject_sd_start <- sigma_start / 2
  }

  ghq <- gauss_hermite_rule(ghq_order)
  start <- c(beta_start, log(sigma_start), log(subject_sd_start))
  objective <- function(par) {
    ll <- tobit_mixed_loglik(par, y, censor_point, censored, X, subject, ghq)
    if (!is.finite(ll)) {
      return(.Machine$double.xmax^0.5)
    }
    -ll
  }
  fit <- optim(start, objective, method = "BFGS", control = list(maxit = 300))
  if (fit$convergence != 0 || !is.finite(fit$value)) {
    stop("Tobit mixed-effects optimization failed")
  }
  list(
    beta = fit$par[seq_len(ncol(X))],
    logLik = -fit$value,
    convergence = fit$convergence
  )
}

mean_within_subject_sd <- function(y, subject) {
  subject_sd <- tapply(y, subject, function(values) {
    if (length(values) < 2) {
      return(NA_real_)
    }
    sd(values)
  })
  mean(subject_sd, na.rm = TRUE)
}

count_ratio <- function(count_table, design_matrix, subject, reftaxa=NULL, censor=1, test_all=FALSE){
  alltaxa <- colnames(count_table)
  if (is.null(reftaxa)) reftaxa <- alltaxa
  TBDtaxa <- NULL
  if (test_all){ # examine count ratio between other taxa and the subset
    TBDtaxa <- alltaxa
  } else{ # relative abundance within the subset of taxa
    TBDtaxa <- reftaxa
  }

  refcounts <- rowSums(count_table[, reftaxa, drop=FALSE]) # sum of reference taxa as denominator
  null_filter <- refcounts != 0
  refcounts <- refcounts[null_filter]
  subset_designmatrix <- as.matrix(design_matrix[null_filter, , drop=FALSE])
  subset_subject <- droplevels(factor(subject[null_filter]))
  TBD_counts <- count_table[null_filter, TBDtaxa, drop=FALSE]
  existence <- 1*(TBD_counts > 0)
  TBD_counts[TBD_counts == 0] <- censor
  prevalences <- colMeans(existence)

  CR_result <- data.frame(Taxa = TBDtaxa,
                          prevalence=prevalences,
                          log10foldchange=NA_real_,
                          teststat=NA_real_,
                          pval=NA_real_,
                          within_subject_sd=NA_real_,
                          stringsAsFactors = FALSE)

  n_taxa <- ncol(TBD_counts)
  condition_col <- match("condition", colnames(subset_designmatrix))
  if (is.na(condition_col)) {
    stop("The design matrix must contain a 'condition' column")
  }
  null_designmatrix <- subset_designmatrix[, -condition_col, drop=FALSE]

  cat("\nRunning mixed-effects Tobit models across taxa...\n")
  for(i in seq_len(n_taxa)) {
    y_vec <- as.numeric(log(TBD_counts[, i]) - log(refcounts))
    censor_point <- log(censor) - log(refcounts)
    censored <- existence[, i] == 0
    CR_result$within_subject_sd[i] <- mean_within_subject_sd(y_vec, subset_subject)

    fit_result <- tryCatch({
      full_fit <- fit_tobit_mixed(y_vec, censor_point, censored, subset_designmatrix, subset_subject)
      null_fit <- fit_tobit_mixed(y_vec, censor_point, censored, null_designmatrix, subset_subject)
      lrt <- max(0, 2 * (full_fit$logLik - null_fit$logLik))
      list(beta = full_fit$beta[condition_col], teststat = lrt)
    }, error = function(e) {
      NULL
    })

    if (!is.null(fit_result)) {
      CR_result$log10foldchange[i] <- fit_result$beta / log(10)
      CR_result$teststat[i] <- fit_result$teststat
      CR_result$pval[i] <- 1 - pchisq(fit_result$teststat, 1)
    }
  }

  failed_models <- sum(is.na(CR_result$pval))
  if (failed_models > 0) {
    warning(sprintf("%d taxa could not be fit and were returned with NA p-values", failed_models))
  }

  return(CR_result)
}
