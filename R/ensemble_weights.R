#' Ensemble TWAS Weights via Stacked Regression
#'
#' Given cross-validated predictions from multiple TWAS weight methods, learns
#' non-negative combination coefficients (summing to 1) via constrained least
#' squares. Returns ensemble weights and per-method performance metrics.
#'
#' This implements the stacked regression approach of SR-TWAS (Dai et al.,
#' Nature Communications, 2024, \doi{10.1038/s41467-024-50983-w}). The ensemble
#' provides a principled way to combine predictions from many TWAS weight
#' methods without requiring the user to pick one method a priori or pay a
#' multiple-testing penalty for running several.
#'
#' For single-dataset usage, pass one \code{twas_weights_cv()} result directly.
#' For multi-dataset ensemble (e.g., combining cell types or reference panels
#' such as CUMC1 + MIT), pass a list of \code{twas_weights_cv()} results along
#' with a list of observed Y vectors — this learns a single joint set of
#' coefficients.
#'
#' @param cv_results Output of \code{\link{twas_weights_cv}}, with \code{$prediction}
#'   (named list of method -> out-of-fold prediction matrix, keys like
#'   \code{"susie_predicted"}). For multi-dataset: a list of such objects.
#' @param Y Observed outcome vector or matrix (samples x contexts). For
#'   multi-dataset: a list of vectors/matrices, one per dataset.
#' @param twas_weight_list Optional named list of weight matrices from
#'   \code{\link{twas_weights}}, with keys like \code{"susie_weights"}. Used to
#'   construct the final combined TWAS weight vector. For multi-dataset: a list
#'   of such lists (the first is used as the weight template).
#' @param context_index Integer indicating which column of Y to use when Y is a
#'   matrix. Default is 1 (univariate).
#'
#' @return A list with components:
#' \describe{
#'   \item{method_coef}{Named numeric vector of combination coefficients
#'     (\eqn{\zeta_k}), non-negative and summing to 1. Names are method
#'     base names (e.g., \code{"susie"}, \code{"enet"}).}
#'   \item{ensemble_twas_weights}{Final combined weight vector
#'     \eqn{w = \sum_k \zeta_k w_k}, or NULL if \code{twas_weight_list}
#'     is not provided. Returned as a vector for univariate Y, matrix otherwise.}
#'   \item{method_performance}{Named numeric vector of per-method R-squared
#'     computed from out-of-fold CV predictions. Preserved so users can still
#'     report individual method performance.}
#'   \item{ensemble_performance}{Named numeric vector with \code{corr}
#'     (Pearson correlation) and \code{rsq} (R-squared) for the ensemble.}
#' }
#'
#' @details
#' The stacked regression solves:
#' \deqn{\min_{\zeta} \|y - P\zeta\|^2 \quad \text{s.t.} \quad \zeta_k \geq 0,\ \sum_k \zeta_k = 1}
#' where P is the \eqn{n \times K} matrix of out-of-fold predictions from K
#' methods. The constrained quadratic program is solved via
#' \code{quadprog::solve.QP}. If the QP solver fails (e.g., due to numerical
#' issues with highly collinear predictions), the function falls back to equal
#' weights with a message.
#'
#' Methods whose CV predictions have zero variance (e.g., when all weights are
#' zero) are excluded from the optimization and assigned \eqn{\zeta_k = 0}.
#'
#' @seealso \code{\link{twas_weights_cv}}, \code{\link{twas_weights}},
#'   \code{\link{twas_weights_pipeline}}
#'
#' @examples
#' \dontrun{
#' # After running twas_weights_pipeline with CV:
#' res <- twas_weights_pipeline(X, y, cv_folds = 5, weight_methods = methods)
#'
#' ens <- ensemble_weights(
#'   cv_results = res$twas_cv_result,
#'   Y = y,
#'   twas_weight_list = res$twas_weights
#' )
#' ens$method_coef           # combination weights, sum to 1
#' ens$ensemble_performance  # ensemble R-squared
#'
#' # Multi-dataset ensemble (e.g., CUMC1 + MIT cell types):
#' ens_multi <- ensemble_weights(
#'   cv_results = list(res_cumc$twas_cv_result, res_mit$twas_cv_result),
#'   Y = list(y_cumc, y_mit),
#'   twas_weight_list = list(res_cumc$twas_weights, res_mit$twas_weights)
#' )
#' }
#'
#' @export
ensemble_weights <- function(cv_results, Y, twas_weight_list = NULL,
                             context_index = 1) {
  # --- Input validation ---
  if (is.null(cv_results)) {
    stop("'cv_results' is required.")
  }
  if (is.null(Y)) {
    stop("'Y' is required.")
  }
  if (!is.numeric(context_index) || length(context_index) != 1 || context_index < 1) {
    stop("'context_index' must be a positive integer scalar.")
  }

  # --- Normalize single vs multi-dataset input ---
  # Single dataset: cv_results has $prediction directly (is a twas_weights_cv() output).
  # Multi-dataset: cv_results is a list of such outputs.
  is_single <- !is.null(cv_results$prediction)
  if (is_single) {
    cv_results <- list(cv_results)
    Y <- list(Y)
    if (!is.null(twas_weight_list)) twas_weight_list <- list(twas_weight_list)
  } else {
    # Multi-dataset: validate list consistency
    if (!is.list(cv_results) || length(cv_results) == 0) {
      stop("For multi-dataset ensemble, 'cv_results' must be a non-empty list of ",
           "twas_weights_cv() outputs.")
    }
    if (!is.list(Y) || length(Y) != length(cv_results)) {
      stop("'Y' must be a list of the same length as 'cv_results' for ",
           "multi-dataset ensemble.")
    }
    if (!is.null(twas_weight_list)) {
      if (!is.list(twas_weight_list) || length(twas_weight_list) != length(cv_results)) {
        stop("'twas_weight_list' must be a list of the same length as 'cv_results'.")
      }
    }
    for (d in seq_along(cv_results)) {
      if (is.null(cv_results[[d]]$prediction)) {
        stop("cv_results[[", d, "]] does not contain '$prediction'. ",
             "Expected a twas_weights_cv() output.")
      }
    }
  }

  # --- Extract and validate method names ---
  pred_names <- names(cv_results[[1]]$prediction)
  if (is.null(pred_names) || any(pred_names == "")) {
    stop("cv_results$prediction must be a named list (output of twas_weights_cv).")
  }
  base_names <- gsub("_predicted$", "", pred_names)
  K <- length(base_names)

  if (K < 2) {
    stop("Ensemble learning requires at least 2 methods. Found: ", K, ".")
  }

  # Consistency: all datasets must report the same methods in the same order
  for (d in seq_along(cv_results)) {
    if (!identical(names(cv_results[[d]]$prediction), pred_names)) {
      stop("All cv_results must have the same method names (in $prediction) ",
           "in the same order. Dataset 1 has: ", paste(pred_names, collapse = ", "),
           "; dataset ", d, " has: ",
           paste(names(cv_results[[d]]$prediction), collapse = ", "))
    }
  }

  # --- Build stacked prediction matrix P and observed y vector ---
  pred_list <- list()
  y_list <- list()

  for (d in seq_along(cv_results)) {
    preds_d <- cv_results[[d]]$prediction
    y_raw <- Y[[d]]
    y_d <- if (is.matrix(y_raw) || is.data.frame(y_raw)) {
      if (context_index > ncol(y_raw)) {
        stop("context_index (", context_index, ") exceeds number of columns in Y[[",
             d, "]] (", ncol(y_raw), ").")
      }
      as.numeric(as.matrix(y_raw)[, context_index])
    } else {
      as.numeric(y_raw)
    }
    n_d <- length(y_d)

    P_d <- matrix(NA_real_, nrow = n_d, ncol = K)
    colnames(P_d) <- base_names
    for (k in seq_along(pred_names)) {
      pred_mat <- preds_d[[pred_names[k]]]
      p_col <- if (is.matrix(pred_mat)) pred_mat[, context_index] else as.numeric(pred_mat)
      if (length(p_col) != n_d) {
        stop("Prediction length for method '", pred_names[k], "' in dataset ", d,
             " (", length(p_col), ") does not match Y length (", n_d, ").")
      }
      P_d[, k] <- p_col
    }
    pred_list[[d]] <- P_d
    y_list[[d]] <- y_d
  }

  P <- do.call(rbind, pred_list)   # (n_total x K)
  y_obs <- unlist(y_list)           # (n_total)

  # Remove rows with any NA (in P or y)
  complete <- stats::complete.cases(P, y_obs)
  n_dropped <- sum(!complete)
  if (n_dropped > 0) {
    message("Dropping ", n_dropped, " observation(s) with NA predictions or outcomes.")
  }
  if (sum(complete) < K + 1) {
    stop("Too few complete observations (", sum(complete), ") for ", K,
         " methods. Need at least ", K + 1, ".")
  }
  P <- P[complete, , drop = FALSE]
  y_obs <- y_obs[complete]

  # --- Identify methods with non-zero variance predictions ---
  method_sds <- apply(P, 2, stats::sd)
  valid_methods <- method_sds > .Machine$double.eps
  n_valid <- sum(valid_methods)

  if (n_valid < 1) {
    stop("All methods have zero-variance predictions. Cannot compute ensemble. ",
         "This typically means all methods returned zero weights — check that ",
         "the input data has sufficient signal.")
  }

  # --- Solve the constrained QP ---
  if (n_valid == 1) {
    # Only one method has signal: assign it full weight
    zeta <- rep(0, K)
    zeta[valid_methods] <- 1
    names(zeta) <- base_names
    message("Only one method ('", base_names[valid_methods],
            "') has non-zero variance predictions. Assigning it full weight.")
  } else {
    if (!requireNamespace("quadprog", quietly = TRUE)) {
      stop("Package 'quadprog' is required for ensemble_weights. ",
           "Install with: install.packages('quadprog')")
    }

    P_valid <- P[, valid_methods, drop = FALSE]
    K_valid <- ncol(P_valid)

    Dmat <- crossprod(P_valid)
    dvec <- as.vector(crossprod(P_valid, y_obs))
    # Ridge term for numerical stability (small relative to trace)
    Dmat <- Dmat + 1e-8 * mean(diag(Dmat)) * diag(K_valid)

    # Constraint matrix: first constraint is equality (sum = 1), then K_valid
    # non-negativity constraints.
    Amat <- cbind(rep(1, K_valid), diag(K_valid))
    bvec <- c(1, rep(0, K_valid))

    qp_sol <- tryCatch(
      quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1),
      error = function(e) {
        warning("QP solver failed: ", conditionMessage(e),
                ". Falling back to equal weights among valid methods.")
        NULL
      }
    )

    if (is.null(qp_sol)) {
      zeta_valid <- rep(1 / K_valid, K_valid)
    } else {
      # Numerical cleanup: clamp to non-negative and renormalize
      zeta_valid <- pmax(qp_sol$solution, 0)
      zeta_sum <- sum(zeta_valid)
      if (zeta_sum <= 0) {
        warning("QP returned all-zero solution. Falling back to equal weights.")
        zeta_valid <- rep(1 / K_valid, K_valid)
      } else {
        zeta_valid <- zeta_valid / zeta_sum
      }
    }

    zeta <- rep(0, K)
    zeta[valid_methods] <- zeta_valid
    names(zeta) <- base_names
  }

  # --- Performance metrics ---
  ensemble_pred <- as.numeric(P %*% zeta)
  ensemble_corr <- if (stats::sd(ensemble_pred) > 0) {
    stats::cor(y_obs, ensemble_pred)
  } else NA_real_
  ensemble_rsq <- if (!is.na(ensemble_corr)) ensemble_corr^2 else NA_real_

  method_rsq <- setNames(vapply(seq_len(K), function(k) {
    if (method_sds[k] > 0) stats::cor(y_obs, P[, k])^2 else NA_real_
  }, numeric(1)), base_names)

  # --- Build ensemble TWAS weight vector (uses first dataset's weights) ---
  ensemble_twas_wt <- NULL
  if (!is.null(twas_weight_list)) {
    wt_list <- twas_weight_list[[1]]
    if (!is.list(wt_list) || length(wt_list) == 0) {
      warning("twas_weight_list[[1]] is empty or not a list; skipping weight combination.")
    } else {
      wt_keys <- paste0(base_names, "_weights")
      matched <- wt_keys %in% names(wt_list)

      if (any(matched)) {
        first_wt <- wt_list[[wt_keys[which(matched)[1]]]]
        if (!is.matrix(first_wt)) first_wt <- matrix(first_wt, ncol = 1)
        p <- nrow(first_wt)
        n_contexts <- ncol(first_wt)

        ensemble_twas_wt <- matrix(0, nrow = p, ncol = n_contexts)
        rownames(ensemble_twas_wt) <- rownames(first_wt)
        colnames(ensemble_twas_wt) <- colnames(first_wt)

        for (i in which(matched)) {
          w_mat <- wt_list[[wt_keys[i]]]
          if (!is.matrix(w_mat)) w_mat <- matrix(w_mat, ncol = 1)
          if (!identical(dim(w_mat), dim(ensemble_twas_wt))) {
            warning("Weight matrix for '", wt_keys[i],
                    "' has inconsistent dimensions; skipping.")
            next
          }
          ensemble_twas_wt <- ensemble_twas_wt + zeta[i] * w_mat
        }

        # For univariate case, return as vector
        if (n_contexts == 1) {
          ensemble_twas_wt <- setNames(
            as.numeric(ensemble_twas_wt),
            rownames(ensemble_twas_wt)
          )
        }
      } else {
        warning("No matching weight keys found in twas_weight_list. ",
                "Expected keys like: ",
                paste(wt_keys[seq_len(min(3, K))], collapse = ", "))
      }
    }
  }

  list(
    method_coef = zeta,
    ensemble_twas_weights = ensemble_twas_wt,
    method_performance = method_rsq,
    ensemble_performance = c(corr = ensemble_corr, rsq = ensemble_rsq)
  )
}
