#' Ensemble TWAS Weights via Stacked Regression
#'
#' Given cross-validated predictions from multiple methods (output of \code{twas_weights_cv})
#' and optionally from multiple datasets, learns non-negative combination weights
#' (summing to 1) via constrained least squares (stacked regression). Returns ensemble
#' weights and per-method performance metrics.
#'
#' For single-dataset usage, pass one \code{twas_weights_cv()} result directly.
#' For multi-dataset ensemble (e.g., combining CUMC1 + MIT cell types), pass a list
#' of \code{twas_weights_cv()} results along with a list of observed Y vectors.
#'
#' @param cv_results Output of \code{twas_weights_cv()}, containing \code{$prediction}
#'   (named list of method -> matrix with keys like \code{"susie_predicted"}) and
#'   \code{$performance}. For multi-dataset: a list of such objects.
#' @param Y Observed outcome vector or matrix. For multi-dataset: a list of vectors/matrices.
#' @param twas_weight_list Optional named list of weight matrices from \code{twas_weights()},
#'   with keys like \code{"susie_weights"}. Used to construct the final ensemble weight vector.
#'   For multi-dataset: a list of such lists (uses the first dataset's weights as base).
#' @param context_index Integer, which column of Y to use (default 1, for univariate).
#'
#' @return A list with:
#' \describe{
#'   \item{method_coef}{Named numeric vector of combination coefficients (zeta_k),
#'     non-negative and summing to 1.}
#'   \item{ensemble_twas_weights}{Final combined weight vector (or NULL if
#'     \code{twas_weight_list} not provided).}
#'   \item{method_performance}{Named numeric vector of per-method R-squared from CV.}
#'   \item{ensemble_performance}{Named vector with \code{corr} and \code{rsq} for
#'     the ensemble.}
#' }
#'
#' @details
#' The stacked regression solves:
#' \deqn{\min_{\zeta} \|y - P\zeta\|^2 \quad \text{s.t.} \quad \zeta \geq 0, \sum \zeta = 1}
#' where P is the n x K matrix of out-of-fold predictions from K methods.
#' This follows the SR-TWAS approach (Dai et al., Nature Communications, 2024).
#'
#' @export
ensemble_weights <- function(cv_results, Y, twas_weight_list = NULL,
                              context_index = 1) {
  # --- Handle single vs multi-dataset input ---
  # Single dataset: cv_results has $prediction directly
  # Multi-dataset: cv_results is a list of cv_result objects
  is_single <- !is.null(cv_results$prediction)
  if (is_single) {
    cv_results <- list(cv_results)
    Y <- list(Y)
    if (!is.null(twas_weight_list)) twas_weight_list <- list(twas_weight_list)
  }

  # --- Extract method names from prediction keys ---
  pred_names <- names(cv_results[[1]]$prediction)
  # Keys are like "susie_predicted", "enet_predicted"
  base_names <- gsub("_predicted$", "", pred_names)
  K <- length(base_names)

  if (K < 2) {
    stop("Ensemble learning requires at least 2 methods. Found: ", K)
  }

  # --- Build stacked prediction matrix P and observed y vector ---
  pred_list <- list()
  y_list <- list()

  for (d in seq_along(cv_results)) {
    preds_d <- cv_results[[d]]$prediction
    y_d <- if (is.matrix(Y[[d]])) Y[[d]][, context_index] else Y[[d]]
    n_d <- length(y_d)

    P_d <- matrix(NA, nrow = n_d, ncol = K)
    colnames(P_d) <- base_names
    for (k in seq_along(pred_names)) {
      pred_mat <- preds_d[[pred_names[k]]]
      P_d[, k] <- if (is.matrix(pred_mat)) pred_mat[, context_index] else pred_mat
    }
    pred_list[[d]] <- P_d
    y_list[[d]] <- y_d
  }

  P <- do.call(rbind, pred_list)   # (n_total x K)
  y_obs <- unlist(y_list)           # (n_total)

  # Remove rows with NA
  complete <- complete.cases(P, y_obs)
  if (sum(complete) < K + 1) {
    stop("Too few complete observations (", sum(complete), ") for ", K, " methods.")
  }
  P <- P[complete, , drop = FALSE]
  y_obs <- y_obs[complete]

  # --- Remove methods with zero-variance predictions ---
  method_sds <- apply(P, 2, sd)
  valid_methods <- method_sds > 0
  if (sum(valid_methods) < 1) {
    stop("All methods have zero-variance predictions. Cannot compute ensemble.")
  }

  # If only 1 valid method, assign it weight 1
  if (sum(valid_methods) == 1) {
    zeta <- ifelse(valid_methods, 1, 0)
    names(zeta) <- base_names
    message("Only one method has non-zero variance predictions. Assigning it full weight.")
  } else {
    # --- Solve constrained QP ---
    # min 0.5 * zeta' * Dmat * zeta - dvec' * zeta
    # s.t. Amat' * zeta >= bvec, with first meq constraints as equalities
    if (!requireNamespace("quadprog", quietly = TRUE)) {
      stop("Package 'quadprog' is required for ensemble_weights. ",
           "Install with: install.packages('quadprog')")
    }

    P_valid <- P[, valid_methods, drop = FALSE]
    K_valid <- ncol(P_valid)

    Dmat <- crossprod(P_valid)                    # K_valid x K_valid
    dvec <- as.vector(crossprod(P_valid, y_obs))  # K_valid

    # Ridge for numerical stability
    Dmat <- Dmat + 1e-8 * diag(K_valid)

    # Constraints: first meq=1 is equality (sum = 1), then K_valid inequalities (>= 0)
    Amat <- cbind(rep(1, K_valid), diag(K_valid))
    bvec <- c(1, rep(0, K_valid))

    qp_sol <- tryCatch(
      quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1),
      error = function(e) {
        message("QP solver failed: ", e$message, ". Falling back to equal weights.")
        NULL
      }
    )

    if (is.null(qp_sol)) {
      zeta_valid <- rep(1 / K_valid, K_valid)
    } else {
      zeta_valid <- pmax(qp_sol$solution, 0)  # numerical cleanup
      zeta_valid <- zeta_valid / sum(zeta_valid)  # re-normalize
    }

    # Map back to full method vector
    zeta <- rep(0, K)
    zeta[valid_methods] <- zeta_valid
    names(zeta) <- base_names
  }

  # --- Compute performance metrics ---
  # Ensemble prediction
  ensemble_pred <- P %*% zeta
  ensemble_corr <- if (sd(ensemble_pred) > 0) cor(y_obs, ensemble_pred) else NA
  ensemble_rsq <- if (!is.na(ensemble_corr)) ensemble_corr^2 else NA

  # Per-method R²
  method_rsq <- setNames(sapply(1:K, function(k) {
    if (sd(P[, k]) > 0) cor(y_obs, P[, k])^2 else NA
  }), base_names)

  # --- Build ensemble TWAS weight vector ---
  ensemble_twas_wt <- NULL
  if (!is.null(twas_weight_list)) {
    # Use first dataset's weights as the base
    wt_list <- twas_weight_list[[1]]

    # Match method names: base_names -> paste0(base_name, "_weights")
    wt_keys <- paste0(base_names, "_weights")
    matched <- wt_keys %in% names(wt_list)

    if (any(matched)) {
      first_wt <- wt_list[[which(matched)[1]]]
      p <- nrow(first_wt)
      n_contexts <- ncol(first_wt)

      ensemble_twas_wt <- matrix(0, nrow = p, ncol = n_contexts)
      rownames(ensemble_twas_wt) <- rownames(first_wt)
      colnames(ensemble_twas_wt) <- colnames(first_wt)

      for (i in which(matched)) {
        w_mat <- wt_list[[wt_keys[i]]]
        ensemble_twas_wt <- ensemble_twas_wt + zeta[i] * w_mat
      }

      # For univariate case, return as vector
      if (n_contexts == 1) {
        ensemble_twas_wt <- ensemble_twas_wt[, 1]
      }
    } else {
      warning("No matching weight keys found in twas_weight_list. ",
              "Expected keys like: ", paste(wt_keys[1:min(3, K)], collapse = ", "))
    }
  }

  return(list(
    method_coef = zeta,
    ensemble_twas_weights = ensemble_twas_wt,
    method_performance = method_rsq,
    ensemble_performance = c(corr = ensemble_corr, rsq = ensemble_rsq)
  ))
}
