context("ensemble_weights")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Build a synthetic twas_weights_cv() output with K methods. Each method's
# prediction is a convex combination of the truth + noise, letting us control
# per-method accuracy. Returns a list shaped exactly like twas_weights_cv()'s
# output (with $prediction, $performance, $sample_partition).
make_cv_result <- function(n = 100, K = 4, seed = 1, method_quality = NULL) {
  set.seed(seed)
  y <- rnorm(n)
  sample_names <- paste0("sample_", seq_len(n))

  if (is.null(method_quality)) {
    # Methods with decreasing quality (noise amounts)
    method_quality <- seq(0.1, 0.9, length.out = K)
  }
  stopifnot(length(method_quality) == K)

  method_names <- paste0("method", seq_len(K))
  pred_names <- paste0(method_names, "_predicted")

  prediction <- setNames(lapply(seq_len(K), function(k) {
    noise_sd <- method_quality[k]
    pred <- y + rnorm(n, sd = noise_sd)
    mat <- matrix(pred, ncol = 1)
    rownames(mat) <- sample_names
    colnames(mat) <- "outcome_1"
    mat
  }), pred_names)

  # Dummy performance (not used by ensemble_weights)
  performance <- setNames(lapply(seq_len(K), function(k) {
    m <- matrix(NA, nrow = 1, ncol = 6)
    colnames(m) <- c("corr", "rsq", "adj_rsq", "pval", "RMSE", "MAE")
    m
  }), paste0(method_names, "_performance"))

  list(
    sample_partition = data.frame(Sample = sample_names,
                                   Fold = rep(1:5, length.out = n),
                                   stringsAsFactors = FALSE),
    prediction = prediction,
    performance = performance,
    time_elapsed = 0,
    .y = y,
    .method_names = method_names
  )
}

# Build synthetic twas_weights() output
make_weight_list <- function(p = 20, method_names, seed = 2) {
  set.seed(seed)
  setNames(lapply(method_names, function(m) {
    w <- matrix(rnorm(p), ncol = 1)
    rownames(w) <- paste0("var_", seq_len(p))
    colnames(w) <- "outcome_1"
    w
  }), paste0(method_names, "_weights"))
}

# ===========================================================================
#  Input validation
# ===========================================================================

test_that("ensemble_weights: NULL cv_results errors", {
  expect_error(ensemble_weights(NULL, Y = rnorm(10)), "cv_results")
})

test_that("ensemble_weights: NULL Y errors", {
  cv <- make_cv_result(n = 20, K = 3)
  expect_error(ensemble_weights(cv, Y = NULL), "'Y' is required")
})

test_that("ensemble_weights: single method errors (need >= 2 for ensemble)", {
  cv <- make_cv_result(n = 20, K = 1)
  expect_error(ensemble_weights(cv, Y = cv$.y),
               "at least 2 methods")
})

test_that("ensemble_weights: invalid context_index errors", {
  cv <- make_cv_result(n = 20, K = 3)
  expect_error(ensemble_weights(cv, Y = cv$.y, context_index = 0),
               "context_index")
  expect_error(ensemble_weights(cv, Y = cv$.y, context_index = "a"),
               "context_index")
})

test_that("ensemble_weights: context_index beyond Y columns errors", {
  cv <- make_cv_result(n = 20, K = 3)
  Y_mat <- matrix(cv$.y, ncol = 1)
  expect_error(ensemble_weights(cv, Y = Y_mat, context_index = 5),
               "context_index")
})

test_that("ensemble_weights: multi-dataset with mismatched lengths errors", {
  cv1 <- make_cv_result(n = 20, K = 3, seed = 1)
  cv2 <- make_cv_result(n = 20, K = 3, seed = 2)
  expect_error(ensemble_weights(list(cv1, cv2), Y = list(cv1$.y)),
               "same length")
})

test_that("ensemble_weights: multi-dataset with different methods errors", {
  cv1 <- make_cv_result(n = 20, K = 3, seed = 1)
  cv2 <- make_cv_result(n = 20, K = 4, seed = 2)
  expect_error(
    ensemble_weights(list(cv1, cv2), Y = list(cv1$.y, cv2$.y)),
    "same method names"
  )
})

# ===========================================================================
#  Core algorithm correctness
# ===========================================================================

test_that("ensemble_weights: coefficients are non-negative and sum to 1", {
  cv <- make_cv_result(n = 100, K = 4, seed = 42)
  res <- ensemble_weights(cv, Y = cv$.y)

  expect_true(all(res$method_coef >= 0))
  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
})

test_that("ensemble_weights: best method receives the largest coefficient", {
  # Method 1 is best (lowest noise), method K is worst
  cv <- make_cv_result(n = 200, K = 4, seed = 7,
                        method_quality = c(0.1, 0.5, 0.8, 1.2))
  res <- ensemble_weights(cv, Y = cv$.y)

  expect_equal(names(which.max(res$method_coef)), "method1")
})

test_that("ensemble_weights: ensemble R^2 >= best single-method R^2 (or close)", {
  cv <- make_cv_result(n = 300, K = 5, seed = 13)
  res <- ensemble_weights(cv, Y = cv$.y)

  best_method_r2 <- max(res$method_performance, na.rm = TRUE)
  ensemble_r2 <- res$ensemble_performance["rsq"]

  # Stacked regression is constrained (sum=1, non-neg), so it can't always
  # beat the best single method, but should be very close.
  expect_true(ensemble_r2 >= best_method_r2 - 0.01)
})

test_that("ensemble_weights: per-method R^2 values are sensible (between 0 and 1)", {
  cv <- make_cv_result(n = 200, K = 4, seed = 21)
  res <- ensemble_weights(cv, Y = cv$.y)

  expect_true(all(res$method_performance >= 0, na.rm = TRUE))
  expect_true(all(res$method_performance <= 1, na.rm = TRUE))
  expect_equal(length(res$method_performance), 4)
})

test_that("ensemble_weights: method names are stripped of _predicted suffix", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  res <- ensemble_weights(cv, Y = cv$.y)

  expect_equal(names(res$method_coef),
               c("method1", "method2", "method3"))
  expect_equal(names(res$method_performance),
               c("method1", "method2", "method3"))
})

# ===========================================================================
#  Zero-variance / edge cases
# ===========================================================================

test_that("ensemble_weights: zero-variance method gets coefficient 0", {
  cv <- make_cv_result(n = 100, K = 3, seed = 5)
  # Force method 2 to have constant predictions
  cv$prediction$method2_predicted[, 1] <- 0.5
  res <- ensemble_weights(cv, Y = cv$.y)

  expect_equal(res$method_coef["method2"], c(method2 = 0))
  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
})

test_that("ensemble_weights: NA predictions in some samples are dropped", {
  cv <- make_cv_result(n = 100, K = 3, seed = 5)
  cv$prediction$method1_predicted[1:5, 1] <- NA
  expect_message(
    res <- ensemble_weights(cv, Y = cv$.y),
    "Dropping"
  )
  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
})

test_that("ensemble_weights: all zero-variance methods errors", {
  cv <- make_cv_result(n = 50, K = 2, seed = 5)
  cv$prediction$method1_predicted[, 1] <- 0
  cv$prediction$method2_predicted[, 1] <- 0
  expect_error(ensemble_weights(cv, Y = cv$.y),
               "zero-variance predictions")
})

# ===========================================================================
#  Weight combination
# ===========================================================================

test_that("ensemble_weights: ensemble_twas_weights is sum of zeta_k * w_k", {
  cv <- make_cv_result(n = 100, K = 3, seed = 42)
  wt <- make_weight_list(p = 10, method_names = cv$.method_names)

  res <- ensemble_weights(cv, Y = cv$.y, twas_weight_list = wt)

  expect_false(is.null(res$ensemble_twas_weights))

  # Verify the combination is correct
  expected <- matrix(0, nrow = 10, ncol = 1)
  for (k in seq_along(cv$.method_names)) {
    m <- cv$.method_names[k]
    expected <- expected + res$method_coef[m] * wt[[paste0(m, "_weights")]]
  }
  expect_equal(as.numeric(res$ensemble_twas_weights),
               as.numeric(expected),
               tolerance = 1e-10)
})

test_that("ensemble_weights: NULL twas_weight_list returns NULL ensemble_twas_weights", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  res <- ensemble_weights(cv, Y = cv$.y, twas_weight_list = NULL)
  expect_null(res$ensemble_twas_weights)
})

test_that("ensemble_weights: weights with no matching keys warns and skips", {
  cv <- make_cv_result(n = 50, K = 2, seed = 1)
  wt <- list(unknown_weights = matrix(1, nrow = 10, ncol = 1))

  expect_warning(
    res <- ensemble_weights(cv, Y = cv$.y, twas_weight_list = wt),
    "No matching weight keys"
  )
  expect_null(res$ensemble_twas_weights)
})

# ===========================================================================
#  Multi-dataset ensemble
# ===========================================================================

test_that("ensemble_weights: multi-dataset combines predictions correctly", {
  cv1 <- make_cv_result(n = 80, K = 3, seed = 1)
  cv2 <- make_cv_result(n = 80, K = 3, seed = 2)

  res <- ensemble_weights(
    cv_results = list(cv1, cv2),
    Y = list(cv1$.y, cv2$.y)
  )

  expect_true(all(res$method_coef >= 0))
  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
  expect_equal(length(res$method_performance), 3)
})

test_that("ensemble_weights: Y as matrix with context_index works", {
  cv <- make_cv_result(n = 50, K = 3, seed = 1)
  Y_mat <- matrix(cv$.y, ncol = 1)
  colnames(Y_mat) <- "ctx1"

  res <- ensemble_weights(cv, Y = Y_mat, context_index = 1)
  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
})

# ===========================================================================
#  End-to-end with twas_weights_cv (integration)
# ===========================================================================

test_that("ensemble_weights: end-to-end with twas_weights_cv output", {
  skip_if_not_installed("glmnet")

  set.seed(42)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("var_", seq_len(p))
  rownames(X) <- paste0("sample_", seq_len(n))

  beta <- c(1.5, -1.0, 0.8, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.5))

  cv <- suppressMessages(twas_weights_cv(
    X, y, fold = 3,
    weight_methods = list(
      lasso_weights = list(),
      enet_weights = list()
    )
  ))

  res <- ensemble_weights(cv, Y = y)

  expect_equal(sum(res$method_coef), 1, tolerance = 1e-6)
  expect_true(all(res$method_coef >= 0))
  expect_equal(names(res$method_coef), c("lasso", "enet"))
  expect_true(res$ensemble_performance["rsq"] > 0)
})
