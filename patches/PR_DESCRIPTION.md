# Add ensemble TWAS weights via stacked regression

## Motivation

pecotmr currently offers 10+ TWAS weight methods. Each is best for some
genetic architectures (sparse, oligogenic, polygenic) but worse for others.
This leaves users with two unattractive options:

1. **Pick one method** — an arbitrary choice that's sub-optimal for most genes.
2. **Run several and test all** — pay a multiple-testing penalty that hurts
   power, and produces a conceptually awkward "which method is best?" story.

Ensemble learning offers a statistically principled alternative:
**one weight vector per gene**, learned as a convex combination of methods,
enabling a single downstream TWAS test per gene with no method-selection bias.

## What this PR adds

1. **`ensemble_weights()`** — new function in `R/twas_weights.R`. Given
   out-of-fold CV predictions from K methods, solves a constrained quadratic
   program to learn coefficients ζ₁,…,ζₖ with ζₖ ≥ 0 and Σζₖ = 1 that
   maximize prediction R² (stacked regression). Returns the combined weight
   vector along with per-method R² for reporting.

2. **`ensemble` parameter in `twas_weights_pipeline()`** — opt-in
   (`ensemble = FALSE` by default). Fully backward-compatible.

3. **Multi-dataset ensemble** — `ensemble_weights()` accepts a list of
   `twas_weights_cv()` outputs, enabling joint weight learning across
   reference panels or cell types (one set of coefficients for multiple
   datasets, instead of per-dataset weights).

4. **Tests** — `tests/testthat/test_ensemble_weights.R` with 30 tests
   covering input validation, algorithm correctness, edge cases
   (zero-variance methods, NA predictions), weight combination, and
   multi-dataset.

## Reference

Dai, Fang, et al. (2024). **SR-TWAS: leveraging multiple reference panels
to improve transcriptome-wide association study power by ensemble machine
learning.** *Nature Communications* 15, 6486.
[doi:10.1038/s41467-024-50983-w](https://www.nature.com/articles/s41467-024-50983-w)

## API

```r
# Opt-in via pipeline parameter:
res <- twas_weights_pipeline(X, y, cv_folds = 5,
                              weight_methods = methods,
                              ensemble = TRUE)
res$ensemble$method_coef           # zeta_k, sum to 1, >= 0
res$ensemble$ensemble_twas_weights # combined weight vector
res$ensemble$method_performance    # per-method R^2 (preserved)
res$ensemble$ensemble_performance  # ensemble R^2 and corr

# Direct call (also supports multi-dataset):
ens <- ensemble_weights(
  cv_results = res$twas_cv_result,
  Y = y,
  twas_weight_list = res$twas_weights
)

# Multi-dataset (e.g., CUMC1 + MIT cell types):
ens <- ensemble_weights(
  cv_results = list(res_cumc$twas_cv_result, res_mit$twas_cv_result),
  Y = list(y_cumc, y_mit),
  twas_weight_list = list(res_cumc$twas_weights, res_mit$twas_weights)
)
```

## Backward compatibility

- `ensemble = FALSE` is the default; no change in behavior for existing users.
- No new hard dependencies — `quadprog` is already in `Suggests`.
- No changes to `NAMESPACE` beyond an auto-generated `@export` entry.
- No changes to `DESCRIPTION`.

## Verification

- [x] Patch applies cleanly to `StatFunGen/pecotmr@main`
- [x] R syntax check passes on patched `R/twas_weights.R`
- [x] All 30 unit tests pass (1 skipped when `glmnet` unavailable)
- [x] Sanity checks:
  - Coefficients non-negative and sum to 1
  - Best method receives largest coefficient (asymptotically)
  - `ensemble_twas_weights` equals Σₖ ζₖ·wₖ to 1e-10 precision
  - Multi-dataset ensemble produces valid joint coefficients

## Out of scope (future PRs)

- LD sketch integration for summary-statistics methods (block-boundary gap fix)
- Ensemble-as-default configuration in `univariate_analysis_pipeline`
