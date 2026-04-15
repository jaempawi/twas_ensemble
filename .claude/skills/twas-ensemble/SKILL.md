---
name: twas-ensemble
description: >
  Run ensemble TWAS weight computation using pecotmr and simxQTL.
  Guides through simulation setup, 11-method weight training, stacked
  regression ensemble learning (non-negative weights summing to 1),
  cross-validation, and multi-dataset combination.
when_to_use: >
  When computing TWAS weights, comparing prediction methods, training
  ensemble weights, combining weights across datasets or cell types,
  or debugging TWAS weight pipelines in pecotmr.
argument-hint: [n_samples] [n_variants]
allowed-tools: "Bash(Rscript *) Bash(R *) Read Write Edit Grep"
---

# Ensemble TWAS Weight Computation

Guide the user through computing TWAS weights using 11 methods and learning
optimal ensemble combination weights via stacked regression (SR-TWAS approach).

## Overview

The workflow has 3 phases:
1. **Simulate or load data** (X genotypes, Y phenotypes)
2. **Run `twas_weights_pipeline()` with `ensemble=TRUE`** (11 methods + stacked regression)
3. **Inspect results** (per-method R², ensemble R², combination coefficients)

## Package Requirements

```r
# Core packages
remotes::install_github("StatFunGen/pecotmr")
remotes::install_github("StatFunGen/simxQTL")
# Suggested (for full method coverage)
install.packages(c("glmnet", "ncvreg", "L0Learn", "quadprog", "qgg"))
remotes::install_github("stephenslab/susieR")
```

## Phase 1: Data Setup

### Option A: Simulate with simxQTL

Use `$ARGUMENTS` for customization. Default: n=500 samples, p=200 variants.

```r
library(simxQTL)
set.seed(42)

n_samples <- $0   # first argument, or 500
n_variants <- $1  # second argument, or 200

# Simulate genotypes (independent or with LD)
G <- sim_geno_indep(n = n_samples, p = n_variants, min_maf = 0.05, max_maf = 0.4)
# Or with LD: G <- sim_geno_LD(n = n_samples, LD = your_LD_matrix)

# Simulate phenotype with realistic architecture
# (sparse + oligogenic + infinitesimal effects)
sim <- generate_cis_qtl_data(G, h2g = 0.20,
                              n_sparse = 3, n_oligogenic = 5, n_inf = 15)
X <- scale(sim$G)
y <- sim$y
```

### Option B: Load real data

```r
# X: samples x variants genotype matrix (scaled)
# y: phenotype vector matching samples in X
X <- readRDS("your_genotype_matrix.rds")
y <- readRDS("your_phenotype_vector.rds")
```

## Phase 2: Run Ensemble Pipeline

### 11 Default Methods

| # | Method | Function | Type |
|---|--------|----------|------|
| 1 | SuSiE | `susie_weights` | Fine-mapping (Bayesian) |
| 2 | LASSO | `lasso_weights` | Penalized regression |
| 3 | Elastic Net | `enet_weights` | Penalized regression |
| 4 | mr.ash | `mrash_weights` | Bayesian adaptive shrinkage |
| 5 | BayesR | `bayes_r_weights` | 4-component Gaussian mixture |
| 6 | BayesL | `bayes_l_weights` | Bayesian LASSO (Laplace prior) |
| 7 | BayesA | `bayes_a_weights` | t-distribution prior |
| 8 | BayesC | `bayes_c_weights` | Spike-and-slab prior |
| 9 | SCAD | `scad_weights` | Non-convex penalty |
| 10 | MCP | `mcp_weights` | Minimax concave penalty |
| 11 | L0Learn | `l0learn_weights` | L0-regularized regression |

### Run the pipeline

```r
library(pecotmr)

weight_methods <- list(
  susie_weights   = list(refine = FALSE, init_L = 5, max_L = 10),
  lasso_weights   = list(),
  enet_weights    = list(),
  mrash_weights   = list(init_prior_sd = TRUE, max.iter = 100),
  bayes_r_weights = list(),
  bayes_l_weights = list(),
  bayes_a_weights = list(),
  bayes_c_weights = list(),
  scad_weights    = list(),
  mcp_weights     = list(),
  l0learn_weights = list()
)

res <- twas_weights_pipeline(X, y,
  cv_folds = 5,
  weight_methods = weight_methods,
  ensemble = TRUE        # enables stacked regression
)
```

### Inspect results

```r
# Combination coefficients (zeta_k >= 0, sum = 1)
print(round(res$ensemble$method_coef, 4))

# Per-method R-squared (kept for each method)
print(round(res$ensemble$method_performance, 4))

# Ensemble R-squared and correlation
print(res$ensemble$ensemble_performance)

# Final combined TWAS weight vector
str(res$ensemble$ensemble_twas_weights)
```

## Phase 3: Multi-Dataset Ensemble

Combine across datasets or cell types (e.g., CUMC1 + MIT) to learn joint weights:

```r
# Run pipeline on each dataset separately
res_cumc <- twas_weights_pipeline(X_cumc, y_cumc, cv_folds = 5,
                                   weight_methods = weight_methods)
res_mit  <- twas_weights_pipeline(X_mit, y_mit, cv_folds = 5,
                                   weight_methods = weight_methods)

# Combine via ensemble (stacks predictions from both datasets)
ens <- ensemble_weights(
  cv_results = list(res_cumc$twas_cv_result, res_mit$twas_cv_result),
  Y = list(y_cumc, y_mit),
  twas_weight_list = list(res_cumc$twas_weights, res_mit$twas_weights)
)

# Joint learned weights
print(ens$method_coef)
print(ens$ensemble_performance)
```

**Key point:** This produces ONE set of ensemble weights from both datasets,
so you don't need to keep two separate weight vectors.

## Memory Management

When running 11 methods on large datasets, memory can be an issue. Strategies:

### 1. Sequential execution (default)
Use `cv_threads = 1` (the default) so methods run one at a time, freeing
memory between each.

### 2. Limit variants during CV
```r
res <- twas_weights_pipeline(X, y,
  cv_folds = 5,
  weight_methods = weight_methods,
  max_cv_variants = 500,   # subsample variants for CV speed
  ensemble = TRUE
)
```

### 3. Split-and-merge for extreme cases
If 11 methods don't fit in memory together, split into batches:

```r
methods_batch1 <- weight_methods[1:6]
methods_batch2 <- weight_methods[7:11]

cv1 <- twas_weights_cv(X, y, fold = 5, weight_methods = methods_batch1)
gc()  # free memory
cv2 <- twas_weights_cv(X, y, fold = 5, weight_methods = methods_batch2)

# Merge predictions and call ensemble_weights
merged_cv <- cv1
merged_cv$prediction <- c(cv1$prediction, cv2$prediction)
merged_cv$performance <- c(cv1$performance, cv2$performance)

# Compute full-data weights separately too
wt1 <- twas_weights(X, y, weight_methods = methods_batch1)
wt2 <- twas_weights(X, y, weight_methods = methods_batch2)
merged_wt <- c(wt1, wt2)

ens <- ensemble_weights(merged_cv, y, twas_weight_list = merged_wt)
```

## Troubleshooting

### Low ensemble R² (< 0.01)
- Increase heritability in simulation (`h2g`)
- Try more samples (`n_samples >= 300`)
- Check that variants have sufficient MAF variation

### Method failures
- Some methods require specific packages (e.g., `qgg` for BayesR)
- Failed methods return zero weights and are excluded from CV
- Check warnings: "Methods removed from CV because all weights are zeros"

### QP solver failure
- Falls back to equal weights automatically
- Usually caused by highly collinear predictions
- Try removing redundant methods (e.g., keep only one of LASSO/Elastic Net)

### Memory errors
- Use `max_cv_variants` to reduce variant count
- Use split-and-merge approach above
- Set `cv_threads = 1` (avoid parallel overhead)

## Algorithm Reference

The ensemble solves the constrained quadratic program:

```
min_zeta ||y - P * zeta||^2
subject to: zeta_k >= 0 for all k, sum(zeta) = 1
```

where P is the n x K matrix of out-of-fold CV predictions from K methods.
This is the SR-TWAS stacked regression approach
(Dai et al., Nature Communications, 2024, doi:10.1038/s41467-024-50983-w).

The final TWAS weight vector is: `w_ensemble = sum_k(zeta_k * w_k)`
where w_k is the full-data weight vector from method k.
