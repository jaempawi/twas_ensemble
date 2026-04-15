# Ensemble TWAS Weights via Stacked Regression

Minimal extension to [pecotmr](https://github.com/StatFunGen/pecotmr) for learning optimal ensemble combination weights across multiple TWAS prediction methods and datasets.

## Why

pecotmr offers 10+ TWAS weight methods. Each is best for some genetic architectures but worse for others. Currently users face two bad options:

1. **Pick one method** — arbitrary, sub-optimal for most genes
2. **Run several and test all** — pay a multiple-testing penalty that hurts power

Ensemble learning is a statistically principled answer: produce **one weight vector per gene** by learning the optimal convex combination of methods, and run a single downstream test.

## What this does

Given TWAS weight predictions from K methods (e.g., SuSiE, LASSO, Elastic Net, mr.ash, BayesR, ...), learn non-negative combination coefficients that sum to 1, maximizing cross-validated prediction R². This follows the [SR-TWAS](https://www.nature.com/articles/s41467-024-50983-w) stacked regression approach.

**Key features:**
- Uses 11 prediction methods by default
- Learns optimal method combination via constrained QP (quadprog)
- Preserves per-method R² for comparison
- Supports multi-dataset ensemble (e.g., CUMC1 + MIT → one combined weight)
- Minimum changes to pecotmr (~100 lines new, 1 line modified)

## Repository structure

```
.claude/skills/twas-ensemble/SKILL.md    # Claude Code skill (/twas-ensemble)
R/ensemble_weights.R                      # The ensemble_weights() function
tests/testthat/test_ensemble_weights.R    # Unit tests (30 tests, pecotmr style)
scripts/test_ensemble_twas.R              # Simulation test using simxQTL
patches/pecotmr_ensemble.patch            # Patch for pecotmr integration
```

## Status

- Patch applies cleanly to `StatFunGen/pecotmr@main`
- R syntax check: passes
- Unit tests: 30/30 pass (1 skipped without glmnet)
- Sanity check (3 scenarios): coefficients non-negative and sum to 1,
  best method gets most weight, multi-dataset combination verified,
  weight combination numerically correct

## Quick start

### 1. Install dependencies

```r
remotes::install_github("StatFunGen/pecotmr")
remotes::install_github("StatFunGen/simxQTL")
install.packages(c("glmnet", "ncvreg", "L0Learn", "quadprog", "qgg"))
```

### 2. Single-dataset ensemble

```r
library(pecotmr)
source("R/ensemble_weights.R")  # or apply the patch to pecotmr

# After running twas_weights_pipeline with CV:
res <- twas_weights_pipeline(X, y, cv_folds = 5, weight_methods = methods)

# Learn ensemble weights
ens <- ensemble_weights(
  cv_results = res$twas_cv_result,
  Y = y,
  twas_weight_list = res$twas_weights
)
ens$method_coef           # combination weights (sum to 1)
ens$method_performance    # per-method R²
ens$ensemble_performance  # ensemble R² and correlation
ens$ensemble_twas_weights # final combined weight vector
```

### 3. Multi-dataset ensemble

```r
# Run pipeline on each dataset
res1 <- twas_weights_pipeline(X1, y1, cv_folds = 5, weight_methods = methods)
res2 <- twas_weights_pipeline(X2, y2, cv_folds = 5, weight_methods = methods)

# Learn joint weights across datasets
ens <- ensemble_weights(
  cv_results = list(res1$twas_cv_result, res2$twas_cv_result),
  Y = list(y1, y2),
  twas_weight_list = list(res1$twas_weights, res2$twas_weights)
)
```

### 4. Claude Code skill

From this repo directory, use the `/twas-ensemble` skill:
```
/twas-ensemble 500 200
```

## 11 Default Methods

| # | Method | pecotmr function | Type |
|---|--------|-------------------|------|
| 1 | SuSiE | `susie_weights` | Bayesian fine-mapping |
| 2 | LASSO | `lasso_weights` | L1 penalty |
| 3 | Elastic Net | `enet_weights` | L1 + L2 penalty |
| 4 | mr.ash | `mrash_weights` | Adaptive shrinkage |
| 5 | BayesR | `bayes_r_weights` | 4-component mixture |
| 6 | BayesL | `bayes_l_weights` | Laplace prior |
| 7 | BayesA | `bayes_a_weights` | t-distribution prior |
| 8 | BayesC | `bayes_c_weights` | Spike-and-slab |
| 9 | SCAD | `scad_weights` | Non-convex penalty |
| 10 | MCP | `mcp_weights` | Minimax concave |
| 11 | L0Learn | `l0learn_weights` | L0 regularization |

## References

- SR-TWAS: Dai et al. (2024). Nature Communications. [doi:10.1038/s41467-024-50983-w](https://www.nature.com/articles/s41467-024-50983-w)
- pecotmr: [github.com/StatFunGen/pecotmr](https://github.com/StatFunGen/pecotmr)
- simxQTL: [github.com/StatFunGen/simxQTL](https://github.com/StatFunGen/simxQTL)
- SuSiE: [stephenslab.github.io/susieR](https://stephenslab.github.io/susieR/articles/finemapping.html)
