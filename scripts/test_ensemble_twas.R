#!/usr/bin/env Rscript
#
# test_ensemble_twas.R
#
# Simulation test for ensemble TWAS weights using pecotmr + simxQTL.
# Demonstrates:
#   1. Simulating X (genotype) and Y (phenotype) with simxQTL
#   2. Running 11 TWAS weight methods via twas_weights_pipeline()
#   3. Learning ensemble combination weights via stacked regression
#   4. Comparing per-method R² vs ensemble R²
#   5. Multi-dataset ensemble (combining two simulated datasets)
#
# Usage:
#   Rscript test_ensemble_twas.R
#
# Requirements:
#   remotes::install_github("StatFunGen/pecotmr")
#   remotes::install_github("StatFunGen/simxQTL")
#   install.packages(c("glmnet", "ncvreg", "L0Learn", "quadprog", "qgg"))
#

library(simxQTL)
library(pecotmr)

set.seed(42)

# ============================================================
# 1. Simulate genotype data (X)
# ============================================================
cat("=== Phase 1: Simulating genotype data ===\n")

n_samples <- 500
n_variants <- 200

G <- sim_geno_indep(n = n_samples, p = n_variants,
                     min_maf = 0.05, max_maf = 0.4)
cat(sprintf("Genotype matrix: %d samples x %d variants\n", nrow(G), ncol(G)))

# ============================================================
# 2. Simulate phenotype data (Y) with known genetic architecture
# ============================================================
cat("\n=== Phase 2: Simulating phenotype data ===\n")

# Multi-component architecture: sparse + oligogenic + infinitesimal
sim <- generate_cis_qtl_data(G, h2g = 0.20,
                              n_sparse = 3,
                              n_oligogenic = 5,
                              n_inf = 15)
X <- scale(sim$G)
y <- sim$y

cat(sprintf("True heritability: %.2f\n", sim$h2g))
cat(sprintf("Sparse causal SNPs: %d\n", length(sim$sparse_indices)))
cat(sprintf("Oligogenic SNPs: %d\n", length(sim$oligogenic_indices)))
cat(sprintf("Infinitesimal SNPs: %d\n", length(sim$infinitesimal_indices)))

# ============================================================
# 3. Define 11 TWAS weight methods
# ============================================================
cat("\n=== Phase 3: Defining 11 weight methods ===\n")

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

cat(sprintf("Methods: %s\n", paste(names(weight_methods), collapse = ", ")))

# ============================================================
# 4. Run TWAS pipeline with ensemble learning
# ============================================================
cat("\n=== Phase 4: Running TWAS pipeline with ensemble ===\n")

res <- twas_weights_pipeline(X, y,
  cv_folds = 5,
  weight_methods = weight_methods,
  ensemble = TRUE
)

# ============================================================
# 5. Report results
# ============================================================
cat("\n=== Phase 5: Results ===\n")

cat("\n--- Ensemble combination coefficients (zeta_k) ---\n")
cat("Constraint: sum = 1, all >= 0\n")
coef <- res$ensemble$method_coef
print(round(coef, 4))
cat(sprintf("Sum of coefficients: %.6f\n", sum(coef)))

cat("\n--- Per-method R-squared (from CV predictions) ---\n")
method_r2 <- res$ensemble$method_performance
print(round(method_r2, 4))

cat("\n--- Ensemble performance ---\n")
ens_perf <- res$ensemble$ensemble_performance
cat(sprintf("Ensemble correlation: %.4f\n", ens_perf["corr"]))
cat(sprintf("Ensemble R-squared:  %.4f\n", ens_perf["rsq"]))
cat(sprintf("Best single-method R-squared: %.4f (%s)\n",
            max(method_r2, na.rm = TRUE),
            names(which.max(method_r2))))

improvement <- ens_perf["rsq"] - max(method_r2, na.rm = TRUE)
cat(sprintf("Ensemble improvement: %+.4f\n", improvement))

cat("\n--- Ensemble TWAS weight vector ---\n")
ens_wt <- res$ensemble$ensemble_twas_weights
cat(sprintf("Length: %d, Non-zero: %d\n", length(ens_wt), sum(ens_wt != 0)))

# ============================================================
# 6. Multi-dataset ensemble test
# ============================================================
cat("\n=== Phase 6: Multi-dataset ensemble test ===\n")
cat("Simulating second dataset (different phenotype, same X)...\n")

sim2 <- generate_cis_qtl_data(G, h2g = 0.15,
                               n_sparse = 2,
                               n_oligogenic = 3,
                               n_inf = 10)
y2 <- sim2$y

res2 <- twas_weights_pipeline(X, y2,
  cv_folds = 5,
  weight_methods = weight_methods,
  ensemble = FALSE  # we'll do ensemble manually below
)

# Combine two datasets via ensemble_weights()
ens_multi <- ensemble_weights(
  cv_results = list(res$twas_cv_result, res2$twas_cv_result),
  Y = list(y, y2),
  twas_weight_list = list(res$twas_weights, res2$twas_weights)
)

cat("\n--- Multi-dataset ensemble coefficients ---\n")
print(round(ens_multi$method_coef, 4))
cat(sprintf("Multi-dataset ensemble R-squared: %.4f\n",
            ens_multi$ensemble_performance["rsq"]))

cat("\n=== All tests complete ===\n")
