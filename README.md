# `ebal` : Entropy Balancing

Updated version of the `ebal` R package for the entropy balancing method to compute method-of-moments based inverse propensity scores for causal inference and related problems.

Primary function is `ebalance(W, X)` to estimate the ATT (where `W` is the treatment vector).

__New features__: The main user-facing function `ebalance` optionally computes weights using automatic differentiation (via `torch`) when called with `method="AutoDiff"`. This greatly improves stability (at the cost of marginally slower runtime) for small problems and greatly improves speed and stability for large problems. For detailed performance comparisons, see [here](https://nbviewer.org/github/apoorvalal/ebal/blob/master/vignettes/performance_comparisons.pdf).

## Usage

```r
library(ebal)
# lalonde observational data
load(url("https://github.com/apoorvalal/LalRUtils/raw/master/data/lalonde.psid.RData"))
yn = "re78"; wn = "treat"; xn = setdiff(colnames(lalonde.psid), c(yn, wn))

################################################################################
# naive estimate - this is badly biased relative to experimental benchmark of 1794
with(lalonde.psid, mean(re78[treat == 1]) - mean(re78[treat == 0]))
# -15204.7755
```

`ebalance` takes a treatment vector and covariates and returns a weight vector for control units.

```r
# solving for weights
bal_wts = ebalance(lalonde.psid[['treat']], lalonde.psid[,xn], method = "AutoDiff")$w
```

which can then be fed into `weighted.mean`.

```r
################################################################################
# balancing estimate - better
with(lalonde.psid, mean(re78[treat == 1]) - weighted.mean(re78[treat == 0], bal_wts))
# 2422.7965
```

or into a regression.

```r
################################################################################
# augmented balancing estimate

# create weights vector in dataset for regression - uniform weights for treated units
lalonde.psid[['wt']] = 1/sum(lalonde.psid$treat)            # treat weights = 1/Nt
lalonde.psid[['wt']][lalonde.psid$treat == 0] = bal_wts     # ctrl weights = ebal weights

estimatr::lm_robust(re78 ~ . - wt, weights = wt, lalonde.psid) %>%
  summary %>% coefficients %>% .['treat',] %>% round(2)

# Estimate Std. Error    t value   Pr(>|t|)   CI Lower   CI Upper         DF
#  2422.79     750.36       3.23       0.00     953.00    3895.71    2663.00
```


## installation

```
remotes::install_github("apoorvalal/ebal")
```

## lower-level functions for usage with aggregate data, and other reweighting problems
Reweighting may be useful in a wide variety of settings (e.g. survey reweighting, experiment generalization, density ratio estimation). In many cases, you may only have summary statistics and not individual level covariates for the target sample.

In such cases, you may use `eb` or `ebal_torch` functions directly, where you supply a donor matrix `X0`, a target moment vector `X1`, and optionally supply base weights.

