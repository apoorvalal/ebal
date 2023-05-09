# `ebal` : Entropy Balancing

Updated version of the `ebal` R package for the entropy balancing method to compute method-of-moments based inverse propensity scores for causal inference and related problems.

Primary function is `ebalance(W, X)` to estimate the ATT (where `W` is the treatment vector).

__New features__: The main user-facing function `ebalance` optionally computes weights using automatic differentiation (via `torch`) when called with `method="AutoDiff"`. This greatly improves stability (at the cost of marginally slower runtime) for small problems and greatly improves speed and stability for large problems. For detailed performance comparisons, see [here](https://nbviewer.org/github/apoorvalal/ebal/blob/master/vignettes/performance_comparisons.pdf).




## installation:

```
remotes::install_github("apoorvalal/ebal")
```

## lower-level functions for usage with aggregate data, and other reweighting problems 
Reweighting may be useful in a wide variety of settings (e.g. survey reweighting, experiment generalization, density ratio estimation). In many cases, you may only have summary statistics and not individual level covariates for the target sample.

In such cases, you may use `eb` or `ebal_torch` functions directly, where you supply a donor matrix `X0`, a target moment vector `X1`, and optionally supply base weights.


