# Ancillary Arguments for Controlling `coxph_mpl` Fits

Set numeric and algorithmic controls for `coxph_mpl` fits. The function
validates inputs (e.g., number of events per basis element, iteration
limits) to avoid impossible settings.

## Usage

``` r
coxph_mpl.control(
  n.obs = NULL,
  basis = "uniform",
  smooth = NULL,
  max.iter = c(150, 75000, 1e+06),
  tol = 1e-07,
  n.knots = NULL,
  n.events_basis = NULL,
  range.quant = c(0.075, 0.9),
  cover.sigma.quant = 0.25,
  cover.sigma.fixed = 0.25,
  min.theta = 1e-10,
  penalty = 2L,
  order = 3L,
  kappa = 1/0.6,
  epsilon = c(1e-16, 1e-10),
  ties = "epsilon",
  seed = NULL
)
```

## Arguments

- n.obs:

  Number of fully observed (non-censored) outcomes. Required when
  `basis == "uniform"` to derive an acceptable range for
  `n.events_basis`.

- basis:

  Basis used to approximate the baseline hazard. One of `"uniform"`,
  `"gaussian"`, `"msplines"`, or `"epanechikov"`. Defaults to
  `"uniform"`.

- smooth:

  Smoothing parameter value. Defaults to `NULL` (estimated via REML).
  Set to `0` for maximum-likelihood estimates. Must be non-negative.

- max.iter:

  Integer vector of length 3 giving maximum iterations for (1) smoothing
  parameter updates, (2) inner beta/theta updates, and (3) total inner
  iterations. Defaults to `c(150, 7.5e4, 1e6)`.

- tol:

  Convergence tolerance on parameter change between iterations. Defaults
  to `1e-7`.

- n.knots:

  Integer vector of length 2 controlling internal knots for non-uniform
  bases. The first entry sets quantile knots between `range.quant`; the
  second sets equally spaced knots outside that range. Defaults to
  `c(8, 2)` for M-splines and `c(0, 20)` otherwise.

- n.events_basis:

  Integer giving the number of fully observed outcomes per uniform basis
  element. Must lie in `[1, floor(n.obs/2)]`. Defaults to
  `round(3.5 * log(n.obs) - 7.5)` when valid.

- range.quant:

  Length-2 numeric vector giving the quantile range used when setting
  quantile knots for non-uniform bases. Defaults to `c(0.075, 0.9)`.

- cover.sigma.quant:

  Proportion of fully observed outcomes targeted within the 0.025-0.975
  interval of truncated Gaussian bases tied to quantile knots. Defaults
  to `0.25`.

- cover.sigma.fixed:

  Proportion of the outcome range targeted within the 0.025-0.975
  interval of untruncated Gaussian bases tied to fixed knots. Defaults
  to `0.25`.

- min.theta:

  Minimum baseline hazard parameter value reported; estimates below are
  treated as zero (active constraints). Defaults to `1e-10`.

- penalty:

  Integer specifying penalty order. First- and second-order penalties
  are available for `"uniform"` and `"gaussian"` bases; `"epanechikov"`
  uses second-order; `"msplines"` uses `order - 1`. Defaults to `2`.

- order:

  Integer order for `"msplines"` and `"epanechikov"` bases (default
  `3`). Order 1 M-splines correspond to uniform bases; order 2 to
  triangular bases.

- kappa:

  Step-size reduction factor (\>1) used when the penalised likelihood
  fails to increase. Defaults to `1 / 0.6`.

- epsilon:

  Length-2 numeric vector giving safeguards for survival and baseline
  hazard values to avoid logarithm issues. Defaults to
  `c(1e-16, 1e-10)`.

- ties:

  Strategy for handling duplicated fully observed outcomes when defining
  knot sequences. Use `"epsilon"` to jitter duplicates with small random
  noise; use `"unique"` to drop duplicates. Defaults to `"epsilon"`.

- seed:

  Optional seed (integer vector compatible with `.Random.seed`) used
  when `ties == "epsilon"`; preserves the current RNG state when set.

## Value

A list with validated control settings (class `"coxph_mpl.control"`).

## See also

\[coxph_mpl()\]
