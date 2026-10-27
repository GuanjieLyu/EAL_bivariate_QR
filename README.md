# CJQR: Copula-Based Joint Quantile Regression

**CJQR** is an R package for analysing paired continuous outcomes using
enriched asymmetric-Laplace (EAL) marginal distributions and a bivariate
copula. Additive predictors accommodate parametric effects and nonlinear
smooth effects, with different covariates and spline bases allowed in the
two margins.

The package provides conditional quantile predictions, smooth-effect plots,
auxiliary-bootstrap pointwise confidence intervals, numerical convergence
checks and Rosenblatt residual diagnostics.

## Downloads

- [R source package: CJQR 0.22.2](CJQR_0.22.2.tar.gz)
- [Reference manual](CJQR_0.22.2_reference_manual.pdf)

## Installation

Download the source package above. In R, install the dependencies and then
the downloaded archive:

```r
install.packages(c("mgcv", "qgam", "survival"))
install.packages("CJQR_0.22.2.tar.gz", repos = NULL, type = "source")

library(CJQR)
```

Use the full path to the archive if it is not in your working directory.
CJQR requires R >= 4.1.0 and compatible versions of its dependencies.
The package itself contains no compiled code and does not require LaTeX;
dependencies installed from source may require additional build tools.

## Quick start

Generate a reproducible example and fit a Gaussian-copula model with one
parametric effect and one smooth effect in each margin:

```r
library(CJQR)

dat <- cjqr_example_data(n = 200, seed = 1)

fit <- cjqr_eal(
  formulas = list(
    y1 ~ z + s(x, k = 6),
    y2 ~ z + s(x, k = 6)
  ),
  data = dat,
  copula = "gaussian",
  anchors = 0.5,
  seed = 1,
  progress = TRUE
)

summary(fit)
cjqr_fit_diagnostics(fit)
```

This example fixes the fitting level at 0.50. Joint estimation uses multiple
starting values and convergence checks, so fitting time depends on the data
and model complexity.

### Smooth effects

```r
plot(fit, eq = 1, select = 1)
plot(fit, eq = 2, select = 1)

# Plot an estimate without confidence limits.
plot(fit, eq = 1, select = 1, se = FALSE)
```

Here `eq` selects the margin and `select` selects the smooth within that
margin. For a one-dimensional smooth, the default `se = TRUE` draws the
estimate as a solid line and the 95% pointwise confidence limits as two
dashed lines, without shading. The default axes use the covariate name and
smooth label, for example `x` and `s(x)`.

### Conditional quantiles

```r
grid <- data.frame(x = seq(-1, 1, length.out = 101), z = 0)

q <- predict(fit, newdata = grid, probs = c(0.1, 0.5, 0.9))
q[, 1, ]  # First margin: rows by quantile level.
q[, 2, ]  # Second margin: rows by quantile level.

# Request different quantile levels for the two margins.
q_different <- predict(
  fit, newdata = grid,
  probs = list(y1 = 0.1, y2 = 0.9)
)
```

The `anchors` argument determines the level used to parameterize the fitted
EAL distribution. A single value fixes that level; the default
`c(0.25, 0.5, 0.75)` compares candidate fits using the penalized likelihood.
The `probs` argument selects quantiles from the fitted distribution and does
not refit the model. Within one fitted margin, quantiles share the same
non-intercept effects and differ by constant shifts.

### Confidence intervals

```r
ci <- predict(
  fit, newdata = grid, probs = c(0.1, 0.5, 0.9),
  interval = "confidence", level = 0.95,
  B = 500, seed = 2
)

head(as.data.frame(ci))
```

Quantile intervals and smooth-plot intervals use a paired-subject bootstrap
of auxiliary least-squares regressions. Auxiliary intervals are enlarged,
if necessary, to contain the unchanged joint estimates. The joint model is
not refitted. These are pointwise confidence intervals, not simultaneous
bands, prediction intervals or joint-model standard errors.

The default `summary(fit)` does not perform bootstrap inference. Experimental
auxiliary-bootstrap effect inference can be requested explicitly:

```r
summary(fit, inference = "auxiliary", B = 500, seed = 2)
```

The older option `inference = "reference"` remains supported. Auxiliary
inference is separate from the fitting level selected through `anchors`.
See the reference manual for the assumptions and limitations of coefficient,
smooth-effect and dependence inference.

### Residual diagnostics

```r
diagnostics <- cjqr_residual_diagnostics(fit)
plot(diagnostics)
```

The display combines two Rosenblatt residual histograms and empirical
normal-score density contours. The first component transforms the first
response marginally; the second transforms the second response conditionally
on the first. These are graphical diagnostics, not formal goodness-of-fit
tests.

## Supported copulas

| `copula` | Dependence | Additional argument |
|---|---|---|
| `"gaussian"` | Positive or negative | None |
| `"student_t"` | Positive or negative | Fixed `df > 0`, e.g. `df = 5` |
| `"clayton"` | Positive | None |
| `"frank"` | Positive or negative | None |
| `"gumbel"` | Positive | None |

All families use EAL margins. The Student-t `df` argument controls the
copula, not the marginal distributions. For a Student-t copula, set
`copula = "student_t", df = 5` in the fitting call.

## Model scope

- Observational units are independent pairs of continuous responses.
- Marginal scales, enrichment parameters and copula dependence are constant
  across observations; the EAL polynomial degree is fixed at three.
- Both formulas require an intercept. Missing values are rejected by default;
  `na.action = stats::na.omit` uses a common complete-case sample.
- Observation weights and offsets are not supported. Coordinate-based
  smooths describe covariate effects, not residual spatial dependence.
- Auxiliary-bootstrap inference requires an adequate, full-rank fixed basis,
  finite-variance common errors independent of covariates, and positive
  continuous error density at the target quantiles. Nominal finite-sample
  coverage is not guaranteed; these procedures do not provide general
  inference for clustered, heteroscedastic or complex-survey data.

For compatible joint fits, `AIC(fit)`, `BIC(fit)` and `cjqr_eal_ic(fit)`
provide qualified effective-dimension comparison criteria. These are
conditional approximations, not ordinary unpenalized-model criteria, and
are not comparable with criteria from separate `qgam` fits. Consult
`?cjqr_eal_ic` before using them.

## Documentation and citation

```r
help("CJQR-package", package = "CJQR")
?cjqr_eal
?cjqr_eal_methods
?cjqr_residual_diagnostics

citation("CJQR")
```

The [reference manual](CJQR_0.22.2_reference_manual.pdf) documents the
functions, arguments and returned objects. Please cite the package when
using it in published work and record the package version and
`sessionInfo()` for reproducibility.

## Contact and license

**Maintainer:** Guanjie Lyu — <guanjielyu@gmail.com>

For bug reports, include a minimal reproducible example, the error message
and the output of `sessionInfo()`. Do not include confidential participant
data.

CJQR is distributed under the MIT license. The license files are included
in the source archive.
