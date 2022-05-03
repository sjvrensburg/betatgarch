# betatgarch

An R package to estimate and simulate from the Beta-t-GARCH(1,1) model.

⚠ **THIS PACKAGE IS STILL ALPHA AND MISSING MANY FEATURES.** ⚠️

## Installation

Ensure that you are running a recent version of R and have the necessary dependencies installed. Install the package directly from this repository:

```r
# install.packages("remotes")
remotes::install_github("sjvrensburg/betatgarch")
```

## Basic Usage

At this time, the package contains functions to estimate and simulate from the Beta-t-GARCH(1,1) model with leverage.

### Estimation Example

The package includes a dataset `nasdaq` that contains monthly returns on the NASDAQ composite index.

Import the library, create a `BetaTGARCHfit` object and use the `fit` method to estimate the model.

```r
library(betatgarch)

fit_obj <- BetaTGARCHfit$new(nasdaq, n = 12)

# Estimate the model over a restricted parameter space that
# satisfies an empirical version of the Lyapunov condition.
fit_obj$fit(restrict = TRUE)
```

Consult the help for `BetaTGARCHfit` for additional information.

### Simulation Example

One can simulate from either an estimated model or by supplying the initial value for the conditional variance and the model coefficients/parameters.

```r
# Simulate 120 observations using an existing model.
sim_obj01 <- BetaTGARCHsim$new(model = fit_obj$fit)
sim_obj01$simulate(120)

# Simulate 120 observations by passing in values for
# the conditional variance and the model
# coefficients/parameters.
sim_obj02 <- BetaTGARCHsim$new(
  f_0 = fit_obj$f_0(), coef = fit_obj$coef())
sim_obj02$simulate(120)
```

Note that, occasionally, `simulate` may produce series with `NaN`s.

You can also simulate values by supplying your own vector of innovations. Consult the help for `BetaTGARCHsim` for additional information.

## A Note on Estimation

The package estimates parameters by minimising the negative of the log-likelihood. Additionally, it treats the first `n` observations as a pre-sample. By default, `n = 5`. This pre-sample is used to set the initial conditional variance used in the model recursion. Pre-sample observations do not enter the log-likelihood calculation.
