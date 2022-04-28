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

Download the file `HangSeng.csv` from [https://www.econ.cam.ac.uk/DCS/data.html](https://www.econ.cam.ac.uk/DCS/data.html) and import it into R:

```r
# Using the library readr
hang_seng <- readr::read_csv(
  "HangSeng.csv", col_types = cols(
    ...1 = col_skip(),
    date = col_date(format = "%d/%m/%Y")))
hang_seng <- zoo::zoo(hang_seng$returns, order.by = hang_seng$date)
```

Import the library, create a `BetaTGARCHfit` object and use the `fit` method to estimate the model.

```r
library(betatgarch)

fit_obj <- BetaTGARCHfit$new(hang_seng, n = 5)

# Estimate the model without restricting the parameter space
# to a region that satisfies an empirical version of the 
# Lyapunov condition.
fit_obj$fit(restrict = FALSE)

# Estimate the model over a restricted parameter space that
# satisfies an empirical version of the Lyapunov condition.
fit_obj$fit(restrict = TRUE)
```

Consult the help for `BetaTGARCHfit` for additional information.

### Simulation Example

One can simulate from either an estimated model or by supplying the initial value for the conditional variance and the model coefficients/parameters.

```r
# Simulate 252 observations using an existing model.
sim_obj01 <- BetaTGARCHsim$new(model = fit_obj$fit)
sim_obj01$simulate(252)

# Simulate 252 observations by passing in values for
# the conditional variance and the model
# coefficients/parameters.
sim_obj02 <- BetaTGARCHsim$new(
  f_0 = fit_obj$f_0(), coef = fit_obj$coef())
sim_obj02$simulate(252)
```

You can also simulate values by supplying your own vector of innovations. Consult the help for `BetaTGARCHsim` for additional information.

## A Note on Estimation

The package estimates parameters by minimising the negative of the log-likelihood. Additionally, it treats the first `n` observations as a pre-sample. By default, `n = 5`. This pre-sample is used to set the initial conditional variance used in the model recursion. Pre-sample observations do not enter the log-likelihood calculation.
