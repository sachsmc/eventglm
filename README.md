[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![CRAN Status](https://www.r-pkg.org/badges/version/eventglm)](https://cran.r-project.org/package=eventglm)
[![Travis build status](https://travis-ci.org/sachsmc/eventglm.svg?branch=master)](https://travis-ci.org/sachsmc/eventglm)



# eventglm: Regression Models for Event History Outcomes

A user friendly, easy to understand way of doing event history regression for marginal estimands of interest, including the cumulative incidence and the restricted mean survival, using the pseudo observation framework for estimation.
The interface uses the well known formulation of a generalized linear model and allows for features including plotting of residuals, the use of sampling weights, and corrected variance estimation.

## Development status

This package is in stable development. The interface is unlikely to have major changes at this time. New features may be added or changed over time.

## Installation

```{r}
remotes::install_github("sachsmc/eventglm")
```

## Usage

The main functions users will use are `cumincglm` and `rmeanglm`. These are generalized linear regression models for the cumulative incidence and restricted mean of a censored time to event outcome, with or without competing risks. The models are specified just like `glm`, but the outcome must be a call to `Surv` (like in `coxph`), and you must specify the `time` argument (the fixed time at which the cumulative incidence or restricted mean is computed).

```{r}
library(survival)
library(eventglm)

colon.cifit <- cumincglm(Surv(time, status) ~ rx, time = 2500, data = colon)
summary(colon.cifit)
se.ci <- sqrt(diag(vcov(colon.cifit, method = "robust")))
b.ci <- colon.cifit$coefficients
```

Check out the vignette for more examples.

## References

Per Kragh Andersen and Maja Pohar Perme. Pseudo-observations in survival analysis. Statistical Methods in Medical Research, 19(1):71–99, February 2010. doi: 10.1177/0962280209105020. http://journals.sagepub.com/doi/10.1177/0962280209105020
