#' Regression Models for Event History Outcomes
#'
#' A user friendly, easy to understand way of doing event history regression
#' for marginal estimands of interest, including the cumulative incidence and the
#' restricted mean survival, using the pseudo observation framework for estimation.
#' The interface uses the well known formulation of a generalized linear model and allows
#' for features including plotting of residuals, the use of sampling weights, and corrected
#' variance estimation.
#'
#' @docType package
#' @name eventglm
#' @useDynLib eventglm, .registration=TRUE
#' @importFrom stats coef model.frame model.matrix model.response .getXlevels model.offset model.weights naprint pnorm predict quasi residuals residuals.glm terms update update.formula vcov qnorm setNames
NULL
