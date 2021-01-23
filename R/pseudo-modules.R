#' Compute pseudo observations under independent censoring
#'
#' Assuming completely independent censoring, i.e., censoring does not
#' depend on the survival time nor any covariates in the model, the pseudo
#' observations are calculated with the standard jackknife approach
#'
#' @param formula A formula specifying the model. The left hand side must be a
#'   \link[survival]{Surv} object specifying a right censored survival or
#'   competing risks outcome. The status indicator, normally 0=alive, 1=dead.
#'   Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). For competing
#'   risks, the event variable will be a factor, whose first level is treated as
#'   censoring. The right hand side is the usual linear combination of
#'   covariates.
#' @param time Numeric constant specifying the time at which the cumulative
#'   incidence or survival probability effect estimates are desired.
#' @param cause Numeric or character constant specifying the cause indicator of
#'   interest.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param type One of "survival", "cuminc", or "rmean"
#' @param formula.censoring Not used with this method, see \link{aareg} or \link{coxph}
#' @param ipcw.method Not used with this method
#'
#' @return A vector of jackknife pseudo observations
#'
#' @examples
#' POi <- pseudo_independent(Surv(time, status) ~ 1, 1500, cause = 1, data = colon, type = "survival")
#' mean(POi)
#'

pseudo_independent <- function(formula, time, cause = 1, data,
                        type = c("cuminc", "survival", "rmean"),
                        formula.censoring = NULL, ipcw.method = NULL){

  margformula <- update.formula(formula, . ~ 1)
  mr <- model.response(model.frame(margformula, data = data))
  marginal.estimate <- survival::survfit(margformula, data = data)

  if(type == "cuminc") {

    POi <- get_pseudo_cuminc(marginal.estimate, time, cause, mr)

  } else if(type == "survival"){

    if(marginal.estimate$type != "right") {
      stop("Survival estimand not available for outcome with censoring type", marginal.estimate$type)
    }

    POi <- 1 - get_pseudo_cuminc(marginal.estimate, time, cause, mr)

  } else if(type == "rmean") {

    POi <- get_pseudo_rmean(marginal.estimate, time, cause, mr)

  }

  POi

}


#' Compute pseudo observations using stratified jackknife
#'
#' Assuming that the censoring depends on covariates with a finite set of levels,
#' the pseudo observations are calculated with the jackknife approach stratified
#' on those covariates.
#'
#' @param formula A formula specifying the model. The left hand side must be a
#'   \link[survival]{Surv} object specifying a right censored survival or
#'   competing risks outcome. The status indicator, normally 0=alive, 1=dead.
#'   Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). For competing
#'   risks, the event variable will be a factor, whose first level is treated as
#'   censoring. The right hand side is the usual linear combination of
#'   covariates.
#' @param time Numeric constant specifying the time at which the cumulative
#'   incidence or survival probability effect estimates are desired.
#' @param cause Numeric or character constant specifying the cause indicator of
#'   interest.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param type One of "survival", "cuminc", or "rmean"
#' @param formula.censoring A right-sided formula specifying which variables to
#'   stratify on.
#' @param ipcw.method Not used with this method
#'
#' @return A vector of jackknife pseudo observations
#'
#' @examples
#' POi <- pseudo_stratified(Surv(time, status) ~ 1, 1500, cause = 1,
#'   data = colon, formula.censoring = ~ sex, type = "rmean")
#' mean(POi)
#'

pseudo_stratified <- function(formula, time, cause = 1, data,
                               type = c("cuminc", "survival", "rmean"),
                               formula.censoring = NULL, ipcw.method = NULL){

  margformula <- update.formula(formula, . ~ 1)
  mr <- model.response(model.frame(margformula, data = data))

  if(type == "cuminc") {

    thisfun <- get_pseudo_cuminc

  } else if(type == "survival"){

    if(attr(mr, "type") != "right") {
      stop("Survival estimand not available for outcome with censoring type", attr(mr, "type"))
    }

    thisfun <- function(me, tt, ca, mmr) {
      1 - get_pseudo_cuminc(me, tt, ca, mmr)
    }

  } else if(type == "rmean") {

    thisfun <- get_pseudo_rmean

  }


  strata <- interaction(model.frame(formula.censoring, data = data))
  orig.order <- 1:length(strata)
  new.order <- integer(length(strata))
  stratified.jacks <- numeric(length(strata))
  chunk <- 1
  for(i in levels(strata)) {

    thisset <- orig.order[strata == i]
    new.order[chunk:(chunk + length(thisset) - 1)] <- thisset
    mest.i <- survival::survfit(margformula, data = data[thisset, ])
    jres.i <- thisfun(mest.i, time, cause, mr[thisset, ])

    stratified.jacks[chunk:(chunk + length(thisset) - 1)] <- jres.i
    chunk <- chunk + length(thisset)
  }

  POi <- stratified.jacks[order(new.order)]

  POi

}


#' Compute censoring weighted pseudo observations
#'
#' Assuming that the censoring depends on covariates,
#' the pseudo observations are calculated with the inverse probability of
#' censoring weighted approach, where the censoring probabilities are estimated
#' using Aalen's additive hazards model.
#'
#' @param formula A formula specifying the outcome model. The left hand side must be a
#'   \link[survival]{Surv} object specifying a right censored survival or
#'   competing risks outcome. The status indicator, normally 0=alive, 1=dead.
#'   Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). For competing
#'   risks, the event variable will be a factor, whose first level is treated as
#'   censoring. The right hand side is the usual linear combination of
#'   covariates.
#' @param time Numeric constant specifying the time at which the cumulative
#'   incidence or survival probability effect estimates are desired.
#' @param cause Numeric or character constant specifying the cause indicator of
#'   interest.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param type One of "survival", "cuminc", or "rmean"
#' @param formula.censoring A right-sided formula specifying which variables to
#'   use in the model for the censoring distribution.
#' @param ipcw.method Which method to use for calculation of inverse
#'   probability of censoring weighted pseudo observations. "binder" the
#'   default, uses the number of observations as the denominator, while the
#'   "hajek" method uses the sum of the weights as the denominator.
#'
#' @return A vector of pseudo observations
#'
#' @examples
#' POi <- pseudo_aareg(Surv(time, status) ~ 1, 1500, cause = 1,
#'   data = colon, type = "rmean", formula.censoring = ~ sex + age,
#'   ipcw.method = "binder")
#'
#' mean(POi)
#'

pseudo_aareg <- function(formula, time, cause = 1, data,
                              type = c("cuminc", "survival", "rmean"),
                              formula.censoring = NULL, ipcw.method = NULL){

  margformula <- update.formula(formula, . ~ 1)
  mr <- model.response(model.frame(margformula, data = data))

  matcau <- match_cause(mr, cause)
  causen <- matcau$causen
  causec <- matcau$causec

  data$.Ci <- as.numeric(mr[, "status"] == 0)
  data$.Tci <- mr[, "time"]

  if(is.null(formula.censoring)) {
    cens.formula <- update.formula(formula, survival::Surv(.Tci, .Ci) ~ .)
    formula.censoring <- formula[-2]
  } else {
    cens.formula <- update.formula(formula.censoring, survival::Surv(.Tci, .Ci) ~ .)
  }

  predmat <- model.matrix(cens.formula, data = data)

  fitcens <- survival::aareg(cens.formula, data = data)

  tdex <- sapply(pmin(data$.Tci, time), function(t) max(c(1, which(fitcens$times <= t))))
  Gi <- numeric(length(tdex))
  for(i in 1:length(tdex)) {

    Gi[i] <- prod(1 - c(fitcens$coefficient[1:tdex[i], ] %*% t(predmat[i, , drop = FALSE])))

  }

  POi <- calc_ipcw_pos(mr, time, causen, type, ipcw.method, Gi)

  POi


}



#' Compute censoring weighted pseudo observations
#'
#' Assuming that the censoring depends on covariates,
#' the pseudo observations are calculated with the inverse probability of
#' censoring weighted approach, where the censoring probabilities are estimated
#' using Cox's proportional hazards model.
#'
#' @param formula A formula specifying the outcome model. The left hand side must be a
#'   \link[survival]{Surv} object specifying a right censored survival or
#'   competing risks outcome. The status indicator, normally 0=alive, 1=dead.
#'   Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). For competing
#'   risks, the event variable will be a factor, whose first level is treated as
#'   censoring. The right hand side is the usual linear combination of
#'   covariates.
#' @param time Numeric constant specifying the time at which the cumulative
#'   incidence or survival probability effect estimates are desired.
#' @param cause Numeric or character constant specifying the cause indicator of
#'   interest.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param type One of "survival", "cuminc", or "rmean"
#' @param formula.censoring A right-sided formula specifying which variables to
#'   use in the model for the censoring distribution.
#' @param ipcw.method Which method to use for calculation of inverse
#'   probability of censoring weighted pseudo observations. "binder" the
#'   default, uses the number of observations as the denominator, while the
#'   "hajek" method uses the sum of the weights as the denominator.
#'
#' @return A vector of pseudo observations
#'
#' @examples
#' POi <- pseudo_coxph(Surv(time, status) ~ 1, 1500, cause = 1,
#'   data = colon, type = "survival", formula.censoring = ~ sex + age,
#'   ipcw.method = "hajek")
#'
#' mean(POi)
#'

pseudo_coxph <- function(formula, time, cause = 1, data,
                         type = c("cuminc", "survival", "rmean"),
                         formula.censoring = NULL, ipcw.method = NULL){

  margformula <- update.formula(formula, . ~ 1)
  mr <- model.response(model.frame(margformula, data = data))

  matcau <- match_cause(mr, cause)
  causen <- matcau$causen
  causec <- matcau$causec

  data$.Ci <- as.numeric(mr[, "status"] == 0)
  data$.Tci <- mr[, "time"]

  if(is.null(formula.censoring)) {
    cens.formula <- update.formula(formula, survival::Surv(.Tci, .Ci) ~ .)
    formula.censoring <- formula[-2]
  } else {
    cens.formula <- update.formula(formula.censoring, survival::Surv(.Tci, .Ci) ~ .)
  }

  predmat <- model.matrix(cens.formula, data = data)

  fitcens <- survival::coxph(cens.formula, data = data, x = TRUE)
  coxsurv <- survival::survfit(fitcens, newdata = data)
  tdex <- sapply(pmin(data$.Tci, time), function(t) max(c(1, which(coxsurv$time <= t))))
  Gi <- coxsurv$surv[cbind(tdex,1:ncol(coxsurv$surv))]

  POi <- calc_ipcw_pos(mr, time, causen, type, ipcw.method, Gi)

  POi


}

