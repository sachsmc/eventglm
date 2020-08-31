#' Generalized linear models for cumulative incidence
#'
#' Using pseudo observations for the cumulative incidence, this function then
#' runs a generalized linear model and estimates the variance correctly
#' according to Overgaard et al (2018). The link function can be "identity" for
#' estimating differences in the cumulative incidence, "log" for estimating
#' ratios, and any of the other link functions supported by \link[stats]{quasi}.
#'
#' @return A pseudoglm object, with its own methods for print, summary, and
#'   vcov. It inherits from glm, so predict and other glm methods are supported.
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
#' @param link Link function for the cumulative incidence regression model.
#' @param model.censoring Type of model for the censoring distribution. Options
#'   are "stratified", which computes the pseudo-observations stratified on a
#'   set of categorical covariates, "aareg" for Aalen's additive hazards model,
#'   and "coxph" for Cox's proportional hazards model. With those options, we
#'   assume that the time to event and event indicator are conditionally
#'   independent of the censoring time, and that the censoring model is
#'   correctly specified. If "independent", we assume completely independent
#'   censoring, i.e., that the time to event and covariates are independent of
#'   the censoring time. the censoring time is independent of the covariates in
#'   the model.
#' @param formula.censoring A one sided formula (e.g., \code{~ x1 + x2}) specifying the model for the censoring
#'   distribution. If NULL, uses the same mean model as for the outcome.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param ... Other arguments passed to \link[stats]{glm} such as weights,
#'   subset, etc.
#'
#' @export
#'
#' @examples
#'     cumincipcw <- cumincglm(survival::Surv(etime, event) ~ age + sex,
#'          time = 200, cause = "pcm", link = "identity",
#'          model.censoring = "independent", data = mgus2)
#' # stratified on only the categorical covariate
#'      cumincipcw2 <- cumincglm(survival::Surv(etime, event) ~ age + sex,
#'                          time = 200, cause = "pcm", link = "identity",
#'                          model.censoring = "stratified",
#'                          formula.censoring = ~ sex, data = mgus2)

cumincglm <- function(formula, time, cause = 1, link = "identity",
                      model.censoring = "independent", formula.censoring = NULL, data, ...) {


    stopifnot(length(time) == 1)

    mr <- model.response(model.frame(update.formula(formula, . ~ 1), data = data))

    ## match cause
    if(attr(mr, "type") == "mright") {
    states <- attr(mr, "states")
    if(is.numeric(cause)) {
        stopifnot(cause <= length(states))
        causec <- states[cause]
        causen <- cause
    } else {
        stopifnot(length(match(cause, states)) > 0)
        causen <- match(cause, states)[1]
        causec <- cause
    }
    } else if(attr(mr, "type") == "right") {
        causen <- 1
    } else {
        stop("Survival outcome type ", attr(mr, "type"), " not supported.")
    }

    if(max(mr[mr[,"status"] != 0, "time"]) < time){
        stop("Requested time is greater than largest observed event time")
    }

    newdata <- do.call(rbind, lapply(1:length(time), function(i) data))
    newdata$.Ci <- as.numeric(mr[, "status"] == 0)
    newdata$.Tci <- mr[, "time"]


    if(is.null(formula.censoring)) {
        cens.formula <- update.formula(formula, survival::Surv(.Tci, .Ci) ~ .)
        formula.censoring <- formula[-2]
    } else {
        cens.formula <- update.formula(formula.censoring, survival::Surv(.Tci, .Ci) ~ .)
    }

    if(model.censoring == "independent") {

        marginal.estimate <- survival::survfit(update.formula(formula, . ~ 1), data = data)

       jackk <- get_pseudo_cuminc(marginal.estimate, time, cause, mr)


    } else if(model.censoring == "stratified") {

        strata <- interaction(model.frame(formula.censoring, data = data))
        orig.order <- 1:length(strata)
        new.order <- rep(NA, length(strata))
        stratified.jacks <- rep(NA, length(strata))
        chunk <- 1
        for(i in levels(strata)) {

            thisset <- orig.order[strata == i]
            new.order[chunk:(chunk + length(thisset) - 1)] <- thisset
            mest.i <- survival::survfit(update.formula(formula, . ~ 1), data = data[thisset, ])
            jres.i <- get_pseudo_cuminc(mest.i, time, cause, mr[thisset, ])

            stratified.jacks[chunk:(chunk + length(thisset) - 1)] <- jres.i
            chunk <- chunk + length(thisset)
        }

        jackk <- stratified.jacks[order(new.order)]



    } else if(model.censoring == "aareg") {

        predmat <- model.matrix(cens.formula, data = newdata)

        fitcens <- survival::aareg(cens.formula, data = newdata)

        tdex <- sapply(pmin(newdata$.Tci, time), function(t) max(c(1, which(fitcens$times <= t))))
        Gi <- rep(NA, length(tdex))
        for(i in 1:length(tdex)) {

            Gi[i] <- prod(1 - c(fitcens$coefficient[1:tdex[i], ] %*% t(predmat[i, , drop = FALSE])))

        }
        Vi <- as.numeric(mr[, "time"] < time & mr[, "status"] == causen)
        Ii <- as.numeric(mr[, "time"] >= time | mr[, "status"] != 0)

        nn <- length(Vi)
        theta.n <- sum(Ii * Vi / Gi) / sum(Ii / Gi)

        XXi <- Vi * Ii / Gi
        #jackk <- theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i) mean(XXi[-i])))
        jackk <- nn * theta.n - (nn - 1) *
            (sapply(1:length(XXi), function(i) sum(XXi[-i]) / sum(Ii[-i] / Gi[-i])))


    } else if(model.censoring == "coxph") {

        fitcens <- survival::coxph(cens.formula, data = newdata, x = TRUE)
        coxsurv <- survival::survfit(fitcens, newdata = newdata)
        tdex <- sapply(pmin(newdata$.Tci, time), function(t) max(c(1, which(coxsurv$time <= t))))
        Gi <- coxsurv$surv[cbind(tdex,1:ncol(coxsurv$surv))]
        Vi <- as.numeric(mr[, "time"] < time & mr[, "status"] == causen)
        Ii <- as.numeric(mr[, "time"] >= time | mr[, "status"] != 0)

        nn <- length(Vi)
        theta.n <- sum(Ii * Vi / Gi) / sum(Ii / Gi)

        XXi <- Vi * Ii / Gi
        #jackk <- theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i) mean(XXi[-i])))
        jackk <- nn * theta.n - (nn - 1) *
            (sapply(1:length(XXi), function(i) sum(XXi[-i]) / sum(Ii[-i] / Gi[-i])))


    } else {

        stop("Model censoring type '", model.censoring, "' not supported. Options are 'independent', 'stratified', 'aareg' or 'coxph'.")

    }

    nn <- length(jackk)

    newdata[["pseudo.vals"]] <- c(jackk)
    newdata[["pseudo.time"]] <- rep(time, each = nn)

    newdata[["startmu"]] <- rep(mean(newdata$pseudo.vals), nrow(newdata))

    fit.lin <- stats::glm(update.formula(formula, pseudo.vals ~ .),
                          family = quasi(link = link, variance = "constant"),
                          mustart = startmu,
                          data = newdata, x = TRUE, ...)

    datamat <- cbind(mr[, "time"],
                     mr[, "status"] != 0, ## not censored indicator
                     mr[, "status"] == causen,
                     mr[, "status"] == 0 )## censored indicator

    fit.lin$datamat <- datamat
    fit.lin$time <- time
    fit.lin$cause <- cause
    fit.lin$link <- link
    fit.lin$type <- "cuminc"

    class(fit.lin) <- c("pseudoglm", class(fit.lin))

    fit.lin

}


#' Generalized linear models for the restricted mean survival
#'
#' Using pseudo observations for the restricted mean, or the restricted mean
#' lifetime lost in the competing risks case, this function then runs a
#' generalized linear model to estimate associations with covariates. The link
#' function can be "identity" for estimating differences in the restricted mean,
#' "log" for estimating ratios, and any of the other link functions supported by
#' \link[stats]{quasi}.
#'
#' @return A pseudoglm object, with its own methods for print, summary, and
#'   vcov. It inherits from glm, so predict and other glm methods are supported.
#'
#' @param formula A formula specifying the model. The left hand side must be a
#'   \link[survival]{Surv} object specifying a right censored survival or
#'   competing risks outcome. The status indicator, normally 0=alive, 1=dead.
#'   Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). For competing
#'   risks, the event variable will be a factor, whose first level is treated as
#'   censoring. The right hand side is the usual linear combination of
#'   covariates.
#' @param time Numeric constant specifying the time up to which the restricted
#'   mean effect estimates are desired.
#' @param cause Numeric or character constant specifying the cause indicator of
#'   interest.
#' @param link Link function for the restricted mean regression model.
#' @param model.censoring Type of model for the censoring distribution. Options
#'   are "stratified", which computes the pseudo-observations stratified on a
#'   set of categorical covariates, "aareg" for Aalen's additive hazards model,
#'   and "coxph" for Cox's proportional hazards model. With those options, we
#'   assume that the time to event and event indicator are conditionally
#'   independent of the censoring time, and that the censoring model is
#'   correctly specified. If "independent", we assume completely independent
#'   censoring, i.e., that the time to event and covariates are independent of
#'   the censoring time. the censoring time is independent of the covariates in
#'   the model.
#' @param formula.censoring A one sided formula (e.g., \code{~ x1 + x2})
#'   specifying the model for the censoring distribution. If NULL, uses the same
#'   mean model as for the outcome.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param ... Other arguments passed to \link[stats]{glm} such as weights,
#'   subset, etc.
#'
#' @export
#'
#' @examples
#'     cumincipcw <- rmeanglm(survival::Surv(etime, event) ~ age + sex,
#'          time = 200, cause = "pcm", link = "identity",
#'          model.censoring = "independent", data = mgus2)
#' # stratified on only the categorical covariate
#'      cumincipcw2 <- rmeanglm(survival::Surv(etime, event) ~ age + sex,
#'                          time = 200, cause = "pcm", link = "identity",
#'                          model.censoring = "stratified",
#'                          formula.censoring = ~ sex, data = mgus2)

rmeanglm <- function(formula, time, cause = 1, link = "identity",
                     model.censoring = "independent", formula.censoring = NULL, data, ...) {

    stopifnot(length(time) == 1)

    mr <- model.response(model.frame(update.formula(formula, . ~ 1), data = data))

    ## match cause
    if(attr(mr, "type") == "mright") {
        states <- attr(mr, "states")
        if(is.numeric(cause)) {
            stopifnot(cause <= length(states))
            causec <- states[cause]
            causen <- cause
        } else {
            stopifnot(length(match(cause, states)) > 0)
            causen <- match(cause, states)[1]
            causec <- cause
        }
    } else if(attr(mr, "type") == "right") {
        causen <- 1
    } else {
        stop("Survival outcome type ", attr(mr, "type"), " not supported.")
    }

    if(max(mr[mr[,"status"] != 0, "time"]) < time){
        stop("Requested time is greater than largest observed event time")
    }

    newdata <- do.call(rbind, lapply(1:length(time), function(i) data))
    newdata$.Ci <- as.numeric(mr[, "status"] == 0)
    newdata$.Tci <- mr[, "time"]


    if(is.null(formula.censoring)) {
        cens.formula <- update.formula(formula, survival::Surv(.Tci, .Ci) ~ .)
        formula.censoring <- formula[-2]
    } else {
        cens.formula <- update.formula(formula.censoring, survival::Surv(.Tci, .Ci) ~ .)
    }

    if(model.censoring == "independent") {

        marginal.estimate <- survival::survfit(update.formula(formula, . ~ 1), data = data)
        POi <- get_pseudo_rmean(marginal.estimate, time, cause, mr)


    }  else if(model.censoring == "stratified") {

        strata <- interaction(model.frame(formula.censoring, data = data))
        orig.order <- 1:length(strata)
        new.order <- rep(NA, length(strata))
        stratified.jacks <- rep(NA, length(strata))
        chunk <- 1
        for(i in levels(strata)) {

            thisset <- orig.order[strata == i]
            new.order[chunk:(chunk + length(thisset) - 1)] <- thisset
            mest.i <- survival::survfit(update.formula(formula, . ~ 1), data = data[thisset, ])
            jres.i <- get_pseudo_rmean(mest.i, time, cause, mr[thisset, ])

            stratified.jacks[chunk:(chunk + length(thisset) - 1)] <- jres.i
            chunk <- chunk + length(thisset)
        }

        POi <- stratified.jacks[order(new.order)]

    } else if(model.censoring == "aareg") {

        predmat <- model.matrix(cens.formula, data = newdata)

        fitcens <- survival::aareg(cens.formula, data = newdata)

        tdex <- sapply(pmin(newdata$.Tci, time), function(t) max(c(1, which(fitcens$times <= t))))
        Gi <- rep(NA, length(tdex))
        for(i in 1:length(tdex)) {

            Gi[i] <- prod(1 - c(fitcens$coefficient[1:tdex[i], ] %*% t(predmat[i, , drop = FALSE])))

        }
        if(attr(mr, "type") == "mright") {
            Vi <- (time - pmin(mr[, "time"], time)) * as.numeric(mr[, "status"] == causen)

        } else {
            Vi <- pmin(mr[, "time"], time)

        }

        Ii <- as.numeric(mr[, "time"] >= time | mr[, "status"] != 0)

        nn <- length(Vi)
        theta.n <- sum(Ii * Vi / Gi) / sum(Ii / Gi)

        XXi <- Vi * Ii / Gi
        #jackk <- theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i) mean(XXi[-i])))
        POi <- nn * theta.n - (nn - 1) *
            (sapply(1:length(XXi), function(i) sum(XXi[-i]) / sum(Ii[-i] / Gi[-i])))


    } else if(model.censoring == "coxph") {

        fitcens <- survival::coxph(cens.formula, data = newdata, x = TRUE)
        coxsurv <- survival::survfit(fitcens, newdata = newdata)
        tdex <- sapply(pmin(newdata$.Tci, time), function(t) max(c(1, which(coxsurv$time <= t))))
        Gi <- coxsurv$surv[cbind(tdex,1:ncol(coxsurv$surv))]

        if(attr(mr, "type") == "mright") {
            Vi <- (time - pmin(mr[, "time"], time)) * as.numeric(mr[, "status"] == causen)

        } else {
            Vi <- pmin(mr[, "time"], time)

        }

        Ii <- as.numeric(mr[, "time"] >= time | mr[, "status"] != 0)

        nn <- length(Vi)
        theta.n <- sum(Ii * Vi / Gi) / sum(Ii / Gi)

        XXi <- Vi * Ii / Gi
        #jackk <- theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i) mean(XXi[-i])))
        POi <- nn * theta.n - (nn - 1) *
            (sapply(1:length(XXi), function(i) sum(XXi[-i]) / sum(Ii[-i] / Gi[-i])))


    } else {

        stop("Model censoring type '", model.censoring, "' not supported. Options are 'independent', 'stratified', 'aareg' or 'coxph'.")

    }



    # individuals in rows, times in columns

    newdata[["pseudo.vals"]] <- c(POi)
    newdata[["pseudo.time"]] <- rep(time, each = length(POi))

    newdata[["startmu"]] <- rep(mean(newdata$pseudo.vals), nrow(newdata))

    fit.lin <- stats::glm(update.formula(formula, pseudo.vals ~ .),
                          family = quasi(link = link, variance = "constant"),
                          mustart = startmu,
                          data = newdata, x = TRUE, ...)

    ## update variance estimate

    datamat <- cbind(mr[, "time"],
                     mr[, "status"] != 0, ## not censored indicator
                     mr[, "status"] == causen,
                     mr[, "status"] == 0 )## censored indicator


    fit.lin$datamat <- datamat
    fit.lin$time <- time
    fit.lin$cause <- cause
    fit.lin$link <- link
    fit.lin$type <- "rmean"

    class(fit.lin) <- c("pseudoglm", class(fit.lin))

    fit.lin

}

