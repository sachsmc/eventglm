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
#' @param formula.censoring A one sided formula (e.g., \code{~ x1 + x2})
#'   specifying the model for the censoring distribution. If NULL, uses the same
#'   mean model as for the outcome.
#' @param ipcw.method Which method to use for calculation of inverse
#'   probability of censoring weighted pseudo observations. "binder" the
#'   default, uses the number of observations as the denominator, while the
#'   "hajek" method uses the sum of the weights as the denominator.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param weights an optional vector of ‘prior weights’ to be used in the
#'   fitting process. Should be NULL or a numeric vector.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#'   contain \code{NA}s. The default is set by the \code{na.action} setting of
#'   \link[base]{options}, and is \link[stats]{na.fail} if that is unset. The
#'   ‘factory-fresh’ default is \link[stats]{na.omit}. Another possible value is
#'   NULL, no action. Value \link[stats]{na.exclude} can be useful.
#' @param offset this can be used to specify an a priori known component to be
#'   included in the linear predictor during fitting. This should be NULL or a
#'   numeric vector of length equal to the number of cases. One or more
#'   \link[stats]{offset} terms can be included in the formula instead or as
#'   well, and if more than one is specified their sum is used. See
#'   \link[stats]{model.offset}.
#' @param control a list of parameters for controlling the fitting process. This
#'   is passed to \link[stats]{glm.control}.
#' @param model a logical value indicating whether model frame should be
#'   included as a component of the returned value.
#' @param x logical value indicating whether the model matrix used in the
#'   fitting process should be returned as components of the returned value.
#' @param y logical value indicating whether the response vector
#'   (pseudo-observations) used in the fitting process should be returned as
#'   components of the returned value.
#' @param singular.ok logical; if FALSE a singular fit is an error.
#' @param contrasts an optional list. See the contrasts.arg of
#'   \link[stats]{model.matrix.default}.
#' @param ... Other arguments passed to \link[stats]{glm.fit}
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
                      model.censoring = "independent", formula.censoring = NULL,
                      ipcw.method = "binder",
                      data,
                      weights, subset,
                      na.action, offset,
                      control = list(...), model = FALSE,
                      x = TRUE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...) {


    stopifnot(length(time) == 1)
    cal <- match.call()

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

       POi <- get_pseudo_cuminc(marginal.estimate, time, cause, mr)


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

        POi <- stratified.jacks[order(new.order)]



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
        if(ipcw.method == "binder") {

            POi <- theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i) mean(XXi[-i])))


        } else if(ipcw.method == "hajek") {

            POi <- nn * theta.n - (nn - 1) *
                (sapply(1:length(XXi), function(i) sum(XXi[-i]) / sum(Ii[-i] / Gi[-i])))


        } else {
            stop("Weighting method ", ipcw.method, " not available, options are 'binder' or 'hajek'")
        }

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
        if(ipcw.method == "binder") {

            POi <- theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i) mean(XXi[-i])))


        } else if(ipcw.method == "hajek") {

            POi <- nn * theta.n - (nn - 1) *
                (sapply(1:length(XXi), function(i) sum(XXi[-i]) / sum(Ii[-i] / Gi[-i])))


        } else {
            stop("Weighting method ", ipcw.method, " not available, options are 'binder' or 'hajek'")
        }


    } else {

        stop("Model censoring type '", model.censoring, "' not supported. Options are 'independent', 'stratified', 'aareg' or 'coxph'.")

    }

    nn <- length(POi)

    newdata[["pseudo.vals"]] <- c(POi)
    newdata[["pseudo.time"]] <- rep(time, each = nn)


    ## get stuff ready for glm.fit
    formula2 <- update.formula(formula, pseudo.vals ~ .)

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset",
                 "weights", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf[["formula"]] <- formula2
    mf[["data"]] <- newdata
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    control <- do.call("glm.control", control)

    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!stats::is.empty.model(mt)) {
        model.matrix(mt, mf, contrasts)
     } else {
         matrix(, NROW(Y), 0L)
     }
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    if (is.null(weights)) {
        weights <- rep.int(1, nrow(mf))
    }
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                          length(offset), NROW(Y)), domain = NA)
    }
    mustart <- rep(mean(newdata$pseudo.vals), nrow(mf))


    fit <- stats::glm.fit(X, Y, weights = weights,
                          family = quasi(link = link, variance = "constant"),
                          mustart = mustart,
                          intercept = attr(mt, "intercept") > 0L, singular.ok = singular.ok)

    if (model) {
        fit$model <- mf
    }
    fit$na.action <- attr(mf, "na.action")
    if (x) {
        fit$x <- X
    }
    if (!y) {
        fit$y <- NULL
    }
    fit.lin <- structure(c(fit, list(call = cal, formula = formula, terms = mt,
                      data = data, offset = offset, control = list(), method = "glm.fit",
                      contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,
                                                                              mf))),
                      class = c(fit$class, c("glm", "lm")))

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
#' @param ipcw.method Which method to use for calculation of inverse
#'   probability of censoring weighted pseudo observations. "binder" the
#'   default, uses the number of observations as the denominator, while the
#'   "hajek" method uses the sum of the weights as the denominator.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param weights an optional vector of ‘prior weights’ to be used in the
#'   fitting process. Should be NULL or a numeric vector.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#'   contain \code{NA}s. The default is set by the \code{na.action} setting of
#'   \link[base]{options}, and is \link[stats]{na.fail} if that is unset. The
#'   ‘factory-fresh’ default is \link[stats]{na.omit}. Another possible value is
#'   NULL, no action. Value \link[stats]{na.exclude} can be useful.
#' @param offset this can be used to specify an a priori known component to be
#'   included in the linear predictor during fitting. This should be NULL or a
#'   numeric vector of length equal to the number of cases. One or more
#'   \link[stats]{offset} terms can be included in the formula instead or as
#'   well, and if more than one is specified their sum is used. See
#'   \link[stats]{model.offset}.
#' @param control a list of parameters for controlling the fitting process. This
#'   is passed to \link[stats]{glm.control}.
#' @param model a logical value indicating whether model frame should be
#'   included as a component of the returned value.
#' @param x logical value indicating whether the model matrix used in the
#'   fitting process should be returned as components of the returned value.
#' @param y logical value indicating whether the response vector
#'   (pseudo-observations) used in the fitting process should be returned as
#'   components of the returned value.
#' @param singular.ok logical; if FALSE a singular fit is an error.
#' @param contrasts an optional list. See the contrasts.arg of
#'   \link[stats]{model.matrix.default}.
#' @param ... Other arguments passed to \link[stats]{glm.fit}
#'
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
                     model.censoring = "independent", formula.censoring = NULL,
                     ipcw.method = "binder",
                     data,
                     weights, subset,
                     na.action, offset,
                     control = list(...), model = FALSE,
                     x = TRUE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...) {

    stopifnot(length(time) == 1)
    cal <- match.call()

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

        if(ipcw.method == "binder") {

            POi <- theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i) mean(XXi[-i])))


        } else if(ipcw.method == "hajek") {

            POi <- nn * theta.n - (nn - 1) *
                (sapply(1:length(XXi), function(i) sum(XXi[-i]) / sum(Ii[-i] / Gi[-i])))


        } else {
            stop("Weighting method ", ipcw.method, " not available, options are 'binder' or 'hajek'")
        }


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

        if(ipcw.method == "binder") {

            POi <- theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i) mean(XXi[-i])))


        } else if(ipcw.method == "hajek") {

            POi <- nn * theta.n - (nn - 1) *
                (sapply(1:length(XXi), function(i) sum(XXi[-i]) / sum(Ii[-i] / Gi[-i])))


        } else {
            stop("Weighting method ", ipcw.method, " not available, options are 'binder' or 'hajek'")
        }

    } else {

        stop("Model censoring type '", model.censoring, "' not supported. Options are 'independent', 'stratified', 'aareg' or 'coxph'.")

    }



    newdata[["pseudo.vals"]] <- c(POi)
    newdata[["pseudo.time"]] <- rep(time, each = length(POi))

    ## get stuff ready for glm.fit
    formula2 <- update.formula(formula, pseudo.vals ~ .)

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset",
                 "weights", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf[["formula"]] <- formula2
    mf[["data"]] <- newdata
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    control <- do.call("glm.control", control)

    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!stats::is.empty.model(mt)) {
        model.matrix(mt, mf, contrasts)
    } else {
        matrix(, NROW(Y), 0L)
    }
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    if (is.null(weights)) {
        weights <- rep.int(1, nrow(mf))
    }
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                          length(offset), NROW(Y)), domain = NA)
    }
    mustart <- rep(mean(newdata$pseudo.vals), nrow(mf))


    fit <- stats::glm.fit(X, Y, weights = weights,
                          family = quasi(link = link, variance = "constant"),
                          mustart = mustart,
                          intercept = attr(mt, "intercept") > 0L, singular.ok = singular.ok)

    if (model) {
        fit$model <- mf
    }
    fit$na.action <- attr(mf, "na.action")
    if (x) {
        fit$x <- X
    }
    if (!y) {
        fit$y <- NULL
    }
    fit.lin <- structure(c(fit, list(call = cal, formula = formula, terms = mt,
                                     data = data, offset = offset, control = list(), method = "glm.fit",
                                     contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,
                                                                                             mf))),
                         class = c(fit$class, c("glm", "lm")))

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

