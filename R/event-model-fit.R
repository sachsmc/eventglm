#' Generalized linear models for cumulative incidence
#'
#' Using pseudo observations for the cumulative incidence, this function then
#' runs a generalized linear model and estimates the parameters representing
#' contrasts in the cumulative incidence at a particular time (specified by the
#' \code{time} argument) across covariate values. The link function can be
#' "identity" for estimating differences in the cumulative incidence, "log" for
#' estimating ratios, and any of the other link functions supported by
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
#'   the model. Can also be a custom function, see Details and the
#'   "Extending eventglm" vignette.
#' @param formula.censoring A one sided formula (e.g., \code{~ x1 + x2})
#'   specifying the model for the censoring distribution. If NULL, uses the same
#'   mean model as for the outcome.
#' @param ipcw.method Which method to use for calculation of inverse
#'   probability of censoring weighted pseudo observations. "binder" the
#'   default, uses the number of observations as the denominator, while the
#'   "hajek" method uses the sum of the weights as the denominator.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param survival Set to TRUE to use survival (one minus the cumulative incidence)
#'   as the outcome. Not available for competing risks models.
#' @param weights an optional vector of 'prior weights' to be used in the
#'   fitting process. Should be NULL or a numeric vector.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#'   contain \code{NA}s. The default is set by the \code{na.action} setting of
#'   \link[base]{options}, and is \link[stats]{na.fail} if that is unset. The
#'   'factory-fresh' default is \link[stats]{na.omit}. Another possible value is
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
#' @details The argument "model.censoring" determines how the pseudo observations
#' are calculated. This can be the name of a function or the function itself, which
#' must have arguments "formula", "time", "cause", "data", "type",
#' "formula.censoring", and "ipcw.method". If it is the name of a function, this code
#' will look for a function with the prefix "pseudo_" first, to avoid clashes with
#' related methods such as coxph. The function then must return a vector
#' of pseudo observations, one for each subject in data which are used in subsequent
#' calculations. For examples of the implementation, see the "pseudo-modules.R"
#' file, or the vignette "Extending eventglm".
#'
#'
#' @export
#'
#' @examples
#'     cumincipcw <- cumincglm(Surv(etime, event) ~ age + sex,
#'          time = 200, cause = "pcm", link = "identity",
#'          model.censoring = "independent", data = mgus2)
#' # stratified on only the categorical covariate
#'      cumincipcw2 <- cumincglm(Surv(etime, event) ~ age + sex,
#'                          time = 200, cause = "pcm", link = "identity",
#'                          model.censoring = "stratified",
#'                          formula.censoring = ~ sex, data = mgus2)

cumincglm <- function(formula, time, cause = 1, link = "identity",
                      model.censoring = "independent", formula.censoring = NULL,
                      ipcw.method = "binder",
                      data,
                      survival = FALSE,
                      weights, subset,
                      na.action, offset,
                      control = list(...), model = FALSE,
                      x = TRUE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...) {


    stopifnot(length(time) == 1 & is.numeric(time))
    cal <- match.call()

    mr <- model.response(model.frame(update.formula(formula, . ~ 1), data = data))

    newdata <- do.call(rbind, lapply(1:length(time), function(i) data))

    matcau <- match_cause(mr, cause)
    causec <- matcau$causec
    causen <- matcau$causen

    otype <- if(survival) "survival" else "cuminc"

    if(is.character(model.censoring)) {
        pseudo_function <- tryCatch(get(paste0("pseudo_", model.censoring),
                               mode = "function", envir = parent.frame()),
                               error = function(e) NULL)

        if(is.null(pseudo_function)) { # try without the pseudo prefix
            pseudo_function <- get(model.censoring,
                                   mode = "function", envir = parent.frame())
        }
    }

    check_pseudo_function <- names(formals(pseudo_function)) ==
        c("formula", "time", "cause", "data", "type",
                                         "formula.censoring", "ipcw.method")

    if(!all(check_pseudo_function)) {

        missargs <- c("formula", "time", "cause", "data", "type",
                      "formula.censoring", "ipcw.method")[!check_pseudo_function]
        stop("Arguments", paste(missargs, collapse = ", "), "are missing from model.censoring function")

    }

    POi <- pseudo_function(formula = formula, time = time, cause = cause,
                           data = data, type = otype, formula.censoring = formula.censoring,
                           ipcw.method = ipcw.method)

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
    fit.lin$type <- if(survival) "survival" else "cuminc"

    class(fit.lin) <- c("pseudoglm", class(fit.lin))

    fit.lin

}


#' Generalized linear models for the restricted mean survival
#'
#' Using pseudo observations for the restricted mean, or the restricted mean
#' lifetime lost in the competing risks case, this function then runs a
#' generalized linear model to estimate associations of the restricted
#' mean/lifetime lost up to a particular time (specified by the \code{time}
#' argument) with covariates. The link function can be "identity" for estimating
#' differences in the restricted mean, "log" for estimating ratios, and any of
#' the other link functions supported by \link[stats]{quasi}.
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
#'   the model. Can also be a custom function, see Details and the
#'   "Extending eventglm" vignette.
#' @param formula.censoring A one sided formula (e.g., \code{~ x1 + x2})
#'   specifying the model for the censoring distribution. If NULL, uses the same
#'   mean model as for the outcome.
#' @param ipcw.method Which method to use for calculation of inverse
#'   probability of censoring weighted pseudo observations. "binder" the
#'   default, uses the number of observations as the denominator, while the
#'   "hajek" method uses the sum of the weights as the denominator.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param weights an optional vector of 'prior weights' to be used in the
#'   fitting process. Should be NULL or a numeric vector.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#'   contain \code{NA}s. The default is set by the \code{na.action} setting of
#'   \link[base]{options}, and is \link[stats]{na.fail} if that is unset. The
#'   'factory-fresh' default is \link[stats]{na.omit}. Another possible value is
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
#' @details The argument "model.censoring" determines how the pseudo observations
#' are calculated. This can be the name of a function or the function itself, which
#' must have arguments "formula", "time", "cause", "data", "type",
#' "formula.censoring", and "ipcw.method". If it is the name of a function, this code
#' will look for a function with the prefix "pseudo_" first, to avoid clashes with
#' related methods such as coxph. The function then must return a vector
#' of pseudo observations, one for each subject in data which are used in subsequent
#' calculations. For examples of the implementation, see the "pseudo-modules.R"
#' file, or the vignette "Extending eventglm".

#'
#' @examples
#'     cumincipcw <- rmeanglm(Surv(etime, event) ~ age + sex,
#'          time = 200, cause = "pcm", link = "identity",
#'          model.censoring = "independent", data = mgus2)
#' # stratified on only the categorical covariate
#'      cumincipcw2 <- rmeanglm(Surv(etime, event) ~ age + sex,
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

    newdata <- do.call(rbind, lapply(1:length(time), function(i) data))

    matcau <- match_cause(mr, cause)
    causec <- matcau$causec
    causen <- matcau$causen

    otype <- "rmean"

    if(is.character(model.censoring)) {
        pseudo_function <- tryCatch(get(paste0("pseudo_", model.censoring),
                                        mode = "function", envir = parent.frame()),
                                    error = function(e) NULL)

        if(is.null(pseudo_function)) { # try without the pseudo prefix
            pseudo_function <- get(model.censoring,
                                   mode = "function", envir = parent.frame())
        }
    }

    check_pseudo_function <- names(formals(pseudo_function)) == c("formula", "time", "cause", "data", "type",
                                                                  "formula.censoring", "ipcw.method")

    if(!all(check_pseudo_function)) {

        missargs <- c("formula", "time", "cause", "data", "type",
                      "formula.censoring", "ipcw.method")[!check_pseudo_function]
        stop("Arguments", paste(missargs, collapse = ", "), "are missing from model.censoring function")

    }

    POi <- pseudo_function(formula = formula, time = time, cause = cause,
                           data = data, type = otype, formula.censoring = formula.censoring,
                           ipcw.method = ipcw.method)

    nn <- length(POi)

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

