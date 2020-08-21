#' Generalized linear models for cumulative incidence
#'
#' Using pseudo observations for the cumulative incidence, this function then runs a generalized
#' linear model and estimates the variance correctly according to Overgaard et al (2018). The
#' link function can be "identity" for estimating differences in the cumulative incidence, "log"
#' for estimating ratios, and any of the other link functions supported by \link[stats]{quasi}.
#'
#' @return A pseudoglm object, with its own methods for print, summary, and vcov. It inherits from glm, so predict and other glm methods are supported.
#'
#' @param formula A formula specifying the model. The left hand side must be a \link[prodlim]{Hist} object. The right hand side is the usual linear combination of covariates.
#' @param time Numeric constant specifying the time at which the cumulative incidence or survival probability effect estimates are desired.
#' @param cause Numeric or character constant specifying the cause indicator of interest.
#' @param link Link function for the cumulative incidence regression model.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param ... Other arguments passed to \link[stats]{glm} such as weights, subset, etc.
#'
#' @export
#'
cumincglm <- function(formula, time, cause = 1, link = "identity", data, ...) {


    stopifnot(length(time) == 1)

    marginal.estimate <- survival::survfit(update.formula(formula, . ~ 1), data = data)

    if(max(marginal.estimate$time) < time){
        stop("Requested time is greater than largest observed time")
    }

    newdata <- do.call(rbind, lapply(1:length(time), function(i) data))
    mr <- eval(terms(update(formula, . ~ 1))[[2]], envir = data)

    if(marginal.estimate$type == "mright") {
        ## match cause
        states <- marginal.estimate$states[-1] # remove censoring
        if(is.numeric(cause)) {
            stopifnot(cause <= length(states))
            causec <- states[cause]
            causen <- cause
        } else {
            stopifnot(length(match(cause, states)) > 0)
            causen <- match(cause, states)[1]
            causec <- cause
        }

        jackk <- jackknife.competing.risks2(marginal.estimate, times = time,
                                            cause = causec, mr)

        dmatframe <- as.matrix(model.frame(update.formula(formula, .~1), data = newdata)[, 1])
        colnames(dmatframe) <- c("time", "status")

        datamat <- cbind(dmatframe[, "time"],
                         dmatframe[, "status"] != 0,
                         dmatframe[, "status"] == causen,
                         dmatframe[, "status"] == 0)

    } else if(marginal.estimate$type == "right") {

        causen <- 1
        jackk <- 1 - jackknife.survival2(marginal.estimate, times = time,
                                     mr)

        datamat <- cbind(mr[, "time"],
                         mr[, "status"] != 0, ## not censored indicator
                         mr[, "status"] == causen,
                         mr[, "status"] == 0 ## censored indicator
        )


    } else {
        stop("Survival outcome type", marginal.estimate$type, "not supported.")
    }


    nn <- length(jackk)

    newdata[["pseudo.vals"]] <- c(jackk)
    newdata[["pseudo.time"]] <- rep(time, each = nn)

    newdata[["startmu"]] <- rep(mean(newdata$pseudo.vals), nrow(newdata))

    fit.lin <- stats::glm(update.formula(formula, pseudo.vals ~ .),
                          family = quasi(link = link, variance = "constant"),
                          mustart = startmu,
                          data = newdata, x = TRUE, ...)


    fit.lin$datamat <- datamat
    fit.lin$time <- time
    fit.lin$cause <- cause
    fit.lin$link <- link
    fit.lin$type <- "cuminc"

    class(fit.lin) <- c("pseudoglm", class(fit.lin))

    fit.lin

}



#' Generalized linear models for restricted mean
#'
#' Using pseudo observations, this function then runs a generalized
#' linear model. The link function can be "identity" for estimating
#' differences in the restricted mean, "log"
#' for estimating ratios, and any of the other link functions supported by \link[stats]{quasi}.
#'
#' @return A pseudoglm object, with its own methods for print, summary, and vcov. It inherits from glm, so predict and other glm methods are supported.
#'
#' @param formula A formula specifying the model. The left hand side must be a \link[prodlim]{Hist} object. The right hand side is the usual linear combination of covariates.
#' @param time Numeric constant specifying the time at which the cumulative incidence or survival probability effect estimates are desired.
#' @param cause Character constant specifying the cause indicator of interest.
#' @param link Link function for the cumulative incidence regression model.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param ... Other arguments passed to \link[stats]{glm} such as weights, subset, etc.
#'
#' @export
#'
rmeanglm <- function(formula, time, cause = "1", link = "identity", data, ...) {


    stopifnot(length(time) == 1)

    marginal.estimate <- prodlim::prodlim(update.formula(formula, . ~ 1), data = data)

    outcome <- model.response(model.frame(update.formula(formula, . ~ 1), data = data))
    thistype <- ifelse(attr(outcome, "model") == "survival", "surv", "cuminc")
    newdata <- do.call(rbind, lapply(1:length(time), function(i) data))

    uptimes <- c(0, sort(unique(outcome[outcome[, 1] <= time, 1])))

    sfit <- predict(marginal.estimate,
                    type = thistype,
                    cause = cause, times  = uptimes)

    jackk <- prodlim::jackknife(marginal.estimate, times = uptimes, cause = cause)

    smat <- matrix(rep(sfit, each = nrow(jackk)), nrow = nrow(jackk), ncol = ncol(jackk))
    Smi <- (jackk - nrow(outcome) * smat) / (1 - nrow(outcome))

    POi <- pseudo_rmst2(sfit, Smi, uptimes, time, type = thistype)
    # individuals in rows, times in columns

    newdata[["pseudo.vals"]] <- c(POi)
    newdata[["pseudo.time"]] <- rep(time, each = nrow(jackk))

    newdata[["startmu"]] <- rep(mean(newdata$pseudo.vals), nrow(newdata))

    fit.lin <- stats::glm(update.formula(formula, pseudo.vals ~ .),
                          family = quasi(link = link, variance = "constant"),
                          mustart = startmu,
                          data = newdata, x = TRUE, ...)

    ## update variance estimate

    if(marginal.estimate$model == "survival") {

        datamat <- cbind(marginal.estimate$model.response[, "time"],
                         marginal.estimate$model.response[, "status"] != 0, ## not censored indicator
                         as.character(marginal.estimate$model.response[, "status"]) == cause,
                         marginal.estimate$model.response[, "status"] == 0 ## censored indicator
        )

    } else {

        dmatframe <- as.matrix(model.frame(update.formula(formula, .~1), data = newdata)[, 1])
        colnames(dmatframe) <- c("time", "status", "event")

        datamat <- cbind(dmatframe[, "time"],
                         dmatframe[, "status"] != 0,
                         dmatframe[, "event"] == as.numeric(cause),
                         dmatframe[, "status"] == 0)

    }


    #    fit.lin$corrected.vcov <- ovg.vcov
    fit.lin$datamat <- datamat
    fit.lin$time <- time
    fit.lin$cause <- cause
    fit.lin$link <- link
    fit.lin$type <- "rmean"

    class(fit.lin) <- c("pseudoglm", class(fit.lin))

    fit.lin

}


#' Compute restricted mean survival
#'
#'
#' @param fit Survfit
#' @param tmax Max time
pseudo_rmst2 <- function(sfit, jacks, times, tmax, type = "cuminc") {
    # extract the RMST, and the leverage of each subject on the RMST

    rsum <- function(y, x) {  # sum of rectangles
        keep <- which(x < tmax)
        width <- diff(c(x[keep], tmax))
        sum(width * y[keep])
    }

    if (type == "survival") { # ordinary survival
        rmst <- rsum(sfit, times)
        ijack <- apply(jacks, 1, rsum)
        length(ijack) * rmst - (length(ijack) -1)* ijack
    }
    else {
        rmst <- rsum(sfit, times)
        ijack <- apply(jacks, 1, rsum, times)
        length(ijack) * rep(rmst, each= length(ijack)) - (length(ijack)-1)*ijack
    }
}


