#' Utility to get jackknife pseudo observations of cumulative incidence
#'
#' @param marginal.estimate A survfit object with no covariates
#' @param time Time at which to calculate the obs
#' @param cause which cause
#' @param mr Model response of the survival object
#' @return A vector of pseudo-observations

get_pseudo_cuminc <- function(marginal.estimate, time, cause, mr) {

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



    } else if(marginal.estimate$type == "right") {

        causen <- 1
        jackk <- 1 - jackknife.survival2(marginal.estimate, times = time,
                                         mr)




    } else {
        stop("Survival outcome type ", marginal.estimate$type, " not supported.")
    }

    jackk

}

#' Utility to get jackknife pseudo observations of restricted mean
#'
#' @param marginal.estimate A survfit object with no covariates
#' @param time Time at which to calculate the obs
#' @param cause which cause
#' @param mr Model response of the survival object
#' @return A vector of pseudo-observations

get_pseudo_rmean <- function(marginal.estimate, time, cause, mr) {

    if(attr(mr, "type") == "mright") {
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

        Smi <- leaveOneOut.competing.risks(marginal.estimate,
                                           times = time, cause = causec, mr)

    } else if(attr(mr, "type") == "right") {
        Smi <- leaveOneOut.survival(marginal.estimate, times = time, mr)
    }

    thistype <- ifelse(attr(mr, "type") == "right", "surv", "cuminc")

    uptimes <- c(0, sort(unique(mr[mr[, "time"] <= time, "time"])))

    sfit0 <- summary(marginal.estimate, times  = uptimes)

    if(thistype == "surv") {
        sfit <- sfit0$surv
    } else {
        sfit <- sfit0$pstate[, causen + 1]
    }

    pseudo_rmst2(sfit, Smi, uptimes, time, type = thistype)

}


#' Compute pseudo-observations for the restricted mean survival
#'
#'
#' @param sfit A survfit object
#' @param jacks A matrix of leave-one-out jackknife values, subjects in the rows, times in the columns
#' @param times Times at which the survival is calculated
#' @param tmax Max time
#' @param type "cuminc" or "survival"
#' @return A vector of pseudo observations for the restricted mean or lifetime lost

pseudo_rmst2 <- function(sfit, jacks, times, tmax, type = "cuminc") {
    # extract the RMST, and the leverage of each subject on the RMST

    rsum <- function(y, x) {  # sum of rectangles
        keep <- which(x < tmax)
        width <- diff(c(x[keep], tmax))
        sum(width * y[keep])
    }

    if (type == "surv") { # ordinary survival
        rmst <- rsum(sfit, times)
        ijack <- apply(cbind(1, jacks), 1, rsum, times)
        length(ijack) * rmst - (length(ijack) -1)* ijack
    }
    else {
        rmst <- rsum(sfit, times)
        ijack <- apply(cbind(0, jacks), 1, rsum, times)
        length(ijack) * rep(rmst, each= length(ijack)) - (length(ijack)-1)*ijack
    }
}


#' Compute inverse probability of censoring weights pseudo observations
#'
#' @param mr Model response object returned by \link{Surv}
#' @param time Max time
#' @param causen Cause of interest (numeric)
#' @param type Outcome type, "cuminc", "survival", or "rmean"
#' @param ipcw.method "binder" or "hajek"
#' @param Gi vector of estimated censoring probabilities

calc_ipcw_pos <- function(mr, time, causen, type, ipcw.method, Gi) {
    stopifnot(length(Gi) == nrow(mr))

    if (type == "cuminc") {
        Vi <- as.numeric(mr[, "time"] < time & mr[, "status"] == causen)

    } else if (type == "survival") {
        if (attr(mr, "type") != "right") {
            stop(
                "Survival estimand not available for outcome with censoring type",
                attr(mr, "type")
            )
        }

        Vi <-
            1 - as.numeric(mr[, "time"] < time & mr[, "status"] == causen)

    } else if (type == "rmean") {
        if (attr(mr, "type") == "mright") {
            Vi <-
                (time - pmin(mr[, "time"], time)) * as.numeric(mr[, "status"] == causen)

        } else {
            Vi <- pmin(mr[, "time"], time)

        }
    }

    Ii <- as.numeric(mr[, "time"] >= time | mr[, "status"] != 0)

    nn <- length(Vi)
    theta.n <- sum(Ii * Vi / Gi) / sum(Ii / Gi)

    XXi <- Vi * Ii / Gi
    if (ipcw.method == "binder") {
        POi <-
            theta.n + (nn - 1) * (theta.n - sapply(1:length(XXi), function(i)
                mean(XXi[-i])))


    } else if (ipcw.method == "hajek") {
        POi <- nn * theta.n - (nn - 1) *
            (sapply(1:length(XXi), function(i)
                sum(XXi[-i]) / sum(Ii[-i] / Gi[-i])))


    } else {
        stop(
            "Weighting method ",
            ipcw.method,
            " not available, options are 'binder' or 'hajek'"
        )
    }

    POi
}


#' Match cause specification against model response
#'
#' @param mr model.response as returned by \link{surv}
#' @param cause Numeric or string indicating the cause of interest
#'

match_cause <- function(mr, cause) {
    if (!is.null(attr(mr, "states"))) {
        states <- attr(mr, "states")
        if (is.numeric(cause)) {
            stopifnot(cause <= length(states))
            causec <- states[cause]
            causen <- cause
        } else {
            stopifnot(length(match(cause, states)) > 0)
            causen <- match(cause, states)[1]
            causec <- cause
        }
    } else {
        causen <- 1
        causec <- "dead"
    }
    list(causen = causen, causec = causec)

}
