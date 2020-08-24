#' Utility to get jackknife pseudo observations of cumulative incidence
#'
#' @param marginal.estimate A survfit object with no covariates
#' @param time Time at which to calculate the obs
#' @param cause which cause
#' @param mr Model response of the survival object
#' @return A vector of pseudo-observations

get_jackknife <- function(marginal.estimate, time, cause, mr) {

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



