#' Compute jackknife pseudo-observations of the survival function
#'
#' @param object A survfit object, with a single event (no competing risks)
#' @param times Times at which the survival is computed, must be length 1
#' @param mr Model response, the result of a call to Surv, or a matrix with two columns: "time" (observed follow up time) and "status" (0 = censored, 1 = event)
#' @return A vector of jackknifed estimates of survival at time times
#' @export
#' @examples
#'
#' sfit.surv <- survival::survfit(survival::Surv(time, status) ~ 1, data = colon)
#' mrs <- with(colon, survival::Surv(time, status))
#' pseudo.obs <- jackknife.survival2(sfit.surv, times = 1000, mrs)
#' mean(pseudo.obs)
#' # agrees with
#' summary(sfit.surv, times = 1000)


jackknife.survival2 <- function(object,times,mr){
    stopifnot(length(times) == 1)

    event.time.order <- order(mr[,"time"],-1.0 * (mr[,"status"] != 0))
    smary <- summary(object,times=times)
    S <- smary$surv

    Sk <- leaveOneOut.survival2(object,times,mr)
    N <- length(Sk)
    Jk <- N*S - (N-1) * Sk

    ## re-order the pseudo-values
    Jk <- Jk[order(event.time.order)]
    Jk
}

#' Compute jackknife pseudo-observations of the cause-specific cumulative incidence for competing risks
#'
#' @param object A survfit object, with competing events
#' @param times Times at which the cumulative incidence is computed, must be length 1
#' @param cause Value indicating for which cause the cumulative incidence is to be computed, it must match one of the values available in object (see example)
#' @param mr Model response, the result of a call to Surv, or a matrix with two columns: "time" (observed follow up time) and "status" (0 = censored, 1, ..., k = event types)
#' @return A vector of jackknifed pseudo-observations of the cause-specific cumulative incidence at time times
#' @export
#' @examples
#'
#' sfit.cuminc <- survival::survfit(survival::Surv(etime, event) ~ 1, data = mgus2)
#' mrs <- with(mgus2, survival::Surv(etime, event))
#' pseudo.obs <- jackknife.competing.risks2(sfit.cuminc, times = 200, cause = "pcm", mrs)
#' mean(pseudo.obs)
#' # agrees with
#' summary(sfit.cuminc, times = 200)


jackknife.competing.risks2 <- function(object,times,cause,mr){

    stopifnot(length(times) == 1)

    event.time.order <- order(mr[,"time"],-1.0 * (mr[,"status"] != 0))
    smary <- summary(object,times=times)
    F <- smary$pstate[, match(cause, smary$states)]

        Fk <- leaveOneOut.competing.risks2(object,times,cause=cause,mr)
        N <- length(Fk)
        Jk <- N*F-(N-1)*Fk
        Jk[order(event.time.order)]

}

#' Compute leave one out jackknife contributions of the survival function
#'
#' For each subject, the survival function is recomputed leaving that subject out.
#' This is the workhorse for \link{jackknife.survival2} and will generally not
#' be called by the user.
#'
#' @param object A survfit object, with a single event (no competing risks)
#' @param times Times at which the survival is computed, must be length 1
#' @param mr Model response, the result of a call to Surv, or a matrix with two columns: "time" (observed follow up time) and "status" (0 = censored, 1 = event)
#' @return A vector of jackknifed values of survival at time times
#' @export
#' @examples
#'
#' sfit.surv <- survival::survfit(survival::Surv(time, status) ~ 1, data = colon)
#' mrs <- with(colon, survival::Surv(time, status))
#' jackvals <- leaveOneOut.survival2(sfit.surv, times = 1000, mrs)


leaveOneOut.survival2 <- function(object,times,mr){
    stopifnot(length(times)==1)
    event.time.order <- order(mr[,"time"],-1.0*(mr[,"status"] != 0))
    mr <- mr[event.time.order,]

    time <- object$time
    Y <- object$n.risk
    D <- object$n.event
    Y <- Y[D>0]
    time <- time[D>0]
    D <- D[D>0]
    NU <- length(time)
    obstimes <- mr[,"time"]
    status <- mr[,"status"]
    N <- length(obstimes)
    ##
    #S <- predict(object,times=time,newdata=mr)
    ## idea: find the at-risk set for pseudo-value k by
    ##       substracting 1 in the period where subj k is
    ##       at risk. need the position of obstime.k in time ...

    tdex <- max(which(time <= times))

    loo <- .C("loo_surv2",
              Y = as.double(Y),
              D=as.double(D),
              time=as.double(time),
              obsT=as.double(obstimes),
              status=as.double(status),
              S=double(N),
              N=as.integer(N),
              NT=as.integer(NU),
              Tdex=as.integer(tdex - 1),
              PACKAGE="eventglm")$S
    loo
}


#' Compute jackknife pseudo-observations of the cause-specific cumulative incidence for competing risks
#'
#' @param object A survfit object, with competing events
#' @param times Times at which the cumulative incidence is computed, must be length 1
#' @param cause Value indicating for which cause the cumulative incidence is to be computed, it must match one of the values available in object (see example)
#' @param mr Model response, the result of a call to Surv, or a matrix with two columns: "time" (observed follow up time) and "status" (0 = censored, 1, ..., k = event types)
#' @return A vector of jackknifed values of the cause-specific cumulative incidence at time times
#' @export
#' @examples
#'
#' sfit.cuminc <- survival::survfit(survival::Surv(etime, event) ~ 1, data = mgus2)
#' mrs <- with(mgus2, survival::Surv(etime, event))
#' jackvals <- leaveOneOut.competing.risks2(sfit.cuminc, times = 200, cause = "pcm", mrs)


leaveOneOut.competing.risks2 <- function(object, times, cause, mr){
    stopifnot(length(times) == 1)

    event.time.order <- order(mr[,"time"],-1.0*(mr[,"status"] != 0))
    mr <- mr[event.time.order,]
    states <- object$states
    if (missing(cause)) {
        C <- 2
        cause <- states[C]
    }   else{
        C <- match(cause,states,nomatch=0)
        if (length(C)>1 || C==0) stop("Cause must match exactly one of the names of object$n.event.")
    }
    D <- object$n.event[,C]
    #  it is sufficient to consider time points where events occur

    D0 <- rowSums( object$n.event)
    #  it is sufficient to consider time points where events occur
    time <- object$time[D0>0]
    Y <- object$n.risk[D0>0]

    D <- D[D0>0]
    D0 <- D0[D0>0]

    NU <- length(time)
    obstimes <- mr[,"time"]
    status <- as.numeric(mr[,"status"] != 0)
    E <- mr[, "status"] + 1
    N <- length(obstimes)
    tdex <- max(which(time <= times))

    loo2 <- .C("loo_comprisk2",
               Y = as.double(Y),
               D=as.double(D),
               D0 = as.double(D0),
               time=as.double(time),
               obsT=as.double(obstimes),
               status=as.double(status*(E==C)),
               status0=as.double(status),
               F=double(N),
               N=as.integer(N),
               NT=as.integer(NU),
               Tdex=as.integer(tdex - 1),
               PACKAGE="eventglm")$F


    loo2
}

#' Compute leave one out jackknife contributions of the survival function
#'
#' For each subject, the survival function is recomputed leaving that subject out.
#' This one does the calculation for all observed times, for calculation of the restricted mean
#'
#' @param object A survfit object, with a single event (no competing risks)
#' @param times Compute values at observed times up to and including this time
#' @param mr Model response, the result of a call to Surv, or a matrix with two columns: "time" (observed follow up time) and "status" (0 = censored, 1 = event)
#' @return A vector of jackknifed values of survival at time times
#' @export
#' @examples
#'
#' sfit.surv <- survival::survfit(survival::Surv(time, status) ~ 1, data = colon)
#' mrs <- with(colon, survival::Surv(time, status))
#' jackvals <- leaveOneOut.survival(sfit.surv, 1000, mrs)


leaveOneOut.survival <- function(object, times, mr){

    event.time.order <- order(mr[,"time"],-1.0*(mr[,"status"] != 0))
    mr <- mr[event.time.order,]

    time <- object$time

    tdex <- max(which(time <= times))

    Y <- object$n.risk
    D <- object$n.event
# Y <- Y[D>0]
#    time <- time[D>0]
#    D <- D[D>0]
    NU <- length(time)
    obstimes <- mr[,"time"]
    status <- mr[,"status"]
    N <- length(obstimes)
    ##
    #S <- predict(object,times=time,newdata=mr)
    ## idea: find the at-risk set for pseudo-value k by
    ##       substracting 1 in the period where subj k is
    ##       at risk. need the position of obstime.k in time ...

    loo <- .C("loo_surv",
              Y = as.double(Y),
              D=as.double(D),
              time=as.double(time),
              obsT=as.double(obstimes),
              status=as.double(status),
              S=double(N * tdex),
              N=as.integer(N),
              NT=as.integer(NU),
              Tdex=as.integer(tdex-1),
              PACKAGE="eventglm")$S


    out <- matrix(loo, nrow = N, byrow = FALSE)
    out[order(event.time.order),]
}



#' Compute jackknife pseudo-observations of the cause-specific cumulative incidence for competing risks
#'
#' This version computes them for all times up to \code{times}, for the restricted mean lifetime lost
#'
#' @param object A survfit object, with competing events
#' @param times Times at which the cumulative incidence is computed, must be length 1
#' @param cause Value indicating for which cause the cumulative incidence is to be computed, it must match one of the values available in object (see example)
#' @param mr Model response, the result of a call to Surv, or a matrix with two columns: "time" (observed follow up time) and "status" (0 = censored, 1, ..., k = event types)
#' @return A vector of jackknifed values of the cause-specific cumulative incidence at time times
#' @export
#' @examples
#'
#' sfit.cuminc <- survival::survfit(survival::Surv(etime, event) ~ 1, data = mgus2)
#' mrs <- with(mgus2, survival::Surv(etime, event))
#' jackvals <- leaveOneOut.competing.risks(sfit.cuminc, times = 200, cause = "pcm", mrs)


leaveOneOut.competing.risks <- function(object, times, cause, mr){
   # stopifnot(length(times) == 1)

    event.time.order <- order(mr[,"time"],-1.0*(mr[,"status"] != 0))
    mr <- mr[event.time.order,]
    states <- object$states[-1]
    if (missing(cause)) {
        C <- 1
        cause <- states[C]
    }   else{
        C <- match(cause,states,nomatch=0)
        if (length(C)>1 || C==0) stop("Cause must match exactly one of the names of object$n.event.")
    }
    D <- object$n.event[,C + 1]
    #  it is sufficient to consider time points where events occur

    D0 <- rowSums( object$n.event)
    #  it is sufficient to consider time points where events occur
    time <- object$time#[D0>0]
    Y <- object$n.risk#[D0>0]

    tdex <- max(which(time <= times))

    D <- D#[D0>0]
    D0 <- D0#[D0>0]

    NU <- length(time)
    obstimes <- mr[,"time"]
    status <- as.numeric(mr[,"status"] != 0)
    E <- mr[, "status"]
    N <- length(obstimes)

    loo2 <- .C("loo_comprisk",
               Y = as.double(Y),
               D=as.double(D),
               D0 = as.double(D0),
               time=as.double(time),
               obsT=as.double(obstimes),
               status=as.double(status*(E==C)),
               status0=as.double(status),
               F=double(N * tdex),
               N=as.integer(N),
               NT=as.integer(NU),
               Tdex=as.integer(tdex - 1),
               PACKAGE="eventglm")$F


    out <- matrix(loo2, nrow = N, byrow = FALSE)
    out[order(event.time.order),]
}
