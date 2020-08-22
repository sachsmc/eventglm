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
