
jackknife.survival2 <- function(object,times,keepResponse=FALSE,...){
    S <- predict(object,times=times,newdata=object$model.response)
    Sk <- leaveOneOut.survival(object,times,...)
    N <- NROW(Sk)
    Jk <- t(N*S-t((N-1)*Sk))
    colnames(Jk) <- paste("t",times,sep=".")
    if (keepResponse==TRUE){
        Jk <- cbind(object$model.response,Jk)
    }
    ## re-order the pseudo-values
    Jk <- Jk[object$originalDataOrder,,drop=FALSE]
    Jk
}

jackknife.competing.risks2 <- function(object,times,cause,mr){

    stopifnot(length(times) == 1)

    event.time.order <- order(mr[,"time"],-1.0 * (mr[,"status"] != 0))
    smary <- summary(object,times=times)
    F <- smary$pstate[, match(cause, smary$states)]

        Fk <- leaveOneOut.competing.risks2(object,times,cause=cause,mr)
        N <- length(Fk)
        Jk <- N*F-(N-1)*Fk
        matrix(Jk[order(event.time.order)], ncol = 1)

}

leaveOneOut.survival2 <- function(object,times,lag=FALSE,...){
    stopifnot(object$covariate.type==1)
    mr <- object$model.response
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
    S <- predict(object,times=time,newdata=mr)
    ## idea: find the at-risk set for pseudo-value k by
    ##       substracting 1 in the period where subj k is
    ##       at risk. need the position of obstime.k in time ...
    ## pos <- match(obstimes,time)
    ## if (useC==TRUE){
    loo <- .C("loo_surv",
              Y = as.double(Y),
              D=as.double(D),
              time=as.double(time),
              obsT=as.double(obstimes),
              status=as.double(status),
              S=double(NU*N),
              N=as.integer(N),
              NT=as.integer(NU),
              PACKAGE="prodlim")$S
    out <- matrix(loo,nrow=N,ncol=NU,byrow=FALSE)
    ## }
    ## else{
    pos <- sindex(jump.times=time,eval.times=obstimes)
    ## loo2 <- do.call("rbind",lapply(1:N,function(k){
    ## Dk <- D
    ## if (status[k]==1) Dk[pos[k]] <- Dk[pos[k]]-1
    ## Yk <- Y-c(rep(1,pos[k]),rep(0,NU-pos[k]))
    ## cumprod(1-Dk/Yk)}))
    ## }
    ## out <- loo
    if (!missing(times)){
        found <- sindex(jump.times=time,eval.times=times)+1
        if (lag==FALSE)
            out <- cbind(1,out)[,found,drop=TRUE]
        else
            out <- cbind(1,cbind(1,out))[,found,drop=TRUE]
    }
    out
}

leaveOneOut.competing.risks2 <- function(object, times, cause, mr){
    stopifnot(length(times) == 1)

    event.time.order <- order(mr[,"time"],-1.0*(mr[,"status"] != 0))
    mr <- mr[event.time.order,]
    states <- object$states
    if (missing(cause)) {
        C <- 1
        cause <- states[1]
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
    ## browser()
    ## }
    ## else{
    ## pos <- sindex(jump.times=time,eval.times=obstimes)
    ## loo <- do.call("rbind",lapply(1:N,function(k){
    ## Dk <- D
    ## if (status[k]==1 && E[k]==cause) Dk[pos[k]] <- Dk[pos[k]]-1
    ## Yk <- Y-c(rep(1,pos[k]),rep(0,NU-pos[k]))
    ## Sk <- as.numeric(lagSk[k,,drop=TRUE])
    ## Hk <- Dk/Yk
    ## Fk <- cumsum(Sk*Hk)
    ## Fk
    ## }))
    ## out <- loo
    ## }

    loo2
}
