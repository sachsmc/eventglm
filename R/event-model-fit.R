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
#' @param cause Character constant specifying the cause indicator of interest.
#' @param link Link function for the cumulative incidence regression model.
#' @param data Data frame in which all variables of formula can be interpreted.
#' @param ... Other arguments passed to \link[stats]{glm} such as weights, subset, etc.
#'
#' @export
#'
cumincglm <- function(formula, time, cause = "1", link = "identity", data, ...) {


    stopifnot(length(time) == 1)

    marginal.estimate <- survival::survfit(update.formula(formula, . ~ 1), data = data)

    newdata <- do.call(rbind, lapply(1:length(time), function(i) data))
    mr <- with(data, terms(update(formula, . ~ 1))[[2]])
    jackk <- jackknife.competing.risks2(marginal.estimate, times = time,
                                        cause = cause, mr)

    nn <- length(jackk)

    newdata[["pseudo.vals"]] <- c(jackk)
    newdata[["pseudo.time"]] <- rep(time, each = nn)

    if(marginal.estimate$model == "survival") newdata[["pseudo.vals"]] <- 1 - newdata[["pseudo.vals"]]

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
    fit.lin$type <- "cuminc"

    class(fit.lin) <- c("pseudoglm", class(fit.lin))

    fit.lin

}


#' Print method for pseudoglm
#'
#' @export
#'
print.pseudoglm <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
    outcome <- switch(x$type, rmean = "restricted mean", cuminc = "cumulative incidence")
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n",
                           collapse = "\n"), "\n\n", sep = "")
    cat("\nModel for the", x$link, outcome, "of cause", x$cause, "at time", x$time, "\n\n")
    if (length(coef(x))) {
        cat("Coefficients")
        if (is.character(co <- x$contrasts))
            cat("  [contrasts: ", apply(cbind(names(co),
                                              co), 1L, paste, collapse = "="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits),
                      print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
        x$df.residual, "Residual\n")
    if (nzchar(mess <- naprint(x$na.action)))
        cat("  (", mess, ")\n", sep = "")

    cat("\n")
    invisible(x)
}


#' Vcov method
#'
#' @param type one of "corrected", "robust", "naive", or "cluster"
#'
#' @export
vcov.pseudoglm <- function(object, type = "robust", ...) {

    if(type == "corrected") {

        datamat <- object$datamat
        noobs <- len <- nrow(datamat)

        datord <- order(datamat[, 1], -datamat[, 2], -datamat[, 3], -datamat[, 4])
        datamat <- datamat[datord, ]

        datain <- datamat[datamat[, 1] <= object$time, ]

        len2 <- nrow(datain)

        timejump <- datain[1:(len2 - 1), 1] != datain[2:len2, 1]
        times_id <- c(which(timejump == TRUE), len2)

        n_times <- length(times_id)
        times_nr <- cumsum(c(1, timejump))

        Y_lag <- (noobs:1)[times_id]

        N_all <- datain[, 2]
        N_1 <- datain[, 3]
        N_0 <- datain[, 4]

        # Calculating the overall survival function (missing values mean a population die out and should be interpreted as 0)
        H_0 <- cumsum(N_0)/noobs
        H_0_dif = c(H_0[times_id[1]],  ((H_0[times_id])[2:n_times] - (H_0[times_id])[1:(n_times-1)]))


        H_1 = cumsum(N_1)/noobs
        H_1_dif = c(H_1[times_id[1]] , ((H_1[times_id])[2:n_times] - (H_1[times_id])[1:(n_times-1)]))

        H_all = cumsum(N_all)/noobs
        H_all_dif = c(H_all[times_id[1]] , ((H_all[times_id])[2:n_times] - (H_all[times_id])[1:(n_times-1)]))

        H_lag = Y_lag / noobs
        H_lag_inv = 1 / H_lag
        H_lag_inv[is.na(H_lag_inv)] <- 0


        Lambda_0_dif = H_0_dif / H_lag
        Lambda_1_dif = H_1_dif / H_lag
        Lambda_all_dif = H_all_dif / H_lag
        Lambda_0_dif_comp_inv = 1 / (1 - Lambda_0_dif)
        Lambda_0_dif_comp_inv[is.na(Lambda_0_dif_comp_inv)] <- 0

        S = exp(cumsum(log(1 - Lambda_all_dif)))
        S[is.na(S)] <- 0
        S_lag = c(1 , S[1:(n_times-1)])

        G = exp(cumsum(log(1 - Lambda_0_dif)))
        G[is.na(G)] <- 0
        G_lag = c(1 , G[1:(n_times-1)])
        G_lag_inv = 1 / G_lag
        G_lag_inv[is.na(G_lag_inv)] <- 0

        F_1 = cumsum(S_lag * Lambda_1_dif)

        d_phi_1 = c((N_1 * G_lag_inv[times_nr]) , rep(0, noobs - len2))
        d_phi_2 = c(N_0 * (F_1[n_times] - F_1[times_nr]) * Lambda_0_dif_comp_inv[times_nr] * H_lag_inv[times_nr] -
                        cumsum(H_0_dif * (F_1[n_times] - F_1) * Lambda_0_dif_comp_inv * H_lag_inv^2)[times_nr] ,
                    rep(-sum(H_0_dif * (F_1[n_times] - F_1) * Lambda_0_dif_comp_inv * H_lag_inv^2), noobs-len2))

        beta = object$coefficients
        k = length(beta)

        z = object$x[datord,]  ## what if the model contains po at multiple time points?

        muhat = object$family$linkinv(z %*% beta)
        Ahat = z * as.vector(object$family$mu.eta(z %*% beta))
        Ahat_red = Ahat[1:len2,]
        mu_derivhat = z * as.vector(object$family$mu.eta(z %*% beta))

        a_len <- a1_len <- a2_len <- a3_len <-
            b_1 <- b_2 <- b_3 <- b_4 <- b_5 <- b_6 <-
            matrix(NA, nrow = noobs, ncol = k)
        H_0z <- H_0z_dif <- H_1z <- H_1z_dif <- H_z <- H_z_lag <- matrix(NA, nrow = len2, ncol = k)


        ## compute variance
        for(j in 1:k) {
            a_len[,j] = Ahat[,j] * (d_phi_1-muhat+d_phi_2)
            a1_len[,j] = Ahat[,j]  * (d_phi_1)
            a2_len[,j] = Ahat[,j] * (-muhat)
            a3_len[,j] = Ahat[,j] * (d_phi_2)

            H_0z[,j] = cumsum(Ahat_red[,j] * N_0)/noobs
            H_0z_dif[,j] = c(H_0z[times_id[1],j] , ((H_0z[times_id,j])[2:n_times] - (H_0z[times_id,j])[1:(n_times-1)]))

            H_1z[,j] = cumsum(Ahat_red[,j] * N_1)/noobs
            H_1z_dif[,j] = c(H_1z[times_id[1],j] , ((H_1z[times_id,j])[2:n_times] - (H_1z[times_id,j])[1:(n_times-1)]))

            H_z[,j] = (sum(Ahat[,j]) - cumsum(Ahat_red[,j]))/noobs
            H_z_lag[,j] = c(sum(Ahat[,j]) / noobs , (H_z[1:(len2-1),j])[times_id[1:(n_times-1)]])

            temp1 = cumsum( H_0z_dif[,j] * H_lag_inv * Lambda_0_dif_comp_inv - H_z_lag[,j] * H_lag_inv^2 * Lambda_0_dif_comp_inv * H_0_dif )
            temp1_lag = c(0 , temp1[1:(n_times-1)])
            temp1_dif = temp1 - temp1_lag

            temp2 = cumsum( H_0z_dif[,j] * (F_1[n_times] - F_1) * H_lag_inv * Lambda_0_dif_comp_inv -
                                H_z_lag[,j] * (F_1[n_times] - F_1) * H_lag_inv^2 * Lambda_0_dif_comp_inv  * H_0_dif )
            temp3 = cumsum( H_1z_dif[,j] / G_lag )


            b_1[,j] = c(N_1 * G_lag_inv[times_nr] * temp1_lag[times_nr] , rep(0, len-len2))

            b_2[,j] = c(N_0 * (temp3[n_times] - temp3[times_nr]) * H_lag_inv[times_nr] * Lambda_0_dif_comp_inv[times_nr] -
                            cumsum( H_0_dif * (temp3[n_times] - temp3) * H_lag_inv^2 * Lambda_0_dif_comp_inv )[times_nr] ,
                        rep(-sum( H_0_dif * (temp3[n_times] - temp3) * H_lag_inv^2 * Lambda_0_dif_comp_inv ), len-len2))

            temp4 = cumsum( H_1_dif * G_lag_inv * temp1 )

            b_3[,j] = c(N_0 * (temp4[n_times] - temp4[times_nr]) * H_lag_inv[times_nr] * Lambda_0_dif_comp_inv[times_nr] -
                            cumsum( H_0_dif * (temp4[n_times] - temp4) * H_lag_inv^2 * Lambda_0_dif_comp_inv )[times_nr] ,
                        rep(-sum( H_0_dif * (temp4[n_times] -temp4) / H_lag^2 / (1-Lambda_0_dif) ), len-len2 ))


            b_4[,j] = c(- N_0 * (F_1[n_times] - F_1[times_nr]) * H_z_lag[times_nr,j] * H_lag_inv[times_nr]^2 * Lambda_0_dif_comp_inv[times_nr] +
                            cumsum( H_0_dif * (F_1[n_times] - F_1) * H_z_lag[,j] * H_lag_inv^3 * Lambda_0_dif_comp_inv )[times_nr] ,
                        rep(sum( H_0_dif * (F_1[n_times] - F_1) * H_z_lag[,j] * H_lag_inv^3 * Lambda_0_dif_comp_inv ), len-len2 ))


            b_5[,j] = c(-cumsum( H_0z_dif[,j] * (F_1[n_times] - F_1) * H_lag_inv^2 * Lambda_0_dif_comp_inv -
                                     H_z_lag[,j] * (F_1[n_times] - F_1) * H_lag_inv^3 * Lambda_0_dif_comp_inv * H_0_dif )[times_nr] ,
                        rep(- sum( H_0z_dif[,j] * (F_1[n_times] - F_1) * H_lag_inv^2 * Lambda_0_dif_comp_inv -
                                       H_z_lag[,j] * (F_1[n_times] - F_1) * H_lag_inv^3 * Lambda_0_dif_comp_inv * H_0_dif ), len-len2 ))

            b_6[,j] = c(N_0 * (F_1[n_times] - F_1[times_nr]) * temp1_dif[times_nr] * H_lag_inv[times_nr] * Lambda_0_dif_comp_inv[times_nr] -
                            cumsum( H_0_dif * (F_1[n_times] - F_1) * temp1_dif * H_lag_inv^2 * Lambda_0_dif_comp_inv )[times_nr] ,
                        rep(-sum( H_0_dif * (F_1[n_times] - F_1) * temp1_dif * H_lag_inv^2 * Lambda_0_dif_comp_inv ), len-len2 ))


        }

        b <- b_1 + b_2 + b_3 + b_4 + b_5 + b_6
        Sigma <- t(a_len) %*% a_len / noobs + t(a_len) %*% b / noobs +
            t(b) %*% a_len / noobs + t(b) %*% b / noobs


        Minvhat <- solve(t(Ahat) %*% mu_derivhat / nrow(Ahat))

        Minvhat %*% Sigma %*% t(Minvhat) / noobs


    } else if(type == "robust") {

        class(object) <- c("glm", "lm")
        sandwich::sandwich(object)

    } else if(type == "naive") {

        class(object) <- c("glm", "lm")
        stats::vcov(object)


    } else if(type == "cluster") {

        class(object) <- c("glm", "lm")
        sandwich::vcovCL(object, cluster = object$cluster.id)

    } else stop("unknown variance type")

}


#' Summary method
#'
#'
#' @export
#'
summary.pseudoglm <- function (object, correlation = FALSE, symbolic.cor = FALSE,
                               type = "robust",
                               ...)
{
    df.r <- object$df.residual
    aliased <- is.na(coef(object))
    p <- object$rank
    if (p > 0) {
        p1 <- 1L:p
        Qr <- object$qr
        coef.p <- object$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- vcov(object, type = type)[p1, p1]
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err
        dn <- c("Estimate", "Std. Error")
        if (df.r > 0) {
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                          "z value", "Pr(>|z|)"))
        } else {
            coef.table <- cbind(coef.p, NaN, NaN, NaN)
            dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                          "z value", "Pr(>|z|)"))
        }
        df.f <- NCOL(Qr$qr)
    }
    else {
        coef.table <- matrix(, 0L, 4L)
        dimnames(coef.table) <- list(NULL, c("Estimate",
                                             "Std. Error", "z value", "Pr(>|z|)"))
        covmat.unscaled <- covmat <- matrix(, 0L, 0L)
        df.f <- length(aliased)
    }
    keep <- match(c("call", "terms", "family",
                    "deviance", "aic", "contrasts", "df.residual",
                    "null.deviance", "df.null", "iter",
                    "na.action"), names(object), 0L)
    ans <- c(object[keep], list(deviance.resid = residuals(object,
                                                           type = "deviance"), coefficients = coef.table,
                                aliased = aliased, dispersion = 1, df = c(object$rank,
                                                                          df.r, df.f), cov.unscaled = covmat.unscaled, cov.scaled = covmat))
    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.glm"
    return(ans)
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


