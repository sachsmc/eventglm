
#' Print method for pseudoglm
#' @param x A pseudoglm object, as returned by \link{cumincglm} or \link{rmeanglm}
#' @param digits Number of significant digits
#' @param ... Not used
#' @return x, invisibly
#'
#' @export
#'
print.pseudoglm <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
    outcome <- switch(x$type, rmean = "restricted mean",
                      cuminc = "cumulative incidence",
                      survival = "survival")
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


#' Compute covariance matrix of regression coefficient estimates
#'
#'
#' @param object A pseudoglm object, as returned by \link{cumincglm} or
#'   \link{rmeanglm}.
#' @param type The method to use for variance estimation; one of "corrected",
#'   "robust", "naive", or "cluster"
#' @param ... Arguments passed to \link[sandwich]{vcovHC} if type = "robust",
#'   or to \link[sandwich]{vcovCL} if type = "cluster"
#' @return A numeric matrix containing the variance-covariance estimates
#'
#' @details The "corrected" variance estimate is as described in Overgaard et
#'   al. (2017) <doi:10.1214/16-AOS1516>, with code adapted from Overgaard's
#'   Stata program. This method does not handle ties and only has
#'   marginal benefits in reasonable sample sizes. The default is "robust" which
#'   uses a sandwich estimator as implemented in the sandwich package. "cluster"
#'   is another option if you have clustered observations. Finally "naive" uses
#'   the same method as glm to compute the variance, and is known to be
#'   anti-conservative. The bootstrap is another recommended option that can be
#'   implemented using other tools; there is an example in the vignette.
#'
#' @references Overgaard, Morten; Parner, Erik Thorlund; Pedersen, Jan.
#'   Asymptotic theory of generalized estimating equations based on jack-knife
#'   pseudo-observations. Ann. Statist. 45 (2017), no. 5, 1988--2015.
#'   <doi:10.1214/16-AOS1516>.
#' @seealso \link[sandwich]{vcovHC}
#' @export
vcov.pseudoglm <- function(object, type = "robust", ...) {

    if(type == "corrected") {

        if(is.null(object$x)) {
            stop("Corrected variance requires 'x = TRUE' in the model fit.")
        }
        datamat <- object$datamat
        noobs <- len <- nrow(datamat)

        datord <- order(datamat[, 1], -datamat[, 2], -datamat[, 3], -datamat[, 4])
        datamat <- datamat[datord, ]

        datain <- datamat[datamat[, 1] <= object$time, ]

        if(any(table(datain[, 1]) > 1)) {
            stop("Corrected variance not available when there are tied event times")
        }

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
        sandwich::vcovHC(object, ...)

    } else if(type == "naive") {

        class(object) <- c("glm", "lm")
        stats::vcov(object)


    } else if(type == "cluster") {

        class(object) <- c("glm", "lm")
        sandwich::vcovCL(object, cluster = object$cluster.id, ...)

    } else stop("unknown variance type")

}


#' Summary method
#'
#' @param object A pseudoglm object, as returned by \link{cumincglm} or \link{rmeanglm}
#' @param correlation logical; if TRUE, the correlation matrix of the estimated parameters is returned and printed.
#' @param symbolic.cor logical; If TRUE, print the correlations in a symbolic form rather than as numbers.
#' @param type The method to use for variance estimation; one of "corrected", "robust", "naive", or "cluster"
#' @param ... Additional arguments passed to \link{vcov.pseudoglm}
#' @return An object of class \link[stats]{summary.glm}
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
                                                                          df.r, df.f),
                                cov.unscaled = covmat.unscaled, cov.scaled = covmat))
    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.glm"
    return(ans)
}

#' Pseudo-observation scaled residuals
#'
#' Computes residuals according to the recommendations of Pohar-Perme and
#' Andersen (2009) <doi: 10.1002/sim.3401>.
#'
#' @param object A pseudoglm object, as returned by \link{cumincglm} or
#'   \link{rmeanglm}
#' @param type Either "scaled" (the default for cumulative incidence outcomes)
#'   or one of the types available in \link[stats]{residuals.glm} for restricted mean outcomes, with the default being "deviance".
#' @param ... Arguments passed on to \link[stats]{residuals.glm}.
#' @return A numeric vector of residuals
#'
#' @details The scaled residuals are computed as \deqn{\hat{\epsilon}_i =
#'   \frac{\hat{E}(V_i) - \hat{Y}_i}{\sqrt{\hat{Y}_i (1 - \hat{Y}_i)}}} When the
#'   outcome is the cumulative incidence, the denominator corresponds to an
#'   estimate of the standard error of the conditional estimate of the outcome
#'   in the absence of censoring. For the restricted mean, no such rescaling is
#'   done and the computation is passed off to \link[stats]{residuals.glm}.
#' @references Perme MP, Andersen PK. Checking hazard regression models using
#'   pseudo-observations. Stat Med. 2008;27(25):5309-5328.
#'   <doi:10.1002/sim.3401>
#' @export
residuals.pseudoglm <- function(object, type = NULL, ...){

    if(object$type == "rmean") {
        if(is.null(type)) {
            type <- "deviance"
        }
        residuals.glm(object, type = type, ...)
    } else if(object$type == "cuminc") {

        if(is.null(type) || type == "scaled") {

            cond.cuminc <- predict(object, type = "response")
            (object$y - cond.cuminc) / sqrt(cond.cuminc * (1 - cond.cuminc))


        } else {
            residuals.glm(object, type = type, ...)
        }

    }

}


#' Confidence Intervals for pseudoglm Model Parameters
#'
#' Computes Wald confidence intervals for one or more parameters in a fitted model. Users can specify the type of variance estimate used, with the default being the robust sandwich variance estimator.
#'
#' @param object A fitted model object from \link{cumincglm} or \link{rmeanglm}
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param type The type of variance estimate to use, see \link{vcov.pseudoglm}
#' @param ... Not used
#'
#' @return A matrix (or vector) with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 in % (by default 2.5% and 97.5%).
#'
#' @export
#' @examples
#' cumincipcw <- cumincglm(survival::Surv(etime, event) ~ age + sex,
#'          time = 200, cause = "pcm", link = "identity",
#'          model.censoring = "independent", data = mgus2)
#' confint(cumincipcw)
#'
confint.pseudoglm <- function (object, parm, level = 0.95, type = "robust", ...)
{
    cf <- coef(object)
    ses <- sqrt(diag(vcov(object, type = type)))
    pnames <- names(ses)
    if (is.matrix(cf)){
        cf <- setNames(as.vector(cf), pnames)
    }
    if (missing(parm)) {
        parm <- pnames
    } else if (is.numeric(parm)){
        parm <- pnames[parm]
    }
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    fac <- stats::qnorm(a)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
    ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm,
                                                                     pct))
    ci[] <- cf[parm] + ses[parm] %o% fac
    ci
}



