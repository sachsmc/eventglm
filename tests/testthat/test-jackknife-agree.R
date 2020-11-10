test_that("Jackknife agrees with prodlim", {
    library(prodlim)

    sfit <- survival::survfit(survival::Surv(etime, event) ~ 1, data = mgus2)
    marginal.estimate <- prodlim::prodlim(prodlim::Hist(etime, as.numeric(event) - 1) ~ 1, data = mgus2)
    jackk2 <- prodlim:::jackknife.competing.risks(marginal.estimate, times = 200, cause = "1")

    mr <- with(mgus2, survival::Surv(etime, event))
    myest <- eventglm:::jackknife.competing.risks2(sfit, times = 200, cause = "pcm", mr)

    head(cbind(myest, jackk2[, 1]))
    expect_true(all(abs(myest - jackk2[, 1]) < 1e-5))


    sfit.surv <- survival::survfit(survival::Surv(time, status) ~ 1, data = colon)
    me.surv <- prodlim::prodlim(prodlim::Hist(time, status) ~ 1, data = colon)
    jack.s <- prodlim:::jackknife.survival(me.surv, times = 1000)

    mrs <- with(colon, survival::Surv(time, status))
    myests <- eventglm:::jackknife.survival2(sfit.surv, times = 1000, mrs)

    expect_true(all(abs(myests - jack.s[, 1]) < 1e-5))


    ## restricted mean
    times <- sfit.surv$time[sfit.surv$time <= 1000]
    jack.s <- prodlim:::leaveOneOut.survival(me.surv, times = times)
    myests <- eventglm:::leaveOneOut.survival(sfit.surv, 1000, mrs)

    expect_true(all(dim(jack.s) == dim(myests)))

    expect_true(all(abs(jack.s - myests[order(mrs[,"time"],-1.0 * (mrs[,"status"] != 0)),]) < 1e-5))


    times <- sfit$time[sfit$time <= 200]
    jack.s2 <- prodlim:::leaveOneOut.competing.risks(marginal.estimate, times = times)
    myests2 <- eventglm:::leaveOneOut.competing.risks(sfit, 200, cause = "pcm", mr)

    expect_true(all(dim(jack.s2) == dim(myests2)))
    expect_true(all(abs(jack.s2 -
                            myests2[order(mr[,"time"],-1.0 * (mr[,"status"] != 0)),]) < 1e-5))

})
