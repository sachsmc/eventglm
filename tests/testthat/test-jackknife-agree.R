test_that("Jackknife agrees with prodlim", {
    library(prodlim)

    sfit <- survival::survfit(survival::Surv(etime, event) ~ 1, data = mgus2)
    marginal.estimate <- prodlim::prodlim(prodlim::Hist(etime, as.numeric(event) - 1) ~ 1, data = mgus2)
    jackk2 <- prodlim:::jackknife.competing.risks(marginal.estimate, times = 200, cause = "1")

    mr <- with(mgus2, survival::Surv(etime, event))
    myest <- eventglm:::jackknife.competing.risks2(sfit, times = 200, cause = "pcm", mr)

    expect_true(all(myest == jackk2))


    sfit.surv <- survival::survfit(survival::Surv(time, status) ~ 1, data = colon)
    me.surv <- prodlim::prodlim(prodlim::Hist(time, status) ~ 1, data = colon)
    jack.s <- prodlim:::jackknife.survival(me.surv, times = 1000)

    mrs <- with(colon, survival::Surv(time, status))
    myests <- eventglm:::jackknife.survival2(sfit.surv, times = 1000, mrs)

    expect_true(all(myests == jack.s))

})
