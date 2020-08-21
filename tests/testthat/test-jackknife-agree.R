test_that("Jackknife agrees with prodlim", {
    library(prodlim)
    colon <- survival::colon
    colon$event <- factor(ifelse(colon$status == 1, colon$etype, 0),
                          levels = c(0, 1, 2),
                          labels = c("censored", "recurrence", "death"))
    sfit <- survival::survfit(survival::Surv(time, event) ~ 1, data = colon)
    marginal.estimate <- prodlim::prodlim(prodlim::Hist(time, as.numeric(event) - 1) ~ 1, data = colon)
    jackk2 <- prodlim:::jackknife.competing.risks(marginal.estimate, times = 1000, cause = "1")

    mr <- with(colon, survival::Surv(time, event))
    myest <- eventglm:::jackknife.competing.risks2(sfit, times = 1000, cause = "recurrence", mr)

    expect_true(all(myest == jackk2))


    sfit.surv <- survival::survfit(survival::Surv(time, status) ~ 1, data = colon)
    me.surv <- prodlim::prodlim(prodlim::Hist(time, status) ~ 1, data = colon)
    jack.s <- prodlim:::jackknife.survival(me.surv, times = 1000)

    mrs <- with(colon, survival::Surv(time, status))
    myests <- eventglm:::jackknife.survival2(sfit.surv, times = 1000, mrs)

    expect_true(all(myests == jack.s))

})
