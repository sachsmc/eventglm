test_that("Residuals work", {
    cuminctest <- cumincglm(survival::Surv(time, event) ~ rx,
                            time = 1000, cause = "recurrence", link = "identity", colon)

    survtest <- cumincglm(survival::Surv(time, status) ~ rx,
                          time = 1000, link = "identity", data = colon)
})
