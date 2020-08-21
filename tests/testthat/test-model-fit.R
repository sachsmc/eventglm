test_that("Residuals work", {
    cuminctest <- cumincglm(survival::Surv(etime, event) ~ sex,
                            time = 200, cause = "pcm", link = "identity", mgus2)

    cuminclog <- cumincglm(survival::Surv(etime, event) ~ age,
                            time = 200, cause = "pcm", link = "cloglog", mgus2)

    t1 <- residuals(cuminclog)
    t2 <- residuals.glm(cuminclog)
    expect_true(!all(t1 == t2))

    #plot(t1 ~ mgus2$age)

    survtest <- cumincglm(survival::Surv(time, status) ~ age,
                          time = 1000, link = "identity", data = colon)

    s1 <- residuals(survtest)
    s2 <- residuals.glm(survtest)
    expect_true(!all(s1 == s2))

    #plot(s1 ~ colon$age)


    ## restricted mean should be no transformation

})
