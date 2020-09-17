test_that("Residuals work", {
    cuminctest <- cumincglm(survival::Surv(etime, event) ~ sex,
                            time = 200, cause = "pcm", link = "identity", data = mgus2)

    cuminclog <- cumincglm(survival::Surv(etime, event) ~ age,
                            time = 200, cause = "pcm", link = "cloglog", data = mgus2)

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

    rmeantest <- rmeanglm(survival::Surv(time, status) ~ age,
                          time = 1000, link = "identity", data = colon)

    expect_true(is.numeric(confint(rmeantest)))

    expect_error(vcov(rmeantest, type = "corrected"))

    rmeantest2 <- rmeanglm(survival::Surv(etime, event) ~ sex,
                          time = 200, cause = "pcm", link = "identity", data = mgus2)

})

test_that("ipcw works", {


    cumincipcw <- cumincglm(survival::Surv(etime, event) ~ age + sex,
                            time = 200, cause = "pcm", link = "identity",
                            model.censoring = "independent", data = mgus2)

    cumincipcw2 <- cumincglm(survival::Surv(etime, event) ~ age + sex,
                            time = 200, cause = "pcm", link = "identity",
                            model.censoring = "stratified",
                            formula.censoring = ~ sex, data = mgus2)

    cumincipcw3 <- cumincglm(survival::Surv(etime, event) ~ age + sex,
                             time = 200, cause = "pcm", link = "identity",
                             model.censoring = "aareg",
                             data = mgus2)
    cumincipcw4 <- cumincglm(survival::Surv(etime, event) ~ age + sex,
                             time = 200, cause = "pcm", link = "identity",
                             model.censoring = "coxph",
                             data = mgus2)

    expect_true(cumincipcw$coefficients[3] - cumincipcw2$coefficients[3] < 1e-3)
    expect_true(cumincipcw$coefficients[3] - cumincipcw3$coefficients[3] < 1e-3)
    expect_true(cumincipcw$coefficients[3] - cumincipcw4$coefficients[3] < 1e-3)


    rmeanipcw1 <- rmeanglm(survival::Surv(etime, event) ~ age + sex,
                             time = 200, cause = "pcm", link = "identity",
                             model.censoring = "independent",
                             formula.censoring = ~ sex, data = mgus2)

    rmeanipcw2 <- rmeanglm(survival::Surv(etime, event) ~ age + sex,
                             time = 200, cause = "pcm", link = "identity",
                             model.censoring = "stratified",
                             formula.censoring = ~ sex, data = mgus2)

    rmeanipcw3 <- rmeanglm(survival::Surv(etime, event) ~ age + sex,
                           time = 200, cause = "pcm", link = "identity",
                           model.censoring = "aareg",
                           formula.censoring = ~ age + sex, data = mgus2)

    rmeanipcw4 <- rmeanglm(survival::Surv(etime, event) ~ age + sex,
                           time = 200, cause = "pcm", link = "identity",
                           model.censoring = "coxph",
                           formula.censoring = ~ age + sex, data = mgus2)

    expect_true(rmeanipcw1$coefficients[3] - rmeanipcw2$coefficients[3] < 1e-2)
    expect_true(rmeanipcw1$coefficients[3] - rmeanipcw3$coefficients[3] < 1e-2)
    expect_true(rmeanipcw1$coefficients[3] - rmeanipcw4$coefficients[3] < 1e-2)


    expect_error(vcov(rmeanipcw1, type = "corrected"))
    expect_true(all(sqrt(diag(vcov(rmeanipcw1, type = "naive"))) > 0))
    expect_true(all(sqrt(diag(vcov(rmeanipcw1, type = "robust"))) > 0))

})
