test_that("Residuals work", {
    cuminctest <- cumincglm(Surv(etime, event) ~ sex,
                            time = 200, cause = "pcm", link = "identity", data = mgus2)

    cuminclog <- cumincglm(Surv(etime, event) ~ age,
                            time = 200, cause = "pcm", link = "cloglog", data = mgus2)

    t1 <- residuals(cuminclog)
    t2 <- residuals.glm(cuminclog)
    expect_true(!all(t1 == t2))

    #plot(t1 ~ mgus2$age)

    survtest <- cumincglm(Surv(time, status) ~ age,
                          time = 1000, link = "identity", data = colon)

    s1 <- residuals(survtest)
    s2 <- residuals.glm(survtest)
    expect_true(!all(s1 == s2))

    #plot(s1 ~ colon$age)

    ## restricted mean should be no transformation

    rmeantest <- rmeanglm(Surv(time, status) ~ age,
                          time = 1000, link = "identity", data = colon)

    expect_true(is.numeric(confint(rmeantest)))

    expect_error(vcov(rmeantest, type = "corrected"))

    rmeantest2 <- rmeanglm(Surv(etime, event) ~ sex,
                          time = 200, cause = "pcm", link = "identity", data = mgus2)

})

test_that("ipcw works", {


    cumincipcw <- cumincglm(Surv(etime, event) ~ age + sex,
                            time = 200, cause = "pcm", link = "identity",
                            model.censoring = "independent", data = mgus2)

    cumincipcw2 <- cumincglm(Surv(etime, event) ~ age + sex,
                            time = 200, cause = "pcm", link = "identity",
                            model.censoring = "stratified",
                            formula.censoring = ~ sex, data = mgus2)

    cumincipcw3 <- cumincglm(Surv(etime, event) ~ age + sex,
                             time = 200, cause = "pcm", link = "identity",
                             model.censoring = "aareg",
                             data = mgus2)
    cumincipcw4 <- cumincglm(Surv(etime, event) ~ age + sex,
                             time = 200, cause = "pcm", link = "identity",
                             model.censoring = "coxph",
                             data = mgus2)

    expect_true(cumincipcw$coefficients[3] - cumincipcw2$coefficients[3] < 1e-3)
    expect_true(cumincipcw$coefficients[3] - cumincipcw3$coefficients[3] < 1e-3)
    expect_true(cumincipcw$coefficients[3] - cumincipcw4$coefficients[3] < 1e-3)


    rmeanipcw1 <- rmeanglm(Surv(etime, event) ~ age + sex,
                             time = 200, cause = "pcm", link = "identity",
                             model.censoring = "independent",
                             formula.censoring = ~ sex, data = mgus2)

    rmeanipcw2 <- rmeanglm(Surv(etime, event) ~ age + sex,
                             time = 200, cause = "pcm", link = "identity",
                             model.censoring = "stratified",
                             formula.censoring = ~ sex, data = mgus2)

    rmeanipcw3 <- rmeanglm(Surv(etime, event) ~ age + sex,
                           time = 200, cause = "pcm", link = "identity",
                           model.censoring = "aareg",
                           formula.censoring = ~ age + sex, data = mgus2)

    rmeanipcw4 <- rmeanglm(Surv(etime, event) ~ age + sex,
                           time = 200, cause = "pcm", link = "identity",
                           model.censoring = "coxph",
                           formula.censoring = ~ age + sex, data = mgus2)

    expect_true(rmeanipcw1$coefficients[3] - rmeanipcw2$coefficients[3] < 1e-2)
    expect_true(rmeanipcw1$coefficients[3] - rmeanipcw3$coefficients[3] < 1e-2)
    expect_true(rmeanipcw1$coefficients[3] - rmeanipcw4$coefficients[3] < 1e-2)


    expect_error(vcov(rmeanipcw1, type = "corrected"))
    expect_true(all(sqrt(diag(vcov(rmeanipcw1, type = "naive"))) > 0))
    expect_true(all(sqrt(diag(vcov(rmeanipcw1, type = "robust"))) > 0))

    mgus2$age[1:20] <- NA

    expect_error(cumincglm(Surv(etime, event) ~ age + sex,
                             time = 200, cause = "pcm", link = "identity",
                             model.censoring = "coxph",
                             data = mgus2))

    expect_error(cumincglm(Surv(etime, event) ~ age + sex,
                           time = 200, cause = "pcm", link = "identity",
                           model.censoring = "aareg",
                           data = mgus2))

    mgus2$sex[1:20] <- NA
    expect_error(cumincglm(Surv(etime, event) ~ sex,
                           time = 200, cause = "pcm", link = "identity",
                           model.censoring = "stratified",
                           formula.censoring = ~ sex,
                           data = mgus2))



})


test_that("variable names clash", {

    colon$pseudo.vals <- colon$surg
    colon$.Tci <- colon$surg
    colon$.Ci <- colon$surg

    goodest <- cumincglm(Surv(time, status) ~ rx + surg, time = 2500, data = colon)
    clash <- cumincglm(Surv(time, status) ~ rx + pseudo.vals, time = 2500, data = colon)

    expect_equal(unname(coef(goodest)), unname(coef(clash)))

    good2 <- cumincglm(Surv(time, status) ~ rx, time = 2500, data = colon,
                       model.censoring = "coxph", formula.censoring = ~ surg)
    clash2 <- cumincglm(Surv(time, status) ~ rx, time = 2500, data = colon,
                        model.censoring = "coxph", formula.censoring = ~ .Tci)

    expect_equal(unname(coef(good2)), unname(coef(clash2)))

    good3 <- cumincglm(Surv(time, status) ~ rx, time = 2500, data = colon,
                       model.censoring = "aareg", formula.censoring = ~ surg)
    clash3 <- cumincglm(Surv(time, status) ~ rx, time = 2500, data = colon,
                        model.censoring = "aareg", formula.censoring = ~ .Tci)

    expect_equal(unname(coef(good3)), unname(coef(clash3)))

})


test_that("Multiple times work", {

    cuminctest <- cumincglm(Surv(etime, event) ~ 1,
                            time = c(50, 100, 200), cause = "pcm", link = "identity", data = mgus2)
    stest <- survival::survfit(Surv(etime, event) ~ 1, data = mgus2)
    stab <- summary(stest, times = c(50, 100, 200))

    expect_lt(sum(coef(cuminctest)[1] + c(0, coef(cuminctest)[-1]) - stab$pstate[, 2]), 1e-6)

    cuminctest2 <- cumincglm(Surv(etime, event) ~ tve(sex),
                            time = c(50, 100, 200), cause = "pcm", link = "identity", data = mgus2)
    stest2 <- survival::survfit(Surv(etime, event) ~ sex, data = mgus2)
    stab2 <- summary(stest2, times = c(50, 100, 200))
    survests <- stab2$pstate[4:6, 2] - stab2$pstate[1:3, 2]
    expect_lt(sum(coef(cuminctest2)[4:6] - survests), 1e-4)

})

test_that("Missing data", {

  cuminctest <- cumincglm(Surv(etime, event) ~ hgb + creat + sex,
                   time = c(50, 100, 200), cause = "pcm", link = "identity", data = mgus2)

  expect_lt(cuminctest$df.residual, nrow(mgus2) * 3)

  cuminctest <- cumincglm(Surv(etime, event) ~ hgb + creat + sex,
                          time = c(50), cause = "pcm", link = "identity", data = mgus2)

  expect_lt(cuminctest$df.residual, nrow(mgus2))

  mgus2$idtest <- 1:nrow(mgus2)

  rmeantest <- rmeanglm(Surv(etime, event) ~ hgb + creat + sex,
                          time = c(50), cause = "pcm", link = "identity", data = mgus2,
                        id = idtest)
  expect_lt(rmeantest$df.residual, nrow(mgus2))


})



test_that("Clustered data works", {

  cuminctest <- cumincglm(Surv(etime, event) ~ 1,
                          time = c(50, 100, 200), cause = "pcm", link = "identity",
                          id = id, data = mgus2)
  stest <- survival::survfit(Surv(etime, event) ~ 1, data = mgus2)
  stab <- summary(stest, times = c(50, 100, 200))

  expect_lt(sum(coef(cuminctest)[1] + c(0, coef(cuminctest)[-1]) - stab$pstate[, 2]), 1e-6)

  cuminctest2 <- cumincglm(Surv(stop, event) ~ tve(rx), time = c(10, 25),
                           link = "identity", id = id, data = survival::bladder)
  rmeantest2 <- rmeanglm(Surv(stop, event) ~ rx, time = 25,
                           link = "identity", id = id, data = survival::bladder)

  expect_true(rmeantest2$method == "geese")


})


test_that("Glm features work", {

    set.seed(202105)
    colonx <- colon
    colonx$www <- runif(nrow(colonx))
    colonx$ooo <- rnorm(nrow(colonx))

    fitbas <- cumincglm(Surv(time, status) ~ rx, time = 2500, data = colonx,
                        model.censoring = "independent", formula.censoring = ~ surg)
    fit1 <- cumincglm(Surv(time, status) ~ rx + offset(ooo), time = 2500, data = colonx,
              model.censoring = "independent", formula.censoring = ~ surg)
    fit1b <- cumincglm(Surv(time, status) ~ rx, time = 2500, data = colonx,
                       offset = ooo,
                      model.censoring = "independent", formula.censoring = ~ surg)

    expect_true(!is.null(fit1$offset))
    expect_true(sum(abs(fitbas$coefficients - fit1$coefficients)) > .1)
    expect_true(sum(abs(fit1b$coefficients - fit1$coefficients)) < 1e-6)


    fit2 <- cumincglm(Surv(time, status) ~ rx, time = 2500, data = colonx,
                      weights = www,
                      model.censoring = "independent", formula.censoring = ~ surg)

    expect_true(all(fit2$weights != 1))
    expect_true(sum(abs(fitbas$coefficients - fit2$coefficients)) > .05)

    #multitime

    fitbas <- cumincglm(Surv(time, status) ~ tve(rx), time = c(500, 1000, 2500), data = colonx,
                        model.censoring = "independent", formula.censoring = ~ surg)

    fit1 <- cumincglm(Surv(time, status) ~ tve(rx) + offset(ooo),
                      time = c(500, 1000, 2500), data = colonx,
                      model.censoring = "independent", formula.censoring = ~ surg)

    fit1b <- cumincglm(Surv(time, status) ~ tve(rx), time =  c(500, 1000, 2500), data = colonx,
                       offset = ooo,
                       model.censoring = "independent", formula.censoring = ~ surg)

    expect_true(!is.null(fit1$offset))
    expect_true(sum(abs(fitbas$coefficients - fit1$coefficients)) > .05)
    expect_true(sum(abs(fit1b$coefficients - fit1$coefficients)) < .01)


    fit2 <- cumincglm(Surv(time, status) ~ tve(rx), time =  c(500, 1000, 2500), data = colonx,
                      weights = www,
                      model.censoring = "independent", formula.censoring = ~ surg)

    expect_true(all(fit2$weights != 1))
    expect_true(sum(abs(fitbas$coefficients - fit2$coefficients)) > .05)

    fitbass <- cumincglm(Surv(time, status) ~ tve(rx), time = c(500, 1000, 2500), data = colonx,
                        model.censoring = "independent", survival = TRUE, formula.censoring = ~ surg)


})



test_that("Left truncation and infjack", {

  library(survival)
  mdata <- tmerge(myeloid[!is.na(myeloid$crtime),1:2], myeloid[!is.na(myeloid$crtime),],
                  id=id, death= event(futime, death),
                  cr = event(crtime)
                  )

  mdata <- mdata[mdata$cr == 0,]

  sfit <- survfit(Surv(tstart, tstop, death) ~ trt, data = mdata,
                  influence = 1)

  sfitm <- summary(sfit, times = c(750))

  tinf <- cumincglm(Surv(tstart, tstop, death) ~ trt, data = mdata,
                    time = 750, model.censoring = "infjack", survival = TRUE)

  expect_true(abs(diff(sfitm$surv) - tinf$coefficients[2]) < .005)


  cuminctest <- cumincglm(Surv(time, status) ~ rx, 1500, cause = 1,
                          data = colon, survival = TRUE,
                          model.censoring = "stratified",
                          formula.censoring = ~ sex)

  cuminctest2 <- cumincglm(Surv(time, status) ~ rx, 1500, cause = 1,
                          data = colon, survival = TRUE,
                          model.censoring = "infjack",
                          formula.censoring = ~ sex)


  expect_true(mean(abs(cuminctest$coefficients - cuminctest2$coefficients)) < .01)


  cuminctest <- cumincglm(Surv(time, status) ~ rx, 1500, cause = 1,
                          data = colon, survival = TRUE,
                          model.censoring = "independent")

  cuminctest2 <- cumincglm(Surv(time, status) ~ rx, 1500, cause = 1,
                           data = colon, survival = TRUE,
                           model.censoring = "infjack",
                           formula.censoring = ~ 1)


  expect_true(mean(abs(cuminctest$coefficients - cuminctest2$coefficients)) < .005)


  cuminctest <- cumincglm(Surv(time, status) ~ rx, 1500, cause = 1,
                          data = colon, survival = FALSE,
                          model.censoring = "independent")

  cuminctest2 <- cumincglm(Surv(time, status) ~ rx, 1500, cause = 1,
                           data = colon, survival = FALSE,
                           model.censoring = "infjack",
                           formula.censoring = ~ 1)


  expect_true(mean(abs(cuminctest$coefficients - cuminctest2$coefficients)) < .005)

})


test_that("extension after last time", {

  set.seed(12) # it does not always, but quite occasionally happen, so I cherry-picked a seed. set.seed(11) does not reproduce the error, so it might work a second look.
  library(eventglm)

  arm <- sample(c('A', 'B'), 1000, replace=T)
  ev <- sample(c('Event', 'Cen', 'Comrisk'), 1000, replace=T) #competing-risk
  ttev <- rweibull(1000, 1, 1)  # time to event
  ttev <- ifelse(arm=='A' & ev=="Event", pmin(ttev, 2), ttev) # cens everyone of arm A and event "Event" as 2
  ttev <- ifelse(arm=='A' & ev!='Event', pmin(ttev, 1.8), ttev) # cens everyone else of arm A at 1.8. to ensure that 2 is last observed time for arm A and that is the event of interest.

  dt <- data.frame(arm=arm, ev=factor(ev, levels=c('Cen', 'Event', 'Comrisk')), ttev=ttev)

  expect_no_error(
    {
      prev.err <- rmeanglm(Surv(ttev,ev, type='mstate')~arm, time = 6, cause = 'Event', model.censoring = 'stratified',
           formula.censoring = ~ arm,
           data=dt) # Error in stats::glm.fit(X, Y, weights = weights, family = quasi(link = link,  : NA/NaN/Inf in 'y'
    }
  )



})
