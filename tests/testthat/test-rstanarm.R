## test tidy and glance methods from rstanarm_tidiers.R
stopifnot(require("testthat"), require("broom.mixed"))
context("rstanarm tidiers")

if (suppressPackageStartupMessages(require(rstanarm, quietly = TRUE))) {
  load(system.file("extdata", "rstanarm_example.rda", package = "broom.mixed"))

  test_that("tidy works on rstanarm fits", {
    td1 <- tidy(fit)
    td2 <- tidy(fit, effects = "ran_vals")
    td3 <- tidy(fit, effects = "ran_pars")
    td4 <- tidy(fit, effects = "auxiliary")

    check_tidy(td1, exp.names = c("group", "term", "estimate", "std.error"))
    check_tidy(td2, exp.names = c("group", "level", "term", "estimate", "std.error"))
    check_tidy(td3, exp.names = c("group", "term", "estimate"))
    check_tidy(td4, exp.names = c("term", "estimate", "std.error"))
  })

  test_that("tidy with multiple 'effects' selections works on rstanarm fits", {
    td1 <- tidy(fit, effects = c("ran_vals", "auxiliary"))
    expect_true(all(c("sigma", "mean_PPD") %in% td1$term))
    check_tidy(td1, exp.names = c("group", "level", "term", "estimate", "std.error"))
  })

  test_that("conf.int works on rstanarm fits", {
    td1 <- tidy(fit, conf.int = TRUE, conf.level = 0.8)
    td2 <- tidy(fit, effects = "ran_vals", conf.int = TRUE, conf.level = 0.5)
    td3 <- tidy(fit, conf.int = TRUE, conf.level = 0.95)

    check_tidy(td2, exp.names = c("group", "level", "term", "estimate",
                                   "std.error", "conf.low", "conf.high"))

    ## Wider conf.level should give wider intervals (where not NA)
    valid_idx <- !is.na(td3$conf.low) & !is.na(td1$conf.low)
    if (any(valid_idx)) {
      expect_true(all(td3$conf.low[valid_idx] <= td1$conf.low[valid_idx]))
      expect_true(all(td3$conf.high[valid_idx] >= td1$conf.high[valid_idx]))
    }
  })

  test_that("glance works on rstanarm fits", {
    skip_on_cran()
    g1 <- glance(fit)
    g2 <- glance(fit, looic = TRUE, cores = 1, k_threshold = 0.7)

    check_tidy(g1, exp.row = 1, exp.names = c("algorithm", "pss", "nobs", "sigma"))
    check_tidy(g2, exp.names = c("algorithm", "pss", "nobs", "sigma",
                                  "looic", "elpd_loo", "p_loo"))
  })

  test_that("exponentiation works", {
    td1 <- tidy(fit2, conf.int = TRUE, effects = "fixed")
    td1e <- tidy(fit2, conf.int = TRUE, exponentiate = TRUE, effects = "fixed")

    expect_equal(td1e$estimate, exp(td1$estimate))
    expect_equal(td1e$conf.low, exp(td1$conf.low))
    expect_equal(td1e$conf.high, exp(td1$conf.high))
    expect_equal(td1e$std.error, exp(td1$estimate) * td1$std.error)
  })

  ## GH 153
  test_that("mice imputed data works", {
    skip_if_not_installed("mice")
    data(nhanes, package = "mice")
    imp <- mice::mice(nhanes, m = 3, print = FALSE)
    suppressWarnings(
      capture.output(ms <- with(imp, stan_glm(age ~ bmi + chl)))
    )
    expect_s3_class(summary(mice::pool(ms)), "mipo.summary")
  })

  ## GH 156: Increase parallels between tidy.brmsfit and tidy.stanreg
  test_that("robust parameter works", {
    ## Default (robust = TRUE) uses median/mad
    td_robust <- tidy(fit, effects = "fixed")
    ## Non-robust uses mean/sd
    td_nonrobust <- tidy(fit, effects = "fixed", robust = FALSE)

    ## Estimates should differ (median vs mean)
    expect_false(isTRUE(all.equal(td_robust$estimate, td_nonrobust$estimate)))
    ## Std errors should differ (mad vs sd)
    expect_false(isTRUE(all.equal(td_robust$std.error, td_nonrobust$std.error)))
    ## Same columns
    expect_equal(colnames(td_robust), colnames(td_nonrobust))
  })

  test_that("fix.intercept parameter works", {
    ## rstanarm already uses "(Intercept)" by default
    td_fix <- tidy(fit, effects = "fixed", fix.intercept = TRUE)
    td_nofix <- tidy(fit, effects = "fixed", fix.intercept = FALSE)

    ## With fix.intercept = TRUE (default), keep "(Intercept)"
    expect_true("(Intercept)" %in% td_fix$term)
    ## With fix.intercept = FALSE, convert to "Intercept" (no parentheses)
    expect_true("Intercept" %in% td_nofix$term)
    expect_false("(Intercept)" %in% td_nofix$term)
  })

  test_that("rhat and ess parameters work", {
    td_diag <- tidy(fit, effects = "fixed", rhat = TRUE, ess = TRUE)

    expect_true("rhat" %in% colnames(td_diag))
    expect_true("ess" %in% colnames(td_diag))
    ## Values should be numeric and not all NA
    expect_true(is.numeric(td_diag$rhat))
    expect_true(is.numeric(td_diag$ess))
  })

  test_that("robust works with ran_vals", {
    td_robust <- tidy(fit, effects = "ran_vals", robust = TRUE)
    td_nonrobust <- tidy(fit, effects = "ran_vals", robust = FALSE)

    ## Same structure
    expect_equal(colnames(td_robust), colnames(td_nonrobust))
    ## Different values
    expect_false(isTRUE(all.equal(td_robust$estimate, td_nonrobust$estimate)))
  })

  test_that("robust works with auxiliary", {
    td_robust <- tidy(fit, effects = "auxiliary", robust = TRUE)
    td_nonrobust <- tidy(fit, effects = "auxiliary", robust = FALSE)

    ## Same structure
    expect_equal(colnames(td_robust), colnames(td_nonrobust))
    ## Check required columns present
    expect_true(all(c("estimate", "std.error") %in% colnames(td_robust)))
  })

} ## rstanarm available
