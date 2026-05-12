stopifnot(require("testthat"), require("broom.mixed"))
context("brms tidiers")

if (require(brms, quietly = TRUE)) {
  load(system.file("extdata", "brms_example.rda",
    package = "broom.mixed",
    mustWork = TRUE
  ))

  test_that("tidy returns conf.int columns when requested (GH #87)", {
    skip_on_cran()
    tt <- suppressWarnings(tidy(brms_multi, conf.int = TRUE))
    expect_true(all(c("conf.low", "conf.high") %in% names(tt)))
  })

  test_that("glance returns expected columns (GH #101)", {
    gg <- glance(brms_noran)
    check_tidy(gg, exp.row = 1, exp.names = c("algorithm", "pss", "nobs", "sigma"))
  })

  test_that("tidy returns correct structure for model with random effects", {
    skip_on_cran()
    td <- suppressWarnings(tidy(brms_RE))
    check_tidy(td, exp.row = 6, exp.col = 8,
               exp.names = c("effect", "component", "group", "term",
                             "estimate", "std.error"))

    ## Check effect/group/term values
    expect_equal(td$effect, c("fixed", "fixed", rep("ran_pars", 4)))
    expect_equal(td$group[1:2], c(NA_character_, NA_character_))
    expect_equal(td$group[3:6], c("Subject", "Subject", "Subject", "Residual"))
    expect_true("(Intercept)" %in% td$term)
    expect_true("sd__Observation" %in% td$term)
  })

  test_that("tidy returns correct structure for model without random effects", {
    skip_on_cran()
    td <- suppressWarnings(tidy(brms_noran))
    check_tidy(td, exp.row = 3,
               exp.names = c("effect", "component", "group", "term"))

    expect_equal(td$effect, c("fixed", "fixed", "ran_pars"))
    expect_equal(td$group, c(NA_character_, NA_character_, "Residual"))
  })

  test_that("tidy returns correct structure for fixed-effects only model", {
    skip_on_cran()
    td <- suppressWarnings(tidy(brms_brm_fit4))
    check_tidy(td, exp.row = 2,
               exp.names = c("effect", "component", "term"))

    expect_equal(td$effect, c("fixed", "fixed"))
    expect_true(all(is.na(td$group)))
    expect_true("(Intercept)" %in% td$term)
  })

  test_that("tidy errors on invalid effects argument", {
    expect_error(suppressWarnings(tidy(brms_multi, effects = "junk")),
                 "unrecognized effects")
  })

  test_that("tidy respects effects argument", {
    skip_on_cran()
    td_fixed <- suppressWarnings(tidy(brms_RE, effects = "fixed"))
    td_ran_pars <- suppressWarnings(tidy(brms_RE, effects = "ran_pars"))

    expect_true(all(td_fixed$effect == "fixed"))
    expect_true(all(td_ran_pars$effect == "ran_pars"))
  })

  test_that("tidy respects conf.level argument", {
    skip_on_cran()
    td_95 <- suppressWarnings(tidy(brms_RE, conf.int = TRUE, conf.level = 0.95))
    td_50 <- suppressWarnings(tidy(brms_RE, conf.int = TRUE, conf.level = 0.50))

    ## Narrower conf.level should give narrower intervals
    expect_true(all(td_50$conf.low >= td_95$conf.low, na.rm = TRUE))
    expect_true(all(td_50$conf.high <= td_95$conf.high, na.rm = TRUE))
  })

  test_that("robust parameter affects estimates", {
    skip_on_cran()
    td_robust <- suppressWarnings(tidy(brms_RE, effects = "fixed", robust = TRUE))
    td_nonrobust <- suppressWarnings(tidy(brms_RE, effects = "fixed", robust = FALSE))

    ## Same structure, different values
    expect_equal(colnames(td_robust), colnames(td_nonrobust))
    expect_false(isTRUE(all.equal(td_robust$estimate, td_nonrobust$estimate)))
  })

} ## if require(brms)
