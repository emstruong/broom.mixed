#' Tidying methods for a brms model
#'
#' These methods tidy the estimates from
#' \code{\link[brms:brmsfit-class]{brmsfit-objects}}
#' (fitted model objects from the \pkg{brms} package) into a summary.
#'
#' @return All tidying methods return a \code{data.frame} without rownames.
#' The structure depends on the method chosen.
#'
#' @seealso \code{\link[brms]{brms}}, \code{\link[brms]{brmsfit-class}}
#'
#' @name brms_tidiers
#'
#' @param x Fitted model object from the \pkg{brms} package. See
#'   \code{\link[brms]{brmsfit-class}}.
#' @examples
#'  ## original model
#'  \dontrun{
#'     brms_crossedRE <- brm(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
#'            iter = 500, chains = 2)
#'  }
#'  \donttest{
#'    ## too slow for CRAN (>5 seconds)
#'    ## load stored object
#'    if (require("rstan") && require("brms")) {
#'       load(system.file("extdata", "brms_example.rda", package="broom.mixed"))
#'
#'       fit <- brms_crossedRE
#'       tidy(fit)
#'       tidy(fit, parameters = "^sd_", conf.int = FALSE)
#'       tidy(fit, effects = "fixed", conf.method="HPDinterval")
#'       tidy(fit, effects = "ran_vals")
#'       tidy(fit, effects = "ran_pars", robust = TRUE)
#'       if (require("posterior")) {
#'       tidy(fit, effects = "ran_pars", rhat = TRUE, ess = TRUE)
#'
#'    }
#'    }
#'    if (require("rstan") && require("brms")) {
#'    # glance method
#'    glance(fit)
#'    ## this example will give a warning that it should be run with
#'    ## reloo=TRUE; however, doing this will fail
#'    ## because the \code{fit} object has been stripped down to save space
#'    suppressWarnings(glance(fit, looic = TRUE, cores = 1))
#'    head(augment(fit))
#'   }
#' }
#'
NULL
## examples for all methods (tidy/glance/augment) included in the same
##  block so we can surround them with a single "if (require(brms))" block

#' @rdname brms_tidiers
#' @param parameters Names of parameters for which a summary should be
#'   returned, as given by a character vector or regular expressions.
#'   If \code{NA} (the default) summarized parameters are specified
#'   by the \code{effects} argument.
#' @param effects A character vector including one or more of \code{"fixed"},
#'   \code{"ran_vals"}, or \code{"ran_pars"}.
#'   See the Value section for details.
#' @param robust Whether to use median and median absolute deviation of
#' the posterior distribution, rather
#'   than mean and standard deviation, to derive point estimates and uncertainty
#' @param conf.int If \code{TRUE} columns for the lower (\code{conf.low})
#' and upper bounds (\code{conf.high}) of posterior uncertainty intervals are included.
#' @param exponentiate  whether to exponentiate the fixed-effect coefficient estimates and confidence intervals (common for logistic regression); if \code{TRUE}, also scales the standard errors by the exponentiated coefficient, transforming them to the new scale
#' @param conf.level Defines the range of the posterior uncertainty conf.int,
#'  such that \code{100 * conf.level}\% of the parameter's posterior distributio
#'  lies within the corresponding interval.
#'  Only used if \code{conf.int = TRUE}.
#' @param conf.method method for computing confidence intervals
#' ("quantile" or "HPDinterval")
#' @param rhat whether to calculate the *Rhat* convergence metric
#' (\code{FALSE} by default)
#' @param ess whether to calculate the *effective sample size* (ESS) convergence metric
#' (\code{FALSE} by default)
#' @param fix.intercept rename "Intercept" parameter to "(Intercept)", to match
#' behaviour of other model types?
#' @param looic Should the LOO Information Criterion (and related info) be
#'   included? See \code{\link[rstan]{loo.stanfit}} for details. (This
#'   can be slow for models fit to large datasets.)
#' @param ... Extra arguments, not used
#' @return
#' When \code{parameters = NA}, the \code{effects} argument is used
#' to determine which parameters to summarize.
#'
#' Generally, \code{tidy.brmsfit} returns
#' one row for each coefficient, with at least three columns:
#' \item{term}{The name of the model parameter.}
#' \item{estimate}{A point estimate of the coefficient (mean or median).}
#' \item{std.error}{A standard error for the point estimate (sd or mad).}
#'
#' When \code{effects = "fixed"}, only population-level
#' effects are returned.
#'
#' When \code{effects = "ran_vals"}, only group-level effects are returned.
#' In this case, two additional columns are added:
#' \item{group}{The name of the grouping factor.}
#' \item{level}{The name of the level of the grouping factor.}
#'
#' Specifying \code{effects = "ran_pars"} selects the
#' standard deviations and correlations of the group-level parameters.
#'
#' If \code{conf.int = TRUE}, columns for the \code{lower} and
#' \code{upper} bounds of the posterior conf.int computed.
#'
#' @note The names \sQuote{fixed}, \sQuote{ran_pars}, and \sQuote{ran_vals}
#' (corresponding to "non-varying", "hierarchical", and "varying" respectively
#' in previous versions of the package), while technically inappropriate in
#' a Bayesian setting where "fixed" and "random" effects are not well-defined,
#' are used for compatibility with other (frequentist) mixed model types.
#' @note This implementation uses brms built-in functions (fixef, VarCorr, ranef)
#' for extracting posterior summaries when possible, falling back to manual
#' posterior sample processing for HPD intervals and convergence diagnostics.
#' @export
tidy.brmsfit <- function(x, parameters = NA,
                         effects = c("fixed", "ran_pars"),
                         robust = FALSE,
                         conf.int = TRUE, conf.level = 0.95,
                         conf.method = c("quantile", "HPDinterval"),
                         rhat = FALSE, ess = FALSE,
                         fix.intercept = TRUE,
                         exponentiate = FALSE,
                         ...) {

  check_dots(...)
  bad_effects <- setdiff(effects, c("fixed", "ran_pars", "ran_vals", "ran_coefs"))
  if (length(bad_effects) > 0) {
    stop("unrecognized effects: ", paste(bad_effects, collapse = ", "))
  }
  std.error <- NULL ## NSE/code check
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("can't tidy brms objects without brms installed")
  }

  conf.method <- match.arg(conf.method)
  use_effects <- anyNA(parameters)
  is.multiresp <- length(x$formula$forms) > 1

  ## Check for underscores in parameter names
  xr <- brms::restructure(x)
  has_ranef <- nrow(xr$ranef) > 0
  if (any(grepl("_", rownames(brms::fixef(x)))) ||
      (has_ranef && any(grepl("_", names(brms::ranef(x)))))) {
    warning("some parameter names contain underscores: term naming may be unreliable!")
  }

  ## Set up probability quantiles for confidence intervals
  probs <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
  sep <- getOption("broom.mixed.sep1")

  if (!use_effects) {
    ## Custom parameters specified - use sample-based approach
    out <- tidy_brms_parameters(x, parameters, robust, conf.int, conf.level,
                                 conf.method, rhat, ess)
  } else {
    ## Use brms built-in functions for standard effects
    res_list <- list()

    if ("fixed" %in% effects) {
      res_list$fixed <- tidy_brms_fixef(x, robust, conf.int, probs,
                                         conf.method, conf.level, sep)
    }

    if ("ran_pars" %in% effects) {
      ran_pars_result <- tidy_brms_varcorr(x, robust, conf.int, probs,
                                            conf.method, conf.level, sep)
      if (!is.null(ran_pars_result) && nrow(ran_pars_result) > 0) {
        res_list$ran_pars <- ran_pars_result
      }
      if (is.multiresp) {
        warning("ran_pars response/group tidying for multi-response models is currently incorrect")
      }
    }

    if ("ran_vals" %in% effects) {
      ran_vals_result <- tidy_brms_ranef(x, robust, conf.int, probs,
                                          conf.method, conf.level, sep)
      if (!is.null(ran_vals_result) && nrow(ran_vals_result) > 0) {
        res_list$ran_vals <- ran_vals_result
      }
    }

    out <- dplyr::bind_rows(res_list, .id = "effect")

    ## Handle empty data frame case
    if (nrow(out) == 0) {
      stop("No parameter name matches the specified pattern.", call. = FALSE)
    }

    ## Add rhat/ess if requested (requires raw samples)
    if (rhat || ess) {
      out <- add_posterior_metrics(x, out, effects, rhat, ess)
    }
  }

  ## Add component column based on term names
  out$component <- dplyr::case_when(
    grepl("(^|_)zi", out$term) ~ "zi",
    grepl("^disp", out$term) ~ "disp",
    TRUE ~ "cond"
  )

  ## Remove component prefixes from terms
  comp_pref_RE <- "(?<=((^|_)))(zi_|disp_)"
  out$term <- stringr::str_remove(out$term, comp_pref_RE)

  ## Extract response column for multi-response models
  if (is.multiresp) {
    resp_names <- names(x$formula$forms)
    ## Pattern: response_term -> extract response prefix
    resp_pattern <- paste0("^(", paste(resp_names, collapse = "|"), ")_")
    out$response <- stringr::str_extract(out$term, resp_pattern)
    out$response <- stringr::str_remove(out$response, "_$")
    out$term <- stringr::str_remove(out$term, resp_pattern)
  }

  if (exponentiate) {
    vv <- c("estimate", "conf.low", "conf.high")
    out <- (out
      %>% mutate(across(contains(vv), exp))
      %>% mutate(across(std.error, ~ . * estimate))
    )
  }

  if (fix.intercept) {
    out$term <- stringr::str_replace(out$term,
                                      "(?<=(\\b|_))Intercept(?=(\\b|_))",
                                      "(Intercept)")
  }

  out <- reorder_cols(out)
  return(out)
}

#' Tidy fixed effects using brms::fixef()
#' @noRd
tidy_brms_fixef <- function(x, robust, conf.int, probs, conf.method, conf.level, sep) {
  ## Use brms::fixef() for summary statistics
  fe <- brms::fixef(x, summary = TRUE, robust = robust, probs = probs)

  if (is.null(fe) || nrow(fe) == 0) {
    return(dplyr::tibble(
      group = character(),
      term = character(),
      estimate = numeric(),
      std.error = numeric()
    ))
  }

  out <- dplyr::tibble(
    group = NA_character_,
    term = rownames(fe),
    estimate = fe[, "Estimate"],
    std.error = fe[, "Est.Error"]
  )

  if (conf.int) {
    if (conf.method == "HPDinterval") {
      ## Need raw samples for HPD intervals
      fe_names <- rownames(fe)
      samples <- brms::fixef(x, summary = FALSE)
      cc <- coda::HPDinterval(coda::as.mcmc(samples), prob = conf.level)
      out$conf.low <- cc[, 1]
      out$conf.high <- cc[, 2]
    } else {
      ## Use quantiles from fixef() output
      col_names <- colnames(fe)
      q_cols <- grep("^Q|^[0-9]", col_names, value = TRUE)
      if (length(q_cols) >= 2) {
        out$conf.low <- fe[, q_cols[1]]
        out$conf.high <- fe[, q_cols[2]]
      }
    }
  }

  return(out)
}

#' Tidy variance-covariance components using brms::VarCorr()
#' @noRd
tidy_brms_varcorr <- function(x, robust, conf.int, probs, conf.method, conf.level, sep) {
  ## Use brms::VarCorr() for variance components
  ## VarCorr throws an error for models without random effects
  vc <- tryCatch(
    brms::VarCorr(x, summary = TRUE, robust = robust, probs = probs),
    error = function(e) NULL
  )

  if (is.null(vc) || length(vc) == 0) {
    ## Model may still have sigma (residual SD) without random effects
    vc <- list()
  }

  results <- list()

  for (group_name in names(vc)) {
    ## Skip residual__ group - we handle sigma separately below
    if (group_name == "residual__") {
      next
    }

    group_vc <- vc[[group_name]]

    ## Extract standard deviations
    if (!is.null(group_vc$sd)) {
      sd_mat <- group_vc$sd
      if (is.matrix(sd_mat) || is.array(sd_mat)) {
        if (length(dim(sd_mat)) == 2) {
          for (i in seq_len(nrow(sd_mat))) {
            term_name <- rownames(sd_mat)[i]
            row_data <- dplyr::tibble(
              group = group_name,
              term = paste0("sd", sep, term_name),
              estimate = sd_mat[i, "Estimate"],
              std.error = sd_mat[i, "Est.Error"]
            )
            if (conf.int && conf.method != "HPDinterval" && ncol(sd_mat) >= 4) {
              row_data$conf.low <- sd_mat[i, 3]
              row_data$conf.high <- sd_mat[i, 4]
            }
            results[[length(results) + 1]] <- row_data
          }
        }
      }
    }

    ## Extract correlations
    ## Note: VarCorr cor array has dimensions [vars, stats, vars]
    if (!is.null(group_vc$cor)) {
      cor_arr <- group_vc$cor
      if (is.array(cor_arr) && length(dim(cor_arr)) == 3) {
        var_names <- dimnames(cor_arr)[[1]]
        n_vars <- length(var_names)
        if (n_vars > 1) {
          for (i in 2:n_vars) {
            for (j in 1:(i - 1)) {
              term_name <- paste0("cor", sep, var_names[j], ".", var_names[i])
              row_data <- dplyr::tibble(
                group = group_name,
                term = term_name,
                estimate = cor_arr[i, "Estimate", j],
                std.error = cor_arr[i, "Est.Error", j]
              )
              if (conf.int && conf.method != "HPDinterval" && dim(cor_arr)[2] >= 4) {
                row_data$conf.low <- cor_arr[i, 3, j]
                row_data$conf.high <- cor_arr[i, 4, j]
              }
              results[[length(results) + 1]] <- row_data
            }
          }
        }
      }
    }
  }

  ## Check for residual sigma (single or response-specific)
  all_vars <- brms::variables(x)
  sigma_vars <- grep("^sigma($|_)", all_vars, value = TRUE)
  for (sigma_var in sigma_vars) {
    sigma_summary <- brms::posterior_summary(x, variable = sigma_var,
                                              robust = robust, probs = probs)
    ## For multi-response models, sigma_response -> response-specific term
    if (sigma_var == "sigma") {
      term_name <- paste0("sd", sep, "Observation")
    } else {
      ## sigma_response -> sd__response_Observation
      resp_name <- sub("^sigma_", "", sigma_var)
      term_name <- paste0(resp_name, "_sd", sep, "Observation")
    }
    row_data <- dplyr::tibble(
      group = "Residual",
      term = term_name,
      estimate = sigma_summary[, "Estimate"],
      std.error = sigma_summary[, "Est.Error"]
    )
    if (conf.int && conf.method != "HPDinterval") {
      row_data$conf.low <- sigma_summary[, 3]
      row_data$conf.high <- sigma_summary[, 4]
    }
    results[[length(results) + 1]] <- row_data
  }

  if (length(results) == 0) {
    return(NULL)
  }

  out <- dplyr::bind_rows(results)

  ## Handle HPD intervals if requested
  if (conf.int && conf.method == "HPDinterval") {
    out <- compute_varcorr_hpd(x, out, conf.level, sep)
  }

  return(out)
}

#' Compute HPD intervals for variance components
#' @noRd
compute_varcorr_hpd <- function(x, out, conf.level, sep) {
  ## Get all sd_ and cor_ parameters, plus sigma (single or response-specific)
  all_vars <- brms::variables(x)
  sd_vars <- grep("^sd_", all_vars, value = TRUE)
  cor_vars <- grep("^cor_", all_vars, value = TRUE)
  sigma_vars <- grep("^sigma($|_)", all_vars, value = TRUE)

  vars_to_get <- c(sd_vars, cor_vars, sigma_vars)

  if (length(vars_to_get) > 0) {
    samples <- brms::as_draws_matrix(x, variable = vars_to_get)
    cc <- coda::HPDinterval(coda::as.mcmc(samples), prob = conf.level)

    out$conf.low <- NA_real_
    out$conf.high <- NA_real_

    ## Match parameter names to output rows
    for (i in seq_len(nrow(out))) {
      ## Build possible parameter name patterns
      if (out$term[i] == paste0("sd", sep, "Observation")) {
        param_name <- "sigma"
      } else if (grepl(paste0("_sd", sep, "Observation$"), out$term[i])) {
        ## Response-specific sigma: resp_sd__Observation -> sigma_resp
        resp_name <- sub(paste0("_sd", sep, "Observation$"), "", out$term[i])
        param_name <- paste0("sigma_", resp_name)
      } else if (grepl(paste0("^sd", sep), out$term[i])) {
        ## sd parameter: sd__term -> sd_group__term
        term_part <- sub(paste0("^sd", sep), "", out$term[i])
        param_name <- paste0("sd_", out$group[i], "__", term_part)
      } else if (grepl(paste0("^cor", sep), out$term[i])) {
        ## cor parameter: cor__t1.t2 -> cor_group__t1__t2
        term_part <- sub(paste0("^cor", sep), "", out$term[i])
        terms <- strsplit(term_part, "\\.")[[1]]
        param_name <- paste0("cor_", out$group[i], "__", terms[1], "__", terms[2])
      } else {
        next
      }

      if (param_name %in% rownames(cc)) {
        out$conf.low[i] <- cc[param_name, 1]
        out$conf.high[i] <- cc[param_name, 2]
      }
    }
  }

  return(out)
}

#' Tidy random effects values using brms::ranef()
#' @noRd
tidy_brms_ranef <- function(x, robust, conf.int, probs, conf.method, conf.level, sep) {
  ## Use brms::ranef() for random effect values
  re <- brms::ranef(x, summary = TRUE, robust = robust, probs = probs)

  if (is.null(re) || length(re) == 0) {
    return(NULL)
  }

  results <- list()

  for (group_name in names(re)) {
    group_re <- re[[group_name]]
    ## group_re is a 3D array: [levels, statistics, terms]
    if (is.array(group_re) && length(dim(group_re)) == 3) {
      levels <- dimnames(group_re)[[1]]
      terms <- dimnames(group_re)[[3]]

      for (term in terms) {
        for (level in levels) {
          row_data <- dplyr::tibble(
            group = group_name,
            level = level,
            term = term,
            estimate = group_re[level, "Estimate", term],
            std.error = group_re[level, "Est.Error", term]
          )
          if (conf.int && conf.method != "HPDinterval" && dim(group_re)[2] >= 4) {
            row_data$conf.low <- group_re[level, 3, term]
            row_data$conf.high <- group_re[level, 4, term]
          }
          results[[length(results) + 1]] <- row_data
        }
      }
    }
  }

  if (length(results) == 0) {
    return(NULL)
  }

  out <- dplyr::bind_rows(results)

  ## Handle HPD intervals if requested
  if (conf.int && conf.method == "HPDinterval") {
    out <- compute_ranef_hpd(x, out, conf.level)
  }

  return(out)
}

#' Compute HPD intervals for random effects
#' @noRd
compute_ranef_hpd <- function(x, out, conf.level) {
  ## Get random effect samples
  all_vars <- brms::variables(x)
  r_vars <- grep("^r_", all_vars, value = TRUE)

  if (length(r_vars) > 0) {
    samples <- brms::as_draws_matrix(x, variable = r_vars)
    cc <- coda::HPDinterval(coda::as.mcmc(samples), prob = conf.level)

    out$conf.low <- NA_real_
    out$conf.high <- NA_real_

    ## Match parameter names to output rows
    for (i in seq_len(nrow(out))) {
      ## Construct expected parameter name: r_group[level,term]
      param_name <- paste0("r_", out$group[i], "[", out$level[i], ",", out$term[i], "]")
      if (param_name %in% rownames(cc)) {
        out$conf.low[i] <- cc[param_name, 1]
        out$conf.high[i] <- cc[param_name, 2]
      }
    }
  }

  return(out)
}

#' Add posterior metrics (rhat, ess) to tidy output
#' @noRd
add_posterior_metrics <- function(x, out, effects, rhat, ess) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(paste0(c(if (rhat) "rhat", if (ess) "ess"), collapse = ", "),
         " calculation for brmsfit objects requires posterior package")
  }

  ## Build parameter pattern to match effects
  prefs <- list(
    fixed = "b_", ran_vals = "r_",
    ran_pars = c("sd_", "cor_", "sigma")
  )
  mkRE <- function(x) {
    sprintf("(^|_)(%s)", paste(unlist(x), collapse = "|"))
  }
  pref_RE <- mkRE(prefs[effects])

  samples_perchain <- brms::as_draws_array(x, pref_RE, regex = TRUE)

  posterior_metrics <- c()
  if (rhat) {
    posterior_metrics <- c(posterior_metrics, rhat = posterior::rhat)
  }
  if (ess) {
    posterior_metrics <- c(posterior_metrics, ess = posterior::ess_basic)
  }

  metrics_df <- posterior::summarise_draws(samples_perchain, posterior_metrics)

  ## Initialize metric columns
  out[names(posterior_metrics)] <- NA_real_

  sep <- getOption("broom.mixed.sep1")

  ## Match metrics to output rows
  for (i in seq_len(nrow(out))) {
    param_name <- NULL

    if (out$effect[i] == "fixed") {
      ## Fixed effects: b_term
      term_name <- gsub("\\(Intercept\\)", "Intercept", out$term[i])
      param_name <- paste0("b_", term_name)
    } else if (out$effect[i] == "ran_pars") {
      ## Variance parameters
      if (out$term[i] == paste0("sd", sep, "Observation")) {
        param_name <- "sigma"
      } else if (grepl(paste0("^sd", sep), out$term[i])) {
        term_part <- sub(paste0("^sd", sep), "", out$term[i])
        param_name <- paste0("sd_", out$group[i], "__", term_part)
      } else if (grepl(paste0("^cor", sep), out$term[i])) {
        term_part <- sub(paste0("^cor", sep), "", out$term[i])
        terms <- strsplit(term_part, "\\.")[[1]]
        param_name <- paste0("cor_", out$group[i], "__", terms[1], "__", terms[2])
      }
    } else if (out$effect[i] == "ran_vals") {
      ## Random values: r_group[level,term]
      param_name <- paste0("r_", out$group[i], "[", out$level[i], ",", out$term[i], "]")
    }

    if (!is.null(param_name)) {
      match_idx <- which(metrics_df$variable == param_name)
      if (length(match_idx) == 1) {
        for (metric in names(posterior_metrics)) {
          out[[metric]][i] <- metrics_df[[metric]][match_idx]
        }
      }
    }
  }

  return(out)
}

#' Process custom parameter patterns (legacy approach)
#' @noRd
tidy_brms_parameters <- function(x, parameters, robust, conf.int, conf.level,
                                  conf.method, rhat, ess) {
  samples_perchain <- brms::as_draws_array(x, parameters, regex = TRUE)
  if (is.null(samples_perchain) || posterior::nvariables(samples_perchain) == 0) {
    stop("No parameter name matches the specified pattern.", call. = FALSE)
  }
  samples <- brms::as_draws_matrix(samples_perchain)
  terms <- colnames(samples)

  out <- dplyr::tibble(term = terms)

  pointfun <- if (robust) stats::median else base::mean
  stdfun <- if (robust) stats::mad else stats::sd
  out$estimate <- apply(samples, 2, pointfun)
  out$std.error <- apply(samples, 2, stdfun)

  if (conf.int) {
    probs <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
    if (conf.method == "HPDinterval") {
      cc <- coda::HPDinterval(coda::as.mcmc(samples), prob = conf.level)
    } else {
      cc <- t(apply(samples, 2, stats::quantile, probs = probs))
    }
    out$conf.low <- cc[, 1]
    out$conf.high <- cc[, 2]
  }

  if (rhat || ess) {
    if (!requireNamespace("posterior", quietly = TRUE)) {
      stop(paste0(c(if (rhat) "rhat", if (ess) "ess"), collapse = ", "),
           " calculation for brmsfit objects requires posterior package")
    }
    posterior_metrics <- c()
    if (rhat) posterior_metrics <- c(posterior_metrics, rhat = posterior::rhat)
    if (ess) posterior_metrics <- c(posterior_metrics, ess = posterior::ess_basic)
    out[names(posterior_metrics)] <- posterior::summarise_draws(samples_perchain, posterior_metrics)[names(posterior_metrics)]
  }

  return(out)
}


#' @importFrom stats quantile
#' @export
sigma.brmsfit <- function (object, ...)  {
    if (!("sigma" %in% brms::variables(object)))
        return(1)
    brms::posterior_summary(object, variable = "sigma", probs = 0.5)[, "Estimate"]
}

#' @rdname brms_tidiers
#' @export
glance.brmsfit <- function(x, looic = FALSE, ...) {
  ## defined in rstanarm_tidiers.R
  glance_stan(x, looic = looic, type = "brmsfit", ...)
}

#' @rdname brms_tidiers
#' @param data data frame
#' @param newdata new data frame
#' @param se.fit return standard errors of fit?
#' @export
augment.brmsfit <- function(x, data = stats::model.frame(x), newdata = NULL,
                            se.fit = TRUE, ...) {
  ## Use brms::predict() which handles everything internally
  args <- list(x, se.fit = se.fit)
  if (!missing(newdata)) args$newdata <- newdata

  pred <- do.call(stats::predict, args)
  ret <- dplyr::tibble(.fitted = pred[, "Estimate"])
  if (se.fit) ret[[".se.fit"]] <- pred[, "Est.Error"]
  if (is.null(newdata)) {
    ret[[".resid"]] <- stats::residuals(x)[, "Estimate"]
    ret <- dplyr::bind_cols(as_tibble(data), ret)
  } else {
    ret <- dplyr::bind_cols(as_tibble(newdata), ret)
  }
  return(ret)
}
