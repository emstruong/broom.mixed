#' Tidying methods for an rstanarm model
#'
#' These methods tidy the estimates from \code{rstanarm} fits
#' (\code{stan_glm}, \code{stan_glmer}, etc.)
#' into a summary.
#'
#' @return All tidying methods return a \code{data.frame} without rownames.
#' The structure depends on the method chosen.
#'
#' @seealso \code{\link[rstan]{summary,stanfit-method}}
#'
#' @name rstanarm_tidiers
#'
#' @param x Fitted model object from the \pkg{rstanarm} package. See
#'   \code{\link[rstanarm]{stanreg-objects}}.
#' @examples
#'
#' if (require("rstanarm") && require("tibble")) {
#' \dontrun{
#' #'     ## original models
#'   fit <- stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
#'                       iter = 500, chains = 2)
#'   fit2 <- stan_glmer((mpg>20) ~ wt + (1 | cyl) + (1 + wt | gear),
#'                     data = mtcars,
#'                     family = binomial,
#'                     iter = 500, chains = 2
#'   }
#' ## load example data
#'   load(system.file("extdata", "rstanarm_example.rda", package="broom.mixed"))
#'
#'   # non-varying ("population") parameters
#'   tidy(fit, conf.int = TRUE, conf.level = 0.5)
#'   tidy(fit, conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.5)
#'
#'   #  exponentiating (in this case, from log-odds to odds ratios)
#'   (tidy(fit2, conf.int = TRUE, conf.level = 0.5)
#'           %>% dplyr::filter(term != "(Intercept)")
#'   )
#'   (tidy(fit2, conf.int = TRUE, conf.level = 0.5, exponentiate = TRUE)
#'           %>% dplyr::filter(term != "(Intercept)")
#'   )
#' 
#'   # hierarchical sd & correlation parameters
#'   tidy(fit, effects = "ran_pars")
#'
#'   # group-specific deviations from "population" parameters
#'   tidy(fit, effects = "ran_vals")
#'
#'   # glance method
#'    glance(fit)
#'   \dontrun{
#'      glance(fit, looic = TRUE, cores = 1)
#'   }
#' } ## if require("rstanarm")
NULL

#' @rdname rstanarm_tidiers
#' @inheritParams brms_tidiers
#' @param robust Whether to use median and median absolute deviation of
#'   the posterior distribution, rather than mean and standard deviation,
#'   to derive point estimates and uncertainty. Defaults to \code{TRUE}
#'   for consistency with \pkg{rstanarm}'s default behavior.
#' @param conf.level See \code{\link[rstantools]{posterior_interval}}.
#' @param conf.int If \code{TRUE} columns for the lower (\code{conf.low}) and upper (\code{conf.high}) bounds of the
#'   \code{100*prob}\% posterior uncertainty intervals are included. See
#'   \code{\link[rstantools]{posterior_interval}} for details.
#' @param rhat Whether to include the *Rhat* convergence metric
#'   (\code{FALSE} by default). Uses the Rhat values computed by rstanarm
#'   and stored in the model summary.
#' @param ess Whether to include the *effective sample size* (ESS) convergence metric
#'   (\code{FALSE} by default). Uses the n_eff values computed by rstanarm
#'   and stored in the model summary.
#' @param fix.intercept Rename \code{"Intercept"} to \code{"(Intercept)"},
#'   to match behavior of other model types? Defaults to \code{TRUE}.
#'
#' @return
#' When \code{effects="fixed"} (the default), \code{tidy.stanreg} returns
#' one row for each coefficient, with three columns:
#' \item{term}{The name of the corresponding term in the model.}
#' \item{estimate}{A point estimate of the coefficient (posterior median if
#'   \code{robust=TRUE}, posterior mean if \code{robust=FALSE}).}
#' \item{std.error}{A standard error for the point estimate based on
#' \code{\link[stats]{mad}} if \code{robust=TRUE}, or standard deviation
#' if \code{robust=FALSE}. See the \emph{Uncertainty estimates} section in
#' \code{\link[rstanarm]{print.stanreg}} for more details.}
#'
#' For models with group-specific parameters (e.g., models fit with
#' \code{\link[rstanarm]{stan_glmer}}), setting \code{effects="ran_vals"}
#' selects the group-level parameters instead of the non-varying regression
#' coefficients. Addtional columns are added indicating the \code{level} and
#' \code{group}. Specifying \code{effects="ran_pars"} selects the
#' standard deviations and (for certain models) correlations of the group-level
#' parameters.
#'
#' Setting \code{effects="auxiliary"} will select parameters other than those
#' included by the other options. The particular parameters depend on which
#' \pkg{rstanarm} modeling function was used to fit the model. For example, for
#' models fit using \code{\link[rstanarm]{stan_glm}} the overdispersion
#' parameter is included if \code{effects="aux"}, for
#' \code{\link[rstanarm]{stan_lm}} the auxiliary parameters include the residual
#' SD, R^2, and log(fit_ratio), etc.
#'
#' @export
tidy.stanreg <- function(x,
                         effects = c("fixed", "ran_pars"),
                         robust = TRUE,
                         conf.int = FALSE,
                         conf.level = 0.9,
                         conf.method = c("quantile", "HPDinterval"),
                         rhat = FALSE,
                         ess = FALSE,
                         fix.intercept = TRUE,
                         exponentiate = FALSE,
                         ...) {
    ## ignore 'parametric', which may be passed by mice:::summary.mira()
    check_dots(..., .ignore = "parametric")
    conf.method <- match.arg(conf.method)
    std.error <- estimate <- NULL ## fool code checker/NSE
    miss_effects <- missing(effects)
    effects <-
        match.arg(effects,
                  several.ok = TRUE,
                  choices = c(
                      "fixed", "ran_vals",
                      "ran_pars", "auxiliary"
                  )
                  )
    no_ranef <- !inherits(x, "lmerMod")
    if (miss_effects && no_ranef) effects <- setdiff(effects, c("ran_vals", "ran_pars"))
    if (!miss_effects && no_ranef && any(effects %in% c("ran_vals", "ran_pars"))) {
        stop("Model does not have varying ('ran_vals') or hierarchical ('ran_pars') effects.")
    }

    nn <- c("estimate", "std.error")
    ret_list <- list()

    ## Track parameter names for rhat/ess lookup
    pars_for_diag <- character(0)

    ## Get the summary matrix from rstanarm (contains mean, sd, 50%, Rhat, n_eff, etc.)
    stan_summary <- x$stan_summary

    if ("fixed" %in% effects) {
        nv_pars <- names(rstanarm::fixef(x))

        if (robust) {
            ## Use rstanarm's built-in functions (median/mad)
            ret <- cbind(
                rstanarm::fixef(x),
                rstanarm::se(x)[nv_pars]
            )
        } else {
            ## Use rstanarm's summary which already contains mean and sd
            ret <- stan_summary[nv_pars, c("mean", "sd"), drop = FALSE]
        }

        if (inherits(x, "polr")) {
            ## also include cutpoints
            cp_names <- names(x$zeta)
            if (robust) {
                cp <- x$zeta
                se_cp <- rstanarm::se(x)[cp_names]
            } else {
                cp <- stan_summary[cp_names, "mean"]
                se_cp <- stan_summary[cp_names, "sd"]
            }
            ret <- rbind(ret, cbind(cp, se_cp))
            nv_pars <- c(nv_pars, cp_names)
        }

        if (conf.int) {
            cifix <- switch(conf.method,
                            HPDinterval = {
                                m <- as.matrix(x$stanfit)
                                m <- m[, colnames(m) %in% nv_pars, drop = FALSE]
                                coda::HPDinterval(coda::as.mcmc(m),
                                                  prob = conf.level)
                            },
                            quantile = rstanarm::posterior_interval(
                                object = x,
                                pars = nv_pars,
                                prob = conf.level
                            )
            ) ## cifix
            ret <- data.frame(ret, cifix)
            nn <- c(nn, "conf.low", "conf.high")
        }
        ret_list$non_ran_vals <- fix_data_frame(ret, newnames = nn, newcol = "term")
        pars_for_diag <- c(pars_for_diag, nv_pars)
    }
    if ("auxiliary" %in% effects) {
        nn <- c("estimate", "std.error")
        parnames <- rownames(stan_summary)
        auxpars <- c(
            "sigma", "shape", "overdispersion", "R2", "log-fit_ratio",
            grep("mean_PPD", parnames, value = TRUE)
        )
        auxpars <- auxpars[which(auxpars %in% parnames)]

        if (robust) {
            ## Use stan_summary's 50% (median) for estimate
            ## Note: rstanarm's se() uses MAD for robust SE, but stan_summary
            ## only has sd. For auxiliary params, use sd as approximation.
            ret <- stan_summary[auxpars, c("50%", "sd"), drop = FALSE]
        } else {
            ## Use stan_summary's mean and sd
            ret <- stan_summary[auxpars, c("mean", "sd"), drop = FALSE]
        }

        if (conf.int) {
            ints <- rstanarm::posterior_interval(x, pars = auxpars, prob = conf.level)
            ret <- data.frame(ret, ints)
            nn <- c(nn, "conf.low", "conf.high")
        }
        ret_list$auxiliary <-
            fix_data_frame(ret, newnames = nn, newcol = "term")
        pars_for_diag <- c(pars_for_diag, auxpars)
    }
    if ("ran_pars" %in% effects) {
        ret <- (rstanarm::VarCorr(x)
            %>% as.data.frame()
            %>% mutate_if(is.factor,as.character)
        )
        rscale <- "sdcor" # FIXME
        ran_prefix <- c("sd", "cor") # FIXME
        pfun <- function(x) {
            v <- na.omit(unlist(x))
            if (length(v) == 0) v <- "Observation"
            p <- paste(v, collapse = ".")
            if (!identical(ran_prefix, NA)) {
                p <- paste(ran_prefix[length(v)], p, sep = "_")
            }
            return(p)
        }

        rownames(ret) <- paste(apply(ret[c("var1", "var2")], 1, pfun),
                               ret[, "grp"],
                               sep = "."
                               )
        ret_list$hierarchical <- fix_data_frame(ret[c("grp", rscale)],
                                                 newcol="term",
                                                newnames = c("group", "estimate"))
    }

    if ("ran_vals" %in% effects) {
        nn <- c("estimate", "std.error")
        ## Get random effect parameter names (those starting with "b[")
        ran_val_pars <- grep("^b\\[", rownames(stan_summary), value = TRUE)

        if (robust) {
            ## Use stan_summary's 50% (median) and MAD from rstanarm::se()
            ret <- cbind(
                stan_summary[ran_val_pars, "50%"],
                rstanarm::se(x)[ran_val_pars]
            )
        } else {
            ## Use stan_summary's mean and sd
            ret <- stan_summary[ran_val_pars, c("mean", "sd"), drop = FALSE]
        }

        if (conf.int) {
            ciran <- rstanarm::posterior_interval(x,
                                                  regex_pars = "^b\\[",
                                                  prob = conf.level
            )
            ret <- data.frame(ret, ciran)
            nn <- c(nn, "conf.low", "conf.high")
        }

        double_splitter <- function(x, split1, sel1, split2, sel2) {
            y <- unlist(lapply(strsplit(x, split = split1, fixed = TRUE), "[[", sel1))
            unlist(lapply(strsplit(y, split = split2, fixed = TRUE), "[[", sel2))
        }
        vv <- fix_data_frame(ret, newnames = nn, newcol = "term")
        nn <- c("level", "group", "term", nn)
        nms <- vv$term
        vv$term <- NULL
        lev <- double_splitter(nms, ":", 2, "]", 1)
        grp <- double_splitter(nms, " ", 2, ":", 1)
        trm <- double_splitter(nms, " ", 1, "[", 2)
        vv <- data.frame(lev, grp, trm, vv)
        ret_list$ran_vals <- fix_data_frame(vv, newnames = nn, newcol = "term")
        pars_for_diag <- c(pars_for_diag, ran_val_pars)
    }

    if (exponentiate) {
        ret_list$non_ran_vals <- (ret_list$non_ran_vals
            %>% mutate(across(any_of(c("estimate", "conf.low", "conf.high")), exp))
            %>% mutate(std.error = std.error * estimate)
        )
    }

    out <- dplyr::bind_rows(ret_list)

    ## Add rhat and ess if requested
    ## Use rstanarm's built-in Rhat and n_eff from stan_summary
    if (rhat || ess) {
        if (length(pars_for_diag) > 0) {
            available_pars <- intersect(pars_for_diag, rownames(stan_summary))
            if (length(available_pars) > 0) {
                if (rhat && "Rhat" %in% colnames(stan_summary)) {
                    out$rhat <- NA_real_
                    for (par_name in available_pars) {
                        match_idx <- which(out$term == par_name)
                        if (length(match_idx) > 0) {
                            out$rhat[match_idx] <- stan_summary[par_name, "Rhat"]
                        }
                    }
                }
                if (ess && "n_eff" %in% colnames(stan_summary)) {
                    out$ess <- NA_real_
                    for (par_name in available_pars) {
                        match_idx <- which(out$term == par_name)
                        if (length(match_idx) > 0) {
                            out$ess[match_idx] <- stan_summary[par_name, "n_eff"]
                        }
                    }
                }
            }
        }
    }

    ## Apply fix.intercept: rstanarm already uses "(Intercept)" by default.
    ## When fix.intercept = FALSE, convert back to "Intercept" (no parentheses)
    ## to match brms behavior when fix.intercept = FALSE.
    if (!fix.intercept && "term" %in% names(out)) {
        out$term <- gsub("^\\(Intercept\\)$", "Intercept", out$term)
    }

    return(reorder_cols(out))
}


#' @rdname rstanarm_tidiers
#'
#' @param ... For \code{glance}, if \code{looic=TRUE}, optional arguments to
#'   \code{\link[rstan]{loo.stanfit}}.
#' @return \code{glance} returns one row with the columns
#'   \item{algorithm}{The algorithm used to fit the model.}
#'   \item{pss}{The posterior sample size (except for models fit using
#'   optimization).}
#'   \item{nobs}{The number of observations used to fit the model.}
#'   \item{sigma}{The square root of the estimated residual variance, if
#'   applicable. If not applicable (e.g., for binomial GLMs), \code{sigma} will
#'   be given the value \code{1} in the returned object.}
#'
#'   If \code{looic=TRUE}, then the following additional columns are also
#'   included:
#'   \item{looic}{The LOO Information Criterion.}
#'   \item{elpd_loo}{The expected log predictive density (\code{elpd_loo = -2 *
#'   looic}).}
#'   \item{p_loo}{The effective number of parameters.}
#'
#' @export
glance.stanreg <- function(x, looic = FALSE, ...) {
    glance_stan(x, looic = looic, type = "stanreg", ...)
}


glance_stan <- function(x, looic = FALSE, ..., type) {
    sigma <- if (getRversion() >= "3.3.0") {
                 get("sigma", asNamespace("stats"))
             } else {
                 ## FIXME: could fail if old R & called from brms
                 ## & rstanarm not installed ...
                 get("sigma", asNamespace("rstanarm"))
             }
    if (type == "stanreg") {
        algo <- x$algorithm
        sim <- x$stanfit@sim
    } else {
        ## method is recorded for every chain; pick the first
        algo <- x$fit@stan_args[[1]][["method"]]
        sim <- x$fit@sim
    }

    ret <- dplyr::tibble(algorithm = algo)

    if (algo != "optimizing") {
        pss <- sim$n_save
        if (algo %in% c("sample", "sampling")) {
            pss <- pss - sim$warmup2
        }
        ret <- dplyr::mutate(ret, pss = sum(pss))
    }

    ret <- mutate(ret, nobs = stats::nobs(x))
    if (length(sx <- sigma(x)) > 0) {
        ret <- dplyr::mutate(ret, sigma = sx)
    }
    if (looic) {
        if (algo == "sampling") {
            if (type == "stanreg") {
                loo1 <- rstanarm::loo(x, ...)
            } else {
                loo1 <- brms::loo(x, ...)
            }
            loo1_est <- loo1[["estimates"]]
            ret <- data.frame(
                ret,
                rbind(loo1_est[
                    c("looic", "elpd_loo", "p_loo"),
                    "Estimate"
                ])
            )
        } else {
            message("looic only available for models fit using MCMC")
        }
    }
    dplyr::as_tibble(unrowname(ret))
}
