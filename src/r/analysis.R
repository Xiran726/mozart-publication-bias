# analysis.R --------------------------------------------------------------

# Statistical analyses and visualisations for the Mozart-effect
# publication bias simulation.
#
# This script compares meta-analytic estimates from published studies,
# all simulated studies, and the known true effect.
# It also evaluates Egger's test and trim-and-fill across different
# simulation conditions.
#
# Required packages: metafor, ggplot2

# 1. Random-effects meta-analysis -----------------------------------------

#' Fit a random-effects meta-analysis
#'
#' Combines the Hedges' g values from individual studies to estimate 
#' the overall effect. A random-effects model was used, as true
#' effects can differ between studies when tau is greater than zero.
#' The same model is used across all simulation conditions for consistency.
#'
#' @param studies A data frame containing g and var_g.
#' @return A fitted random-effects meta-analysis model.
fit_meta_analysis <- function(studies) {
  
  # Check that the required columns are present
  if (!is.data.frame(studies) ||
      !all(c("g", "var_g") %in% names(studies))) {
    stop("Input must be a data frame containing g and var_g.")
  }
  
  # Check for invalid effect sizes or sampling variances
  if (nrow(studies) < 2 ||
      !is.numeric(studies$g) ||
      !is.numeric(studies$var_g) ||
      any(!is.finite(studies$g)) ||
      any(!is.finite(studies$var_g)) ||
      any(studies$var_g <= 0)) {
    stop("Invalid effect sizes, variances, or number of studies.")
  }
  
  # Use REML to estimate between-study heterogeneity
  # Use smaller steps and more iterations to improve convergence
  metafor::rma(
    yi = g,
    vi = var_g,
    data = studies,
    method = "REML",
    control = list(stepadj = 0.5, maxiter = 1000)
  )
}

# 2. Publication bias detection: Egger's test ------------------------------

# Egger's test checks for funnel-plot asymmetry, which may indicate
# small-study effects. It is used to examine whether publication bias
# can be detected in the simulated published literatures.

#' Run Egger's test for funnel-plot asymmetry
#'
#' Tests whether study effect sizes are associated with their standard
#' errors. We use this test to examine how often small-study effects
#' are detected across the simulated published literatures.
#'
#' @param model A fitted meta-analysis model from fit_meta_analysis().
#' @return The Egger's regression test result, including its p-value.
run_egger_test <- function(model) {
  
  # Check that the input is a fitted meta-analysis model
  if (!inherits(model, "rma")) {
    stop("Input must be a fitted metafor meta-analysis model.")
  }
  
  # Test for funnel-plot asymmetry using standard error as the predictor
  metafor::regtest(
    model,
    model = "rma",
    predictor = "sei"
  )
}

# 3. Publication bias correction: Trim-and-fill ----------------------------

#' Apply trim-and-fill to a meta-analysis
#'
#' Estimates potentially missing studies based on funnel-plot asymmetry
#' and returns an adjusted effect estimate. We compare this estimate
#' with the known true effect to see whether the correction reduces bias.
#'
#' @param model A fitted rma.uni model from fit_meta_analysis().
#' @return A trim-and-fill adjusted meta-analysis model.
run_trimfill <- function(model) {
  
  if (!inherits(model, "rma.uni")) {
    stop("Input must be a fitted metafor rma.uni model.")
  }
  
  metafor::trimfill(model)
}

# 4. Analyse simulated literatures -----------------------------------------

# Each literature is analysed using published studies and all studies.
# We also apply Egger's test and trim-and-fill to the published studies.
# The results are stored in one row per simulated literature.

#' Analyse one simulated literature
#'
#' Compares published-only and all-studies meta-analyses with the
#' known true effect. Also runs Egger's test and trim-and-fill.
#'
#' @param literature Data frame containing one condition and repetition
#'   from the simulated dataset.
#' @return A one-row data frame with effect estimates, estimation bias,
#'   Egger's test result, and number of imputed studies.
analyse_one_literature <- function(literature) {
  
  # Check that the required columns are present
  required <- c("condition_id", "rep", "true_d", "tau",
                "prob_publish_nonsig", "published", "g", "var_g")
  
  if (!is.data.frame(literature) ||
      !all(required %in% names(literature))) {
    stop("The literature is missing required columns.")
  }
  
  # Make sure the data contain only one condition and repetition
  settings <- c("condition_id", "rep", "true_d", "tau",
                "prob_publish_nonsig")
  
  if (any(vapply(literature[settings],
                 function(x) anyNA(x) || length(unique(x)) != 1,
                 logical(1)))) {
    stop("Input must contain exactly one simulated literature.")
  }
  
  if (!is.logical(literature$published) ||
      anyNA(literature$published)) {
    stop("The published column must contain TRUE/FALSE values.")
  }
  
  # Separate the published studies
  published_studies <- literature[literature$published, , drop = FALSE]
  
  # Fit both meta-analyses
  model_published <- fit_meta_analysis(published_studies)
  model_all <- fit_meta_analysis(literature)
  
  # Apply bias detection and correction
  # Keep the meta-analysis estimates if Egger's test fails
  egger <- tryCatch(
    run_egger_test(model_published),
    error = function(e) {
      warning(
        paste("Egger's test failed:", conditionMessage(e)),
        call. = FALSE
      )
      NULL
    }
  )
  
  # Keep the meta-analysis results if trim-and-fill fails
  corrected <- tryCatch(
    run_trimfill(model_published),
    error = function(e) {
      warning(
        paste("Trim-and-fill failed:", conditionMessage(e)),
        call. = FALSE
      )
      NULL
    }
  )
  
  # Extract the effect estimates
  est_published <- as.numeric(coef(model_published))
  est_all <- as.numeric(coef(model_all))
  est_corrected <- if (is.null(corrected)) {
    NA_real_
  } else {
    as.numeric(coef(corrected))
  }
  true_effect <- literature$true_d[1]
  egger_p_value <- if (is.null(egger)) {
    NA_real_
  } else {
    as.numeric(egger$pval)
  }
  # Store the results in one row
  data.frame(
    condition_id = literature$condition_id[1],
    rep = literature$rep[1],
    true_d = true_effect,
    tau = literature$tau[1],
    prob_publish_nonsig = literature$prob_publish_nonsig[1],
    n_all = nrow(literature),
    n_published = nrow(published_studies),
    estimate_published = est_published,
    estimate_all = est_all,
    estimate_trimfill = est_corrected,
    bias_published = est_published - true_effect,
    bias_all = est_all - true_effect,
    bias_trimfill = est_corrected - true_effect,
    egger_p = egger_p_value,
    egger_detected = if (is.na(egger_p_value)) {
      NA
    } else {
      egger_p_value < 0.05
    },
    egger_ok = !is.null(egger),
    trimfill_imputed = if (is.null(corrected)) {
      NA_integer_
    } else {
      corrected$k0
    },
    trimfill_ok = !is.null(corrected)
  )
}

#' Analyse all simulated literatures
#'
#' Separates the data by condition and repetition, then applies
#' analyse_one_literature() to each literature.
#' Errors and warnings are recorded so we can check which analyses
#' completed successfully.
#'
#' @param results Data frame produced by Marina's simulation.
#' @return A list containing the estimates and analysis diagnostics.
analyse_simulations <- function(results) {
  
  # Check the variables needed to identify each literature
  if (!is.data.frame(results) || nrow(results) == 0 ||
      !all(c("condition_id", "rep") %in% names(results)) ||
      anyNA(results$condition_id) || anyNA(results$rep)) {
    stop("Invalid simulation data or missing literature identifiers.")
  }
  
  # Each literature is identified by both condition_id and rep
  groups <- interaction(results$condition_id, results$rep,
                        drop = TRUE)
  literatures <- split(results, groups)
  
  analyses <- lapply(literatures, function(lit) {
    
    warning_messages <- character()
    
    # Keep track of warnings and errors without stopping the whole run
    result <- tryCatch(
      withCallingHandlers(
        analyse_one_literature(lit),
        warning = function(w) {
          warning_messages <<- c(warning_messages, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    
    if (inherits(result, "error")) {
      return(list(
        estimate = NULL,
        log = data.frame(
          condition_id = lit$condition_id[1],
          rep = lit$rep[1],
          status = "error",
          details = paste(c(warning_messages, conditionMessage(result)),
                          collapse = " | ")
        )
      ))
    }
    
    list(
      estimate = result,
      log = data.frame(
        condition_id = lit$condition_id[1],
        rep = lit$rep[1],
        status = if (length(warning_messages)) "warning" else "ok",
        details = if (length(warning_messages)) {
          paste(unique(warning_messages), collapse = " | ")
        } else {
          NA_character_
        }
      )
    )
  })
  
  # Combine successful estimates and keep all diagnostic records
  estimates <- Filter(Negate(is.null),
                      lapply(analyses, `[[`, "estimate"))
  
  estimates <- if (length(estimates)) {
    do.call(rbind, estimates)
  } else {
    NULL
  }
  
  diagnostics <- do.call(rbind, lapply(analyses, `[[`, "log"))
  
  if (!is.null(estimates)) rownames(estimates) <- NULL
  rownames(diagnostics) <- NULL
  
  list(estimates = estimates, diagnostics = diagnostics)
}

# 5. Summarise results by condition ----------------------------------------

# Each condition has 500 simulated literatures. We summarise their
# results to compare estimation bias and the performance of Egger's
# test and trim-and-fill across different simulation settings.

#' Summarise the analyses by simulation condition
#'
#' Calculates mean effect estimates, estimation bias, absolute errors,
#' and the proportion of significant Egger's tests for each condition.
#'
#' @param estimates Data frame from analyse_simulations()$estimates.
#' @return A data frame with one row per simulation condition.
summarise_conditions <- function(estimates) {
  
  # Check that the required columns are present
  required <- c(
    "condition_id", "rep", "true_d", "tau",
    "prob_publish_nonsig", "estimate_published",
    "estimate_all", "estimate_trimfill",
    "bias_published", "bias_all", "bias_trimfill",
    "egger_detected"
  )
  
  if (!is.data.frame(estimates) || nrow(estimates) == 0 ||
      !all(required %in% names(estimates))) {
    stop("Invalid estimates data or missing required columns.")
  }
  
  # Check that each literature appears only once
  if (anyNA(estimates$condition_id) ||
      anyNA(estimates$rep) ||
      anyDuplicated(estimates[c("condition_id", "rep")])) {
    stop("Missing or duplicated literature identifiers.")
  }
  
  # Group the results by condition
  conditions <- split(estimates, estimates$condition_id)
  
  summaries <- lapply(conditions, function(x) {
    
    # Check that the simulation settings are consistent
    settings <- c("true_d", "tau", "prob_publish_nonsig")
    
    if (any(vapply(x[settings], function(y) {
      anyNA(y) || length(unique(y)) != 1
    }, logical(1)))) {
      stop("Inconsistent simulation settings within a condition.")
    }
    
    valid_egger <- !is.na(x$egger_detected)
    valid_trimfill <- !is.na(x$bias_trimfill)
    
    # Calculate averages across repetitions
    data.frame(
      condition_id = x$condition_id[1],
      true_d = x$true_d[1],
      tau = x$tau[1],
      prob_publish_nonsig = x$prob_publish_nonsig[1],
      n_literatures = nrow(x),
      
      mean_est_published = mean(x$estimate_published),
      mean_est_all = mean(x$estimate_all),
      mean_est_trimfill = if (any(valid_trimfill)) {
        mean(x$estimate_trimfill[valid_trimfill])
      } else NA_real_,
      
      mean_bias_published = mean(x$bias_published),
      mean_bias_all = mean(x$bias_all),
      mean_bias_trimfill = if (any(valid_trimfill)) {
        mean(x$bias_trimfill[valid_trimfill])
      } else NA_real_,
      
      mean_abs_error_published = mean(abs(x$bias_published)),
      mean_abs_error_all = mean(abs(x$bias_all)),
      mean_abs_error_trimfill = if (any(valid_trimfill)) {
        mean(abs(x$bias_trimfill[valid_trimfill]))
      } else NA_real_,
      
      egger_significant_rate = if (any(valid_egger)) {
        mean(x$egger_detected[valid_egger])
      } else NA_real_,
      
      n_egger_valid = sum(valid_egger),
      n_trimfill_valid = sum(valid_trimfill)
    )
  })
  
  # Combine the condition summaries
  summary_table <- do.call(rbind, summaries)
  rownames(summary_table) <- NULL
  
  summary_table[order(summary_table$condition_id), ]
}
