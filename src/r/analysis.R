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

# 6. Visualisations --------------------------------------------------------

# The plots show how publication selection affects meta-analytic
# estimates and how well the detection and correction methods perform.

#' Plot estimation bias across simulation conditions
#'
#' Creates a heatmap showing the average bias in published-only
#' meta-analyses across different true effects, publication
#' probabilities, and levels of heterogeneity.
#'
#' @param condition_summary Data frame from summarise_conditions().
#' @return A ggplot heatmap.
plot_bias_heatmap <- function(condition_summary) {
  
  # Prepare labels for the simulation conditions
  plot_data <- condition_summary
  
  plot_data$selection <- factor(
    plot_data$prob_publish_nonsig,
    levels = c(1, 0.33, 0.10),
    labels = c("None", "Moderate", "Strong")
  )
  
  plot_data$true_effect <- factor(
    plot_data$true_d,
    levels = c(0, 0.24),
    labels = c("0.00", "0.24")
  )
  
  plot_data$heterogeneity <- factor(
    plot_data$tau,
    levels = c(0, 0.3),
    labels = c("tau = 0", "tau = 0.3")
  )
  
  # Create the heatmap
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = selection,
      y = true_effect,
      fill = mean_bias_published
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%+.3f", mean_bias_published)),
      size = 4
    ) +
    ggplot2::facet_wrap(~heterogeneity) +
    ggplot2::scale_fill_gradient2(
      low = "#EBDDFF",
      mid = "white",
      high = "#5A189A",
      midpoint = 0,
      name = "Mean bias"
    ) +
    ggplot2::labs(
      title = "Publication selection and estimation bias",
      subtitle = "Mean bias across 500 simulated literatures per condition",
      x = "Publication selection strength",
      y = "True effect (d)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#' Compare meta-analytic estimates with the true effect
#'
#' Shows how published-only, all-studies, and trim-and-fill
#' estimates change across simulation conditions.
#' Dashed lines indicate the known true effect.
#'
#' @param condition_summary Data frame from summarise_conditions().
#' @return A ggplot comparing mean effect estimates.
plot_estimate_comparison <- function(condition_summary) {
  
  # Check that the required columns are present
  required <- c(
    "true_d", "tau", "prob_publish_nonsig",
    "mean_est_published", "mean_est_all", "mean_est_trimfill"
  )
  
  if (!is.data.frame(condition_summary) ||
      !all(required %in% names(condition_summary))) {
    stop("Missing required columns for the comparison plot.")
  }
  
  # Put the three estimation methods into one data frame
  plot_data <- rbind(
    data.frame(
      condition_summary,
      method = "Published only",
      estimate = condition_summary$mean_est_published
    ),
    data.frame(
      condition_summary,
      method = "All studies",
      estimate = condition_summary$mean_est_all
    ),
    data.frame(
      condition_summary,
      method = "Trim-and-fill",
      estimate = condition_summary$mean_est_trimfill
    )
  )
  
  # Prepare labels
  plot_data$selection <- factor(
    plot_data$prob_publish_nonsig,
    levels = c(1, 0.33, 0.10),
    labels = c("None", "Moderate", "Strong")
  )
  
  plot_data$true_effect <- factor(
    plot_data$true_d,
    levels = c(0, 0.24),
    labels = c("True effect = 0", "True effect = 0.24")
  )
  
  plot_data$heterogeneity <- factor(
    plot_data$tau,
    levels = c(0, 0.3),
    labels = c("tau = 0", "tau = 0.3")
  )
  
  # One reference line for each panel
  reference_lines <- unique(
    plot_data[, c("true_effect", "heterogeneity", "true_d")]
  )
  
  # Create the comparison plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = selection,
      y = estimate,
      color = method,
      group = method
    )
  ) +
    ggplot2::geom_hline(
      data = reference_lines,
      ggplot2::aes(yintercept = true_d),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "grey40"
    ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_grid(true_effect ~ heterogeneity) +
    ggplot2::scale_color_manual(
      values = c(
        "Published only" = "#5A189A",
        "All studies" = "#8A8493",
        "Trim-and-fill" = "#B57EDC"
      )
    ) +
    ggplot2::labs(
      title = "Meta-analytic estimates across simulation conditions",
      subtitle = "Dashed lines show the known true effect",
      x = "Publication selection strength",
      y = "Mean effect estimate (Hedges' g)",
      color = "Analysis method"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  return(p)
}

#' Plot Egger's test significance rates
#'
#' Shows the proportion of significant Egger's tests across
#' publication selection conditions, true effects, and heterogeneity.
#' # The dashed line shows a 5% reference rate for comparison.
#'
#' @param condition_summary Data frame from summarise_conditions().
#' @return A ggplot showing Egger's test significance rates.
plot_egger_rates <- function(condition_summary) {
  
  # Check that the required columns are present
  required <- c(
    "true_d", "tau", "prob_publish_nonsig",
    "egger_significant_rate"
  )
  
  if (!is.data.frame(condition_summary) ||
      !all(required %in% names(condition_summary))) {
    stop("Missing required columns for the Egger's test plot.")
  }
  
  # Prepare labels for the simulation conditions
  plot_data <- condition_summary
  
  plot_data$selection <- factor(
    plot_data$prob_publish_nonsig,
    levels = c(1, 0.33, 0.10),
    labels = c("None", "Moderate", "Strong")
  )
  
  plot_data$true_effect <- factor(
    plot_data$true_d,
    levels = c(0, 0.24),
    labels = c("True effect = 0", "True effect = 0.24")
  )
  
  plot_data$heterogeneity <- factor(
    plot_data$tau,
    levels = c(0, 0.3),
    labels = c("tau = 0", "tau = 0.3")
  )
  
  # Create the plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = selection,
      y = egger_significant_rate
    )
  ) +
    ggplot2::geom_col(
      fill = "#7B2CBF",
      width = 0.65
    ) +
    ggplot2::geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      color = "grey40"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.1f%%", egger_significant_rate * 100)
      ),
      vjust = -0.5,
      size = 3.5
    ) +
    ggplot2::facet_grid(true_effect ~ heterogeneity) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(x * 100), "%"),
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      title = "Egger's test across simulation conditions",
      subtitle = "Dashed line shows the 5% reference rate",
      x = "Publication selection strength",
      y = "Proportion of significant tests"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.major.x = ggplot2::element_blank()
    )
  
  return(p)
}