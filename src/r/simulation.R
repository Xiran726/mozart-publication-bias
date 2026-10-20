# =============================================================================
# simulation.R
#
# Simulates a literature of "Mozart effect" experiments, and the selective
# publication of their results, so that we can study how publication bias
# distorts meta-analytic effect estimates.
#
# Why simulate? In a real literature the true effect is unknown, so we can
# never tell how far a meta-analytic estimate is from the truth. In a
# simulation WE choose the true effect, so the distortion can be measured
# exactly: distortion = (meta-analytic estimate) - (true effect).
#
# PIPELINE (each step is one function; later steps call earlier ones)
#
#   1. One experiment         simulate_mozart_study()
#        Mozart vs. control, two groups; returns Hedges' g, SE, t, p.
#   2. Many experiments       simulate_literature()
#        Many labs each run one experiment = the complete literature,
#        including studies that will end up in the "file drawer".
#   3. Journal filter         apply_publication_bias()
#        Significant positive results are always published; all other
#        results only with probability `prob_publish_nonsig`.
#   4. Published literature   simulate_published_literature()
#        Labs keep running studies until `k_published` (default 38) are
#        published. Returns ALL studies run, flagged published TRUE/FALSE.
#   5. Replications           run_condition()
#        Repeats step 4 many times for one combination of settings.
#   6. All conditions         run_simulation()
#        Repeats step 5 for every combination of true effect x strength of
#        publication bias x heterogeneity, and returns one tidy data frame.
#   7. Save for the analysis  save_simulation()
#        Writes the results to data/simulated_literatures.rds.
#
# DESIGN OF ONE EXPERIMENT, modelled on the literature (Rauscher, Shaw & Ky,
# 1993; studies in the meta-analysis of Pietschnig, Voracek & Formann, 2010):
#   - Between-subjects: Mozart sonata K.448 (~10 min) vs. control (silence)
#   - Outcome: spatial reasoning score (e.g. Stanford-Binet Paper Folding &
#     Cutting task) on an IQ-like scale (M = 100, SD = 15)
#   - Sample sizes: drawn from the 38 studies in Pietschnig et al. (2010)
#
# REPRODUCIBILITY: the functions themselves never set a random seed, except
# run_simulation(), which takes an explicit `seed` argument. Run everything
# through run_simulation(seed = ...) and the results are identical every time.
#
# No packages are needed: everything uses base R.
#
# References
#   Franco, A., Malhotra, N., & Simonovits, G. (2014). Publication bias in the
#     social sciences: Unlocking the file drawer. Science, 345, 1502-1505.
#   Pietschnig, J., Voracek, M., & Formann, A. K. (2010). Mozart effect-
#     Shmozart effect: A meta-analysis. Intelligence, 38, 314-323.
#   Rauscher, F. H., Shaw, G. L., & Ky, K. N. (1993). Music and spatial task
#     performance. Nature, 365, 611.
#   Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R.
#     (2009). Introduction to Meta-Analysis. Wiley.
# =============================================================================


# =============================================================================
# Input checks (used by all functions below to give informative errors)
# =============================================================================

#' Check that a value is a single whole number of at least `min`
#'
#' @param x Value to check.
#' @param name Name of the argument (used in the error message).
#' @param min Smallest allowed value.
#' @return Nothing; throws an informative error if the check fails.
check_count <- function(x, name, min = 1) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) ||
      x != round(x) || x < min) {
    stop(sprintf("`%s` must be a single whole number >= %s (got: %s).",
                 name, format(min), paste(x, collapse = ", ")),
         call. = FALSE)
  }
}


#' Check that a value is a single finite number of at least `min`
#'
#' @param x Value to check.
#' @param name Name of the argument (used in the error message).
#' @param min Smallest allowed value (default: no lower limit).
#' @param strict If TRUE, the value must be strictly greater than `min`.
#' @return Nothing; throws an informative error if the check fails.
check_number <- function(x, name, min = -Inf, strict = FALSE) {
  valid <- is.numeric(x) && length(x) == 1 && is.finite(x) &&
    (if (strict) x > min else x >= min)
  if (!valid) {
    limit <- if (is.finite(min)) {
      sprintf(" %s %s", if (strict) ">" else ">=", format(min))
    } else {
      ""
    }
    stop(sprintf("`%s` must be a single number%s (got: %s).",
                 name, limit, paste(x, collapse = ", ")), call. = FALSE)
  }
}


#' Check that a value is a single probability between 0 and 1
#'
#' @param x Value to check.
#' @param name Name of the argument (used in the error message).
#' @return Nothing; throws an informative error if the check fails.
check_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || x > 1) {
    stop(sprintf("`%s` must be a single number between 0 and 1 (got: %s).",
                 name, paste(x, collapse = ", ")), call. = FALSE)
  }
}


#' Check that every value in a vector is a valid probability
#'
#' @param x Vector to check.
#' @param name Name of the argument (used in the error message).
#' @return Nothing; throws an informative error if the check fails.
check_probabilities <- function(x, name) {
  if (!is.numeric(x) || length(x) < 1 || any(!is.finite(x)) ||
      any(x < 0) || any(x > 1)) {
    stop(sprintf("`%s` must contain numbers between 0 and 1.", name),
         call. = FALSE)
  }
}


# =============================================================================
# Step 1: one single experiment
# =============================================================================

#' Generate participant-level data for one Mozart-effect experiment
#'
#' Scores in both groups are drawn from normal distributions with the same
#' SD. The Mozart group's population mean is `true_d` standard deviations
#' higher than the control group's, so `true_d` is the true Cohen's d.
#' Example: with sd = 15, true_d = 0.24 means a true advantage of 3.6 IQ
#' points; Rauscher et al. (1993) reported 8-9 points (d of about 0.55).
#'
#' @param n_per_group Number of participants in each group.
#' @param true_d True population effect size (Cohen's d) of listening to
#'   Mozart vs. control. 0 means "no Mozart effect".
#' @param control_mean Population mean of the control group (IQ scale).
#' @param sd Population SD of scores in both groups (IQ scale).
#' @return A data frame with one row per participant and columns
#'   `participant`, `condition` ("control" / "mozart") and `score`.
generate_study_data <- function(n_per_group = 20, true_d = 0,
                                control_mean = 100, sd = 15) {
  check_count(n_per_group, "n_per_group", min = 2)
  check_number(true_d, "true_d")
  check_number(control_mean, "control_mean")
  check_number(sd, "sd", min = 0, strict = TRUE)
  
  mozart_mean <- control_mean + true_d * sd
  
  data.frame(
    participant = seq_len(2 * n_per_group),
    condition   = factor(rep(c("control", "mozart"), each = n_per_group),
                         levels = c("control", "mozart")),
    score       = c(rnorm(n_per_group, mean = control_mean, sd = sd),
                    rnorm(n_per_group, mean = mozart_mean,  sd = sd))
  )
}


#' Compute the effect size and significance test for two groups of scores
#'
#' This is the core calculation behind every simulated study: everything a
#' meta-analysis needs from one study.
#'   - Cohen's d: difference in means divided by the pooled SD.
#'   - Hedges' g: d multiplied by a small-sample correction J, because d
#'     slightly overestimates the true effect in small samples
#'     (with 20 per group, J is about 0.98).
#'   - var_g / se_g: sampling variance and standard error of g. Meta-analyses
#'     use these to weight studies; funnel plots and Egger's test use se_g.
#'   - t and p: Student's independent-samples t-test (two-sided, equal
#'     variances), identical to t.test(..., var.equal = TRUE) but faster.
#' Formulas: Borenstein et al. (2009), chapter 4.
#'
#' @param mozart Numeric vector of scores in the Mozart group.
#' @param control Numeric vector of scores in the control group.
#' @return A named numeric vector: n_mozart, n_control, mean_mozart,
#'   mean_control, d, g, var_g, se_g, t, df, p_value.
compute_effect_size <- function(mozart, control) {
  if (!is.numeric(mozart) || !is.numeric(control)) {
    stop("`mozart` and `control` must be numeric vectors.", call. = FALSE)
  }
  n1 <- length(mozart)
  n2 <- length(control)
  if (n1 < 2 || n2 < 2) {
    stop("Each group needs at least 2 participants.", call. = FALSE)
  }
  
  df        <- n1 + n2 - 2
  sd_pooled <- sqrt(((n1 - 1) * var(mozart) + (n2 - 1) * var(control)) / df)
  d         <- (mean(mozart) - mean(control)) / sd_pooled
  
  j     <- 1 - 3 / (4 * df - 1)                       # Hedges' correction
  g     <- j * d
  var_g <- j^2 * ((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  
  t_stat <- d / sqrt(1 / n1 + 1 / n2)
  p      <- 2 * pt(-abs(t_stat), df)
  
  c(n_mozart = n1, n_control = n2,
    mean_mozart = mean(mozart), mean_control = mean(control),
    d = d, g = g, var_g = var_g, se_g = sqrt(var_g),
    t = t_stat, df = df, p_value = p)
}


#' Summarise one experiment as a meta-analysis would see it
#'
#' @param data Data frame as returned by generate_study_data().
#' @return A one-row data frame with the results of compute_effect_size().
summarise_study <- function(data) {
  if (!is.data.frame(data) || !all(c("condition", "score") %in% names(data))) {
    stop("`data` must be a data frame with columns `condition` and `score`.",
         call. = FALSE)
  }
  result <- compute_effect_size(
    mozart  = data$score[data$condition == "mozart"],
    control = data$score[data$condition == "control"]
  )
  as.data.frame(as.list(result))
}


#' Simulate one complete Mozart-effect experiment
#'
#' Convenience wrapper: generates participant data, then summarises it.
#' Useful to illustrate a single study; the pipeline below uses the same
#' building blocks.
#'
#' @inheritParams generate_study_data
#' @return A one-row data frame (see summarise_study()) plus a column
#'   `true_d` recording the effect size used to generate the data.
simulate_mozart_study <- function(n_per_group = 20, true_d = 0,
                                  control_mean = 100, sd = 15) {
  data <- generate_study_data(n_per_group, true_d, control_mean, sd)
  cbind(true_d = true_d, summarise_study(data))
}


# =============================================================================
# Step 2: many experiments (the complete literature, file drawer included)
# =============================================================================

#' Total sample sizes of the 38 studies in Pietschnig et al. (2010)
#'
#' Taken from the `mozart` dataset in the R package metaviz (column `n`).
#' Stored here so the simulation needs no extra packages.
#' Median = 40.5 (about 20 per group), range 8-206.
#' Using real sample sizes matters: bias-detection methods such as Egger's
#' test look for small studies reporting larger effects, so the simulated
#' studies need a realistic mix of small and large samples.
MOZART_TOTAL_N <- c(36, 24, 68, 92, 40, 8, 14, 32, 16, 38, 36, 12, 206, 86,
                    28, 28, 32, 67, 99, 90, 90, 145, 40, 40, 36, 26, 26, 51,
                    107, 66, 37, 160, 152, 63, 62, 59, 41, 86)


#' Simulate a literature of Mozart-effect experiments
#'
#' Simulates `n_studies` independent experiments, as if that many labs each
#' ran one study.
#'
#' Sample size: each study's total N is drawn at random from `total_n_pool`
#' (by default the real Mozart literature) and split equally over the groups.
#'
#' Heterogeneity: each study's own true effect is drawn from
#' Normal(`true_d`, `tau`). With tau = 0 every lab studies exactly the same
#' effect. With tau > 0 true effects differ between labs (e.g. because of
#' different tasks, music or participants), as in most real literatures.
#'
#' @param n_studies Number of studies to simulate.
#' @param true_d Average true effect size (Cohen's d) across labs.
#' @param tau Between-study SD of true effects (heterogeneity).
#' @param total_n_pool Numeric vector of total sample sizes to sample from.
#'   Pass a single number to give every study the same total N.
#' @param con