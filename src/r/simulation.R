# =============================================================================
# simulation.R
#
# Simulates "Mozart effect" experiments and selective publication:
#   Step 1  set the truth        (arguments `true_d` and `tau`)
#   Step 2  run many studies     simulate_literature()
#   Step 3  journal filter       apply_publication_bias()
#   Step 4  published literature get_published_studies()
#
# Design modelled on the literature (Rauscher, Shaw & Ky, 1993; studies in
# the meta-analysis of Pietschnig, Voracek & Formann, 2010):
#   - Between-subjects: Mozart sonata K.448 (~10 min) vs. control (silence/rest)
#   - Outcome: spatial reasoning score (e.g. Stanford-Binet Paper Folding &
#     Cutting), expressed on an IQ-like scale (M = 100, SD = 15)
#   - Typical total N ~ 40 (median in Pietschnig et al., 2010), i.e. ~20/group
#
# This file only contains functions. It does not run anything on its own and
# does not set a random seed: call set.seed() once in main.rmd before
# simulating, so results are reproducible.
# =============================================================================


#' Check that a value is a single whole number of at least `min`
#'
#' @param x Value to check.
#' @param name Name of the argument (used in the error message).
#' @param min Smallest allowed value.
#' @return Nothing; throws an informative error if the check fails.
check_count <- function(x, name, min = 2) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) ||
      x != round(x) || x < min) {
    stop(sprintf("`%s` must be a single whole number >= %d (got: %s).",
                 name, min, paste(x, collapse = ", ")), call. = FALSE)
  }
}


#' Check that a value is a single finite number
#'
#' @param x Value to check.
#' @param name Name of the argument (used in the error message).
#' @param positive If TRUE, the value must also be > 0.
#' @return Nothing; throws an informative error if the check fails.
check_number <- function(x, name, positive = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) ||
      (positive && x <= 0)) {
    stop(sprintf("`%s` must be a single %snumber (got: %s).",
                 name, if (positive) "positive " else "",
                 paste(x, collapse = ", ")), call. = FALSE)
  }
}


#' Check that a value is a single probability between 0 and 1
#'
#' @param x Value to check.
#' @param name Name of the argument (used in the error message).
#' @return Nothing; throws an informative error if the check fails.
check_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0 || x > 1) {
    stop(sprintf("`%s` must be a single number between 0 and 1 (got: %s).",
                 name, paste(x, collapse = ", ")), call. = FALSE)
  }
}


# =============================================================================
# Building block: one single experiment
# =============================================================================

#' Generate participant-level data for one Mozart-effect experiment
#'
#' Scores in both groups are drawn from normal distributions with the same
#' SD. The Mozart group's mean is shifted by `true_d` standard deviations.
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
  check_count(n_per_group, "n_per_group")
  check_number(true_d, "true_d")
  check_number(control_mean, "control_mean")
  check_number(sd, "sd", positive = TRUE)
  
  mozart_mean <- control_mean + true_d * sd
  
  data.frame(
    participant = seq_len(2 * n_per_group),
    condition   = factor(rep(c("control", "mozart"), each = n_per_group),
                         levels = c("control", "mozart")),
    score       = c(rnorm(n_per_group, mean = control_mean, sd = sd),
                    rnorm(n_per_group, mean = mozart_mean,  sd = sd))
  )
}


#' Summarise one experiment as a meta-analysis would see it
#'
#' Computes the standardised mean difference (Hedges' g, i.e. Cohen's d with
#' small-sample correction), its standard error, and an independent-samples
#' t-test (equal variances, two-sided).
#'
#' @param data Data frame as returned by generate_study_data().
#' @return A one-row data frame with group sizes and means, Cohen's d,
#'   Hedges' g, the standard error and sampling variance of g, the t
#'   statistic, degrees of freedom and two-sided p-value.
summarise_study <- function(data) {
  required <- c("condition", "score")
  if (!is.data.frame(data) || !all(required %in% names(data))) {
    stop("`data` must be a data frame with columns `condition` and `score`.",
         call. = FALSE)
  }
  
  mozart  <- data$score[data$condition == "mozart"]
  control <- data$score[data$condition == "control"]
  n1 <- length(mozart)
  n2 <- length(control)
  if (n1 < 2 || n2 < 2) {
    stop("Each condition needs at least 2 participants.", call. = FALSE)
  }
  
  df        <- n1 + n2 - 2
  sd_pooled <- sqrt(((n1 - 1) * var(mozart) + (n2 - 1) * var(control)) / df)
  d         <- (mean(mozart) - mean(control)) / sd_pooled
  
  # Hedges' small-sample correction and the sampling variance of g
  # (Borenstein et al., 2009, Introduction to Meta-Analysis, ch. 4)
  j     <- 1 - 3 / (4 * df - 1)
  g     <- j * d
  var_g <- j^2 * ((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  
  # Student's t-test, identical to t.test(..., var.equal = TRUE)
  t_stat <- d / sqrt(1 / n1 + 1 / n2)
  p      <- 2 * pt(-abs(t_stat), df)
  
  data.frame(
    n_mozart     = n1,
    n_control    = n2,
    mean_mozart  = mean(mozart),
    mean_control = mean(control),
    d            = d,
    g            = g,
    se_g         = sqrt(var_g),
    var_g        = var_g,
    t            = t_stat,
    df           = df,
    p_value      = p
  )
}


#' Simulate one complete Mozart-effect experiment
#'
#' Convenience wrapper: generates participant data, then summarises it.
#'
#' @inheritParams generate_study_data
#' @return A one-row data frame (see summarise_study()) plus a column
#'   `true_d` recording the effect size used to generate the data.
simulate_mozart_study <- function(n_per_group = 20, true_d = 0,
                                  control_mean = 100, sd = 15) {
  data   <- generate_study_data(n_per_group, true_d, control_mean, sd)
  result <- summarise_study(data)
  cbind(true_d = true_d, result)
}


# =============================================================================
# Step 2: a whole literature of studies (all studies ever run)
# =============================================================================

#' Total sample sizes of the 38 studies in Pietschnig et al. (2010)
#'
#' Taken from the `mozart` dataset in the R package metaviz (column `n`).
#' Stored here so the simulation needs no extra packages. Median = 40.5,
#' range 8-206.
MOZART_TOTAL_N <- c(36, 24, 68, 92, 40, 8, 14, 32, 16, 38, 36, 12, 206, 86,
                    28, 28, 32, 67, 99, 90, 90, 145, 40, 40, 36, 26, 26, 51,
                    107, 66, 37, 160, 152, 63, 62, 59, 41, 86)


#' Simulate a literature of Mozart-effect experiments
#'
#' Simulates `n_studies` independent experiments, as if that many labs each
#' ran one study. This is the complete literature, including studies that
#' will later end up in the file drawer.
#'
#' Each study's sample size is drawn at random from `total_n_pool` (by
#' default, the sample sizes of the real Mozart literature) and split equally
#' over the two groups. Each study's true effect is drawn from
#' Normal(`true_d`, `tau`): with `tau = 0` every lab studies exactly the same
#' effect; with `tau > 0` true effects differ between labs (heterogeneity).
#'
#' @param n_studies Number of studies to simulate.
#' @param true_d Average true effect size (Cohen's d) across labs.
#' @param tau Between-study SD of true effects (heterogeneity). Default 0.
#' @param total_n_pool Numeric vector of total sample sizes to sample from.
#'   Pass a single number to give every study the same total N.
#' @param control_mean,sd Passed on to generate_study_data().
#' @return A data frame with one row per study: `study_id`, the average
#'   `true_d`, the study's own `study_true_d`, and all columns from
#'   summarise_study().
simulate_literature <- function(n_studies = 100, true_d = 0, tau = 0,
                                total_n_pool = MOZART_TOTAL_N,
                                control_mean = 100, sd = 15) {
  check_count(n_studies, "n_studies", min = 1)
  check_number(true_d, "true_d")
  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau < 0) {
    stop("`tau` must be a single number >= 0.", call. = FALSE)
  }
  if (!is.numeric(total_n_pool) || length(total_n_pool) < 1 ||
      any(!is.finite(total_n_pool)) || any(total_n_pool < 4)) {
    stop("`total_n_pool` must contain sample sizes of at least 4.",
         call. = FALSE)
  }
  
  # Note: sample(x, ...) behaves differently when x has length 1, so a
  # single sample size is repeated instead of sampled.
  total_n <- if (length(total_n_pool) == 1) {
    rep(total_n_pool, n_studies)
  } else {
    sample(total_n_pool, n_studies, replace = TRUE)
  }
  n_per_group  <- round(total_n / 2)
  study_true_d <- rnorm(n_studies, mean = true_d, sd = tau)
  
  studies <- lapply(seq_len(n_studies), function(i) {
    data <- generate_study_data(n_per_group[i], study_true_d[i],
                                control_mean, sd)
    summarise_study(data)
  })
  
  cbind(study_id     = seq_len(n_studies),
        true_d       = true_d,
        study_true_d = study_true_d,
        do.call(rbind, studies))
}


# =============================================================================
# Step 3: the journal filter (publication bias)
# =============================================================================

#' Decide which studies get published
#'
#' Models selective publication: studies with a significant result in the
#' predicted direction (p < alpha AND the Mozart group scored higher) are
#' always published. All other studies (non-significant, or significant in
#' the "wrong" direction) are published only with probability
#' `prob_publish_nonsig`.
#'
#' `prob_publish_nonsig = 1` means no publication bias; `0` means only
#' significant positive results ever appear. As a reference point, Franco,
#' Malhotra & Simonovits (2014, Science) found null results were published
#' about a third as often as strong results (~20% vs ~60%).
#'
#' @param literature Data frame as returned by simulate_literature().
#' @param prob_publish_nonsig Probability (0-1) that a study without a
#'   significant positive result is published.
#' @param alpha Significance threshold.
#' @return `literature` with two added columns: `significant` (TRUE if
#'   p < alpha and g > 0) and `published` (TRUE/FALSE). All studies are
#'   kept, so published and unpublished studies can be compared.
apply_publication_bias <- function(literature, prob_publish_nonsig = 0.33,
                                   alpha = 0.05) {
  if (!is.data.frame(literature) ||
      !all(c("g", "p_value") %in% names(literature))) {
    stop("`literature` must be a data frame from simulate_literature().",
         call. = FALSE)
  }
  check_probability(prob_publish_nonsig, "prob_publish_nonsig")
  check_probability(alpha, "alpha")
  
  significant <- literature$p_value < alpha & literature$g > 0
  lucky       <- runif(nrow(literature)) < prob_publish_nonsig
  
  literature$significant <- significant
  literature$published   <- significant | lucky
  literature
}


# =============================================================================
# Step 4: the published literature
# =============================================================================

#' Keep only the published studies
#'
#' @param literature Data frame as returned by apply_publication_bias().
#' @return The subset of rows with `published == TRUE`: what a
#'   meta-analyst searching the journals would find.
get_published_studies <- function(literature) {
  if (!"published" %in% names(literature)) {
    stop("Run apply_publication_bias() first.", call. = FALSE)
  }
  literature[literature$published, , drop = FALSE]
}