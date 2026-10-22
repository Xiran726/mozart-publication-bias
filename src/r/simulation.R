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
#' @param n_mozart Number of participants in the Mozart group.
#' @param n_control Number of participants in the control group (by default
#'   the same as `n_mozart`). The groups may differ in size, so that studies
#'   with an odd total sample size (e.g. 67 = 34 + 33) can be simulated.
#' @param true_d True population effect size (Cohen's d) of listening to
#'   Mozart vs. control. 0 means "no Mozart effect".
#' @param control_mean Population mean of the control group (IQ scale).
#' @param sd Population SD of scores in both groups (IQ scale).
#' @return A data frame with one row per participant and columns
#'   `participant`, `condition` ("control" / "mozart") and `score`.
generate_study_data <- function(n_mozart = 20, n_control = n_mozart,
                                true_d = 0, control_mean = 100, sd = 15) {
  check_count(n_mozart, "n_mozart", min = 2)
  check_count(n_control, "n_control", min = 2)
  check_number(true_d, "true_d")
  check_number(control_mean, "control_mean")
  check_number(sd, "sd", min = 0, strict = TRUE)
  
  mozart_mean <- control_mean + true_d * sd
  
  data.frame(
    participant = seq_len(n_control + n_mozart),
    condition   = factor(rep(c("control", "mozart"),
                             times = c(n_control, n_mozart)),
                         levels = c("control", "mozart")),
    score       = c(rnorm(n_control, mean = control_mean, sd = sd),
                    rnorm(n_mozart,  mean = mozart_mean,  sd = sd))
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
  if (any(!is.finite(mozart)) || any(!is.finite(control))) {
    stop("All scores must be finite numbers (no NA, NaN or Inf).",
         call. = FALSE)
  }
  
  df        <- n1 + n2 - 2
  sd_pooled <- sqrt(((n1 - 1) * var(mozart) + (n2 - 1) * var(control)) / df)
  # If every score is identical, the pooled SD is 0 and d would be a
  # division by zero, so the effect size is undefined
  if (sd_pooled <= 0) {
    stop("The pooled SD is 0 (all scores are identical), ",
         "so the effect size is undefined.", call. = FALSE)
  }
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
  # Every row must be labelled exactly "mozart" or "control"; otherwise rows
  # with a typo (e.g. "Mozart") would be silently left out of the analysis
  labels <- as.character(data$condition)
  unexpected <- unique(labels[is.na(labels) |
                                !labels %in% c("mozart", "control")])
  if (length(unexpected) > 0) {
    stop(sprintf(paste("`condition` may only contain \"mozart\" and",
                       "\"control\" (found: %s)."),
                 paste(unexpected, collapse = ", ")), call. = FALSE)
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
simulate_mozart_study <- function(n_mozart = 20, n_control = n_mozart,
                                  true_d = 0, control_mean = 100, sd = 15) {
  data <- generate_study_data(n_mozart, n_control, true_d, control_mean, sd)
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
#' (by default the real Mozart literature) and split as evenly as possible
#' over the two groups. The exact total N is kept: an odd total such as 67
#' becomes 34 (Mozart) + 33 (control), so every sample size in the real
#' literature can occur. (Which group gets the extra participant does not
#' matter for the results; we give it to the Mozart group.)
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
#' @param control_mean,sd Passed on to generate_study_data().
#' @return A data frame with one row per study: `study_id`, the study's own
#'   true effect `study_true_d`, and the results of compute_effect_size().
simulate_literature <- function(n_studies = 100, true_d = 0, tau = 0,
                                total_n_pool = MOZART_TOTAL_N,
                                control_mean = 100, sd = 15) {
  check_count(n_studies, "n_studies", min = 1)
  check_number(true_d, "true_d")
  check_number(tau, "tau", min = 0)
  if (!is.numeric(total_n_pool) || length(total_n_pool) < 1 ||
      any(!is.finite(total_n_pool)) || any(total_n_pool < 4) ||
      any(total_n_pool != round(total_n_pool))) {
    stop("`total_n_pool` must contain whole-number sample sizes of at least 4.",
         call. = FALSE)
  }
  
  # Draw each study's sample size. (sample() behaves differently when given
  # a single number, so a single sample size is repeated instead.)
  total_n <- if (length(total_n_pool) == 1) {
    rep(total_n_pool, n_studies)
  } else {
    sample(total_n_pool, n_studies, replace = TRUE)
  }
  n_control    <- floor(total_n / 2)    # e.g. 67 -> 33
  n_mozart     <- total_n - n_control   # e.g. 67 -> 34 (exact total kept)
  study_true_d <- rnorm(n_studies, mean = true_d, sd = tau)
  
  # Run every study: generate participants, then compute its effect size.
  # vapply() returns a matrix with one column per study, so t() turns it
  # into one row per study.
  results <- vapply(seq_len(n_studies), function(i) {
    data <- generate_study_data(n_mozart[i], n_control[i], study_true_d[i],
                                control_mean, sd)
    compute_effect_size(mozart  = data$score[data$condition == "mozart"],
                        control = data$score[data$condition == "control"])
  }, FUN.VALUE = numeric(11))
  
  results <- as.data.frame(t(results))
  names(results) <- c("n_mozart", "n_control", "mean_mozart", "mean_control",
                      "d", "g", "var_g", "se_g", "t", "df", "p_value")
  
  cbind(study_id = seq_len(n_studies), study_true_d = study_true_d, results)
}


# =============================================================================
# Step 3: the journal filter (publication bias)
# =============================================================================

#' Decide which studies get published
#'
#' Models selective publication. A study is "significant_positive" when
#' p < alpha AND the Mozart group scored higher (the direction that supports
#' the Mozart effect). These studies are always published. All other
#' studies (non-significant, or significant in the "wrong" direction) are
#' published only with
#' probability `prob_publish_nonsig`.
#'
#' Choosing `prob_publish_nonsig`:
#'   1    = no publication bias (everything is published)
#'   0.33 = realistic: Franco et al. (2014) found that null results were
#'          published about a third as often as strong results (~20% vs ~60%)
#'   0.10 = strong publication bias
#'   0    = only significant positive results are ever published
#'
#' @param literature Data frame as returned by simulate_literature().
#' @param prob_publish_nonsig Probability (0-1) that a study without a
#'   significant positive result is published.
#' @param alpha Significance threshold.
#' @return `literature` with two added columns: `significant_positive`
#'   (TRUE/FALSE)
#'   and `published` (TRUE/FALSE). All studies are kept, so published and
#'   unpublished studies can be compared.
apply_publication_bias <- function(literature, prob_publish_nonsig = 0.33,
                                   alpha = 0.05) {
  if (!is.data.frame(literature) ||
      !all(c("g", "p_value") %in% names(literature))) {
    stop("`literature` must be a data frame from simulate_literature().",
         call. = FALSE)
  }
  check_probability(prob_publish_nonsig, "prob_publish_nonsig")
  check_probability(alpha, "alpha")
  
  significant_positive <- literature$p_value < alpha & literature$g > 0
  lucky       <- runif(nrow(literature)) < prob_publish_nonsig
  
  literature$significant_positive <- significant_positive
  literature$published            <- significant_positive | lucky
  literature
}


#' Keep only the published studies
#'
#' @param literature Data frame with a `published` column, e.g. from
#'   apply_publication_bias() or run_simulation().
#' @return The subset of rows with `published == TRUE`: what a meta-analyst
#'   searching the journals would find.
get_published_studies <- function(literature) {
  if (!is.data.frame(literature) || !"published" %in% names(literature)) {
    stop("`literature` needs a `published` column; ",
         "run apply_publication_bias() first.", call. = FALSE)
  }
  literature[literature$published, , drop = FALSE]
}


# =============================================================================
# Step 4: keep running studies until k are published
# =============================================================================

#' Simulate a literature until a fixed number of studies is published
#'
#' Labs keep running Mozart experiments one after another, and every
#' finished study goes through the journal filter. This stops as soon as
#' `k_published` studies have been published.
#'
#' Why a fixed number of PUBLISHED studies? It keeps the meta-analyses
#' comparable across conditions: every meta-analysis pools the same number of
#' studies, so differences between conditions reflect publication bias rather
#' than how many studies happened to be published. The default of 38 matches
#' the meta-analysis of Pietschnig et al. (2010). The stronger the bias, the
#' more studies end up in the file drawer before 38 are published.
#'
#' For speed, studies are simulated in batches. Studies run after the
#' `k_published`-th publication are discarded, so the result is the same as
#' running studies one at a time and stopping at the k-th publication.
#' A batch is never larger than the number of studies still allowed by
#' `max_studies`, so `max_studies` is a strict upper limit.
#'
#' @param k_published Number of published studies to collect.
#' @param true_d,tau,total_n_pool,control_mean,sd Passed on to
#'   simulate_literature().
#' @param prob_publish_nonsig,alpha Passed on to apply_publication_bias().
#' @param batch_size Number of studies simulated per batch.
#' @param max_studies Safety limit on the total number of studies run. Stops
#'   with an error if `k_published` is not reached by then (e.g. when almost
#'   no study is ever published). Must be at least `k_published`, since no
#'   more studies can be published than are run.
#' @return A data frame with ALL studies that were run (published AND file
#'   drawer), in the order they were run. Exactly `k_published` rows have
#'   `published == TRUE`.
simulate_published_literature <- function(k_published = 38, true_d = 0,
                                          tau = 0, prob_publish_nonsig = 0.33,
                                          alpha = 0.05,
                                          total_n_pool = MOZART_TOTAL_N,
                                          control_mean = 100, sd = 15,
                                          batch_size = 50,
                                          max_studies = 1e5) {
  check_count(k_published, "k_published", min = 1)
  check_count(batch_size, "batch_size", min = 1)
  check_count(max_studies, "max_studies", min = 1)
  if (k_published > max_studies) {
    stop(sprintf(paste("`k_published` (%d) cannot be larger than",
                       "`max_studies` (%d): you cannot publish more studies",
                       "than are run."),
                 as.integer(k_published), as.integer(max_studies)),
         call. = FALSE)
  }
  
  batches     <- list()
  n_published <- 0
  n_run       <- 0
  
  # Keep simulating batches of studies until enough have been published
  while (n_published < k_published) {
    if (n_run >= max_studies) {
      stop(sprintf(paste("Only %d of %d studies were published after %d",
                         "studies were run. Increase `max_studies` or",
                         "`prob_publish_nonsig`."),
                   as.integer(n_published), as.integer(k_published),
                   as.integer(n_run)), call. = FALSE)
    }
    # Never run more studies than `max_studies` still allows
    this_batch <- min(batch_size, max_studies - n_run)
    batch <- simulate_literature(this_batch, true_d, tau, total_n_pool,
                                 control_mean, sd)
    batch <- apply_publication_bias(batch, prob_publish_nonsig, alpha)
    
    batches[[length(batches) + 1]] <- batch
    n_published <- n_published + sum(batch$published)
    n_run       <- n_run + this_batch
  }
  
  literature <- do.call(rbind, batches)
  
  # Drop the studies that were run after the k-th publication
  last_study <- which(cumsum(literature$published) == k_published)[1]
  literature <- literature[seq_len(last_study), ]
  
  literature$study_id <- seq_len(last_study)   # renumber 1, 2, 3, ...
  rownames(literature) <- NULL
  literature
}


# =============================================================================
# Step 5: many replications of one condition
# =============================================================================

#' Simulate many literatures under the same settings
#'
#' One simulated literature is a single random draw: its meta-analytic
#' estimate could be too high or too low by chance. Repeating the simulation
#' `n_reps` times shows what happens ON AVERAGE under these settings, and how
#' much results vary.
#'
#' @param n_reps Number of replications (simulated literatures).
#' @param ... Settings passed on to simulate_published_literature(), e.g.
#'   `true_d`, `tau`, `prob_publish_nonsig`, `k_published`.
#' @return One data frame with all studies of all replications; the column
#'   `rep` (1, 2, ..., n_reps) says which literature a study belongs to.
run_condition <- function(n_reps = 500, ...) {
  check_count(n_reps, "n_reps", min = 1)
  
  literatures <- lapply(seq_len(n_reps), function(r) {
    literature <- simulate_published_literature(...)
    cbind(rep = r, literature)
  })
  do.call(rbind, literatures)
}


# =============================================================================
# Step 6: all conditions
# =============================================================================

#' Create the grid of simulation conditions
#'
#' Every combination of the three factors becomes one condition. The
#' defaults give 2 x 3 x 2 = 12 conditions:
#'   - true_d: 0 (the Mozart effect does not exist) and 0.24 (a small real
#'     effect: our random-effects estimate for studies outside the Rauscher
#'     and Rideout labs in the Pietschnig et al., 2010, data)
#'   - prob_publish_nonsig: 1 (no bias), 0.33 (realistic), 0.10 (strong)
#'   - tau: 0 (same effect in every lab) and 0.3 (effects differ between
#'     labs; close to the heterogeneity in the real Mozart literature)
#'
#' @param true_d Vector of average true effects.
#' @param prob_publish_nonsig Vector of publication probabilities for
#'   non-significant results.
#' @param tau Vector of heterogeneity values.
#' @return A data frame with one row per condition and columns
#'   `condition_id`, `true_d`, `tau`, `prob_publish_nonsig`.
make_conditions <- function(true_d = c(0, 0.24),
                            prob_publish_nonsig = c(1, 0.33, 0.10),
                            tau = c(0, 0.3)) {
  if (!is.numeric(true_d) || length(true_d) < 1 || any(!is.finite(true_d))) {
    stop("`true_d` must contain finite numbers.", call. = FALSE)
  }
  if (!is.numeric(tau) || length(tau) < 1 || any(!is.finite(tau)) ||
      any(tau < 0)) {
    stop("`tau` must contain numbers >= 0.", call. = FALSE)
  }
  check_probabilities(prob_publish_nonsig, "prob_publish_nonsig")
  
  conditions <- expand.grid(true_d = true_d, tau = tau,
                            prob_publish_nonsig = prob_publish_nonsig)
  cbind(condition_id = seq_len(nrow(conditions)), conditions)
}


#' Columns of each simulated study that are kept in the final results
RESULT_COLUMNS <- c("rep", "study_id", "study_true_d", "n_mozart",
                    "n_control", "g", "var_g", "se_g", "t", "p_value",
                    "significant_positive", "published")


#' Run the full simulation
#'
#' For every condition in `conditions`, simulates `n_reps` literatures with
#' `k_published` published studies each, and stacks everything into one
#' tidy data frame (one row per simulated study).
#'
#' @param conditions Data frame from make_conditions().
#' @param n_reps Number of replications (literatures) per condition.
#' @param k_published Number of published studies per literature.
#' @param seed Random seed, so results are reproducible. Use NULL to not set
#'   a seed.
#' @param verbose If TRUE, prints progress messages.
#' @return A data frame with one row per simulated study and columns:
#'   condition_id, true_d, tau,
#'     prob_publish_nonsig  the condition this study belongs to
#'   rep                    which simulated literature (1 ... n_reps)
#'   study_id               order in which studies in that literature were run
#'   study_true_d           this study's own true effect (= true_d if tau = 0)
#'   n_mozart, n_control    group sizes
#'   g, var_g, se_g         Hedges' g, its sampling variance and standard error
#'   t, p_value             t-test result
#'   significant_positive   p < .05 AND Mozart group higher
#'   published              passed the journal filter (FALSE = file drawer)
#'   The settings used are stored in attr(result, "settings").
run_simulation <- function(conditions = make_conditions(), n_reps = 500,
                           k_published = 38, seed = 2026, verbose = TRUE) {
  needed <- c("condition_id", "true_d", "tau", "prob_publish_nonsig")
  if (!is.data.frame(conditions) || !all(needed %in% names(conditions))) {
    stop("`conditions` must be a data frame from make_conditions().",
         call. = FALSE)
  }
  if (nrow(conditions) == 0) {
    stop("`conditions` is empty: it needs at least one condition.",
         call. = FALSE)
  }
  # condition_id identifies a condition in the results; duplicated or
  # missing IDs would mix up studies from different conditions in the
  # analysis
  if (anyNA(conditions$condition_id) ||
      anyDuplicated(conditions$condition_id) > 0) {
    stop("`condition_id` must be unique and must not be missing.",
         call. = FALSE)
  }
  check_count(n_reps, "n_reps", min = 1)
  check_count(k_published, "k_published", min = 1)
  if (!is.null(seed)) {
    check_count(seed, "seed", min = 0)
    set.seed(seed)
  }
  
  start_time <- Sys.time()
  results <- vector("list", nrow(conditions))
  
  for (i in seq_len(nrow(conditions))) {
    cond <- conditions[i, ]
    if (verbose) {
      message(sprintf(
        "Condition %d of %d: true_d = %.2f, tau = %.2f, prob_publish_nonsig = %.2f",
        i, nrow(conditions), cond$true_d, cond$tau, cond$prob_publish_nonsig))
    }
    studies <- run_condition(
      n_reps              = n_reps,
      k_published         = k_published,
      true_d              = cond$true_d,
      tau                 = cond$tau,
      prob_publish_nonsig = cond$prob_publish_nonsig
    )
    # Put the condition settings in front of every study, and keep only the
    # columns the analysis needs (group means, Cohen's d and df are left out
    # to keep the saved file small; g, var_g and se_g contain what is needed)
    results[[i]] <- cbind(cond[rep(1, nrow(studies)), needed],
                          studies[, RESULT_COLUMNS])
  }
  
  results <- do.call(rbind, results)
  rownames(results) <- NULL
  
  attr(results, "settings") <- list(n_reps = n_reps,
                                    k_published = k_published,
                                    seed = seed,
                                    created = Sys.time(),
                                    r_version = R.version.string)
  if (verbose) {
    message(sprintf("Done: %d studies simulated in %.1f minutes.",
                    nrow(results),
                    as.numeric(difftime(Sys.time(), start_time,
                                        units = "mins"))))
  }
  results
}


#' Describe what happened in each condition (a quick sanity check)
#'
#' Simple descriptive summaries per condition. NOT a meta-analysis: means
#' are unweighted. The proper meta-analytic estimates and bias tests are
#' part of the analysis step.
#'
#' @param results Data frame from run_simulation().
#' @return A data frame with one row per condition: the condition settings,
#'   the average number of studies run per literature, the share of studies
#'   published, the mean g of all studies, and the mean g of published
#'   studies only.
describe_simulation <- function(results) {
  needed <- c("condition_id", "true_d", "tau", "prob_publish_nonsig",
              "rep", "g", "published")
  if (!is.data.frame(results) || !all(needed %in% names(results))) {
    stop("`results` must be a data frame from run_simulation().",
         call. = FALSE)
  }
  
  by_condition <- split(results, results$condition_id)
  rows <- lapply(by_condition, function(x) {
    data.frame(
      condition_id          = x$condition_id[1],
      true_d                = x$true_d[1],
      tau                   = x$tau[1],
      prob_publish_nonsig   = x$prob_publish_nonsig[1],
      mean_studies_run      = nrow(x) / length(unique(x$rep)),
      prop_published        = mean(x$published),
      mean_g_all            = mean(x$g),
      mean_g_published      = mean(x$g[x$published])
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# =============================================================================
# Step 7: save the results for the analysis
# =============================================================================

#' Save simulation results to an .rds file
#'
#' The full simulation takes a few minutes, so it is run once and saved.
#' The analysis can then load the results with readRDS() instead of
#' re-running the simulation.
#'
#' @param results Data frame from run_simulation().
#' @param file Path of the .rds file to write. Missing folders are created.
#' @return The file path (invisibly).
save_simulation <- function(results,
                            file = file.path("data",
                                             "simulated_literatures.rds")) {
  if (!is.data.frame(results)) {
    stop("`results` must be a data frame from run_simulation().",
         call. = FALSE)
  }
  if (!is.character(file) || length(file) != 1 || !nzchar(file)) {
    stop("`file` must be a single file path.", call. = FALSE)
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, file)
  invisible(file)
}