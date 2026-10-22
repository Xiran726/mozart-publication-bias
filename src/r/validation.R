# Validate the study-level output of run_simulation().
# This module uses base R and does not require sourcing simulation.R.

#' Validate simulated literature data before analysis
#'
#' Checks the output structure, numerical values, study identifiers,
#' condition consistency, and optionally the published count per literature.
#' Does not modify results or verify the statistical simulation algorithm.
#'
#' @param results Non-empty data frame returned by run_simulation().
#' @param expected_k Optional positive integer: required number of published
#'   studies in each condition_id x rep. NULL skips this count requirement.
#' @param tolerance Positive numerical tolerance for numerical consistency checks.
#' @param alpha Significance threshold used to generate the publication flags.
#'   Defaults to 0.05, matching run_simulation().
#' @return TRUE invisibly when valid; otherwise stops with an informative error.
validate_simulation_results <- function(results, expected_k = NULL,
                                        tolerance = 1e-8, alpha = 0.05) {
  fail <- function(message) stop(message, call. = FALSE)
  
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) {
    fail("`tolerance` must be a single positive finite number.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    fail("`alpha` must be a single finite number between 0 and 1.")
  }
  if (!is.null(expected_k)) {
    if (!is.numeric(expected_k) || length(expected_k) != 1L ||
        !is.finite(expected_k) || expected_k < 1 ||
        expected_k != floor(expected_k)) {
      fail("`expected_k` must be NULL or a single positive whole number.")
    }
  }
  if (!is.data.frame(results) || nrow(results) == 0L) {
    fail("`results` must be a non-empty data frame.")
  }
  if (anyDuplicated(names(results)) > 0L) {
    fail("`results` must not contain duplicate column names.")
  }
  
  # The expected outcome variables from simulation script
  required <- c(
    "condition_id", "true_d", "tau", "prob_publish_nonsig", "rep",
    "study_id", "study_true_d", "n_mozart", "n_control", "g",
    "var_g", "se_g", "t", "p_value", "significant_positive", "published"
  )
  missing_columns <- setdiff(required, names(results))
  if (length(missing_columns) > 0L) {
    fail(paste0("Missing required columns: ",
                paste(missing_columns, collapse = ", "), "."))
  }
  
  numeric_columns <- setdiff(required, c("significant_positive", "published"))
  for (column in numeric_columns) {
    values <- results[[column]]
    if (!is.numeric(values) || !is.null(dim(values))) {
      fail(sprintf("Column `%s` must be a numeric vector.", column))
    }
    bad <- which(!is.finite(values))
    if (length(bad) > 0L) {
      fail(sprintf("Column `%s` contains a missing or non-finite value at row %d.",
                   column, bad[1L]))
    }
  }
  for (column in c("significant_positive", "published")) {
    values <- results[[column]]
    if (!is.logical(values) || !is.null(dim(values)) || anyNA(values)) {
      fail(sprintf("Column `%s` must be a logical vector without missing values.",
                   column))
    }
  }
  
  integer_columns <- c("condition_id", "rep", "study_id", "n_mozart", "n_control")
  for (column in integer_columns) {
    lower <- if (column %in% c("n_mozart", "n_control")) 2 else 1
    values <- results[[column]]
    bad <- which(values < lower | values != floor(values))
    if (length(bad) > 0L) {
      fail(sprintf("Column `%s` must contain whole numbers >= %d; invalid row %d.",
                   column, lower, bad[1L]))
    }
  }
  for (column in c("p_value", "prob_publish_nonsig")) {
    bad <- which(results[[column]] < 0 | results[[column]] > 1)
    if (length(bad) > 0L) {
      fail(sprintf("Column `%s` must be between 0 and 1; invalid row %d.",
                   column, bad[1L]))
    }
  }
  if (any(results$tau < 0)) fail("Column `tau` must be non-negative.")
  
  # Without heterogeneity, every study has the condition's true effect.
  bad <- which(
    results$tau == 0 &
      abs(results$study_true_d - results$true_d) >
      tolerance * pmax(1, abs(results$true_d))
  )
  if (length(bad) > 0L) {
    fail(sprintf(
      "`study_true_d` does not match `true_d` when `tau = 0` at row %d.",
      bad[1L]
    ))
  }
  for (column in c("var_g", "se_g")) {
    if (any(results[[column]] <= 0)) {
      fail(sprintf("Column `%s` must be strictly positive.", column))
    }
  }
  
  # Compare using relative tolerance, with a small absolute allowance near zero.
  squared_se <- results$se_g^2
  difference <- abs(squared_se - results$var_g)
  allowance <- tolerance * pmax(1, abs(results$var_g))
  bad <- which(!is.finite(squared_se) | difference > allowance)
  if (length(bad) > 0L) {
    fail(sprintf("`se_g^2` does not match `var_g` at row %d.", bad[1L]))
  }
  
  # Hedges' g and the t statistic must have the same direction.
  bad <- which(sign(results$t) != sign(results$g))
  if (length(bad) > 0L) {
    fail(sprintf(
      "`t` and `g` do not have the same sign at row %d.",
      bad[1L]
    ))
  }
  
  # Match the simulation's two-sided, equal-variance t test.
  df <- results$n_mozart + results$n_control - 2
  bad <- which(!is.finite(df))
  if (length(bad) > 0L) {
    fail(sprintf("Non-finite t-test degrees of freedom at row %d.", bad[1L]))
  }
  
  expected_p_value <- 2 * pt(-abs(results$t), df = df)
  bad <- which(abs(results$p_value - expected_p_value) > tolerance)
  if (length(bad) > 0L) {
    fail(sprintf(
      "`p_value` does not match the two-sided t test at row %d.",
      bad[1L]
    ))
  }
  
  # A condition ID may repeat across rows, but must map to one set of settings.
  setting_columns <- c("condition_id", "true_d", "tau", "prob_publish_nonsig")
  settings <- unique(results[, setting_columns, drop = FALSE])
  if (anyDuplicated(settings$condition_id) > 0L) {
    fail("Each `condition_id` must correspond to exactly one set of condition settings.")
  }
  
  key_columns <- c("condition_id", "rep", "study_id")
  duplicated_rows <- which(duplicated(results[, key_columns, drop = FALSE]))
  if (length(duplicated_rows) > 0L) {
    fail(sprintf("Duplicate study identifier within a literature at row %d.",
                 duplicated_rows[1L]))
  }
  
  # The simulation numbers all retained studies consecutively from 1.
  # Sort IDs so validation does not depend on the input row order.
  literature_rows <- split(
    seq_len(nrow(results)),
    paste(results$condition_id, results$rep, sep = ":")
  )
  for (rows in literature_rows) {
    study_ids <- sort(results$study_id[rows])
    if (any(study_ids != seq_along(rows))) {
      i <- rows[1L]
      fail(sprintf(
        paste0("Missing or non-consecutive `study_id` values in condition_id=%s, ",
               "rep=%s. Expected IDs 1 through %s with no gaps; ",
               "provide all published and unpublished studies with their original IDs."),
        results$condition_id[i], results$rep[i], length(rows)
      ))
    }
  }
  
  # Verify the flag before using it to check the publication rules.
  expected_significant_positive <- results$p_value < alpha & results$g > 0
  bad <- which(results$significant_positive != expected_significant_positive)
  if (length(bad) > 0L) {
    fail(sprintf(
      "`significant_positive` does not match `p_value < alpha & g > 0` at row %d (alpha = %s).",
      bad[1L], format(alpha)
    ))
  }
  
  # Do not infer alpha: significant_positive means positive significance in this schema.
  # Every such study must be published under the current simulation contract.
  if (any(results$significant_positive & !results$published)) {
    fail("A study flagged `significant_positive` is unpublished, contrary to the simulation filter.")
  }
  if (any(results$prob_publish_nonsig == 1 & !results$published)) {
    fail("With `prob_publish_nonsig = 1`, all studies must be published.")
  }
  if (any(results$prob_publish_nonsig == 0 &
          results$published != results$significant_positive)) {
    fail("With `prob_publish_nonsig = 0`, publication must match positive significance.")
  }
  
  counts <- aggregate(
    results$published,
    by = results[, c("condition_id", "rep"), drop = FALSE],
    FUN = sum
  )
  names(counts)[3L] <- "n_published"
  bad <- which(counts$n_published < 1)
  if (length(bad) > 0L) {
    i <- bad[1L]
    fail(sprintf("No published studies in condition_id=%s, rep=%s.",
                 counts$condition_id[i], counts$rep[i]))
  }
  if (!is.null(expected_k)) {
    bad <- which(counts$n_published != expected_k)
    if (length(bad) > 0L) {
      i <- bad[1L]
      fail(sprintf(
        "Expected %s published studies in condition_id=%s, rep=%s; found %s.",
        expected_k, counts$condition_id[i], counts$rep[i], counts$n_published[i]
      ))
    }
  }
  
  # The simulation stops at the final required publication.
  # Use study_id rather than row order, since input rows may be shuffled.
  for (rows in literature_rows) {
    last_row <- rows[which.max(results$study_id[rows])]
    if (!results$published[last_row]) {
      fail(sprintf(
        "The last study is unpublished in condition_id=%s, rep=%s.",
        results$condition_id[last_row], results$rep[last_row]
      ))
    }
  }
  
  invisible(TRUE)
}