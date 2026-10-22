# Tests for the statistical analysis -----------------------------

# These tests use small example datasets to check that our
# analysis functions behave as expected.

source("src/r/analysis.R")

# 1. Test the meta-analysis function ----------------------------

example_studies <- data.frame(
  g = c(-0.10, 0.05, 0.20, 0.12, 0.35,
        0.08, 0.25, -0.05, 0.18, 0.30),
  var_g = c(0.05, 0.03, 0.06, 0.07, 0.04,
            0.08, 0.05, 0.09, 0.035, 0.065)
)

model <- fit_meta_analysis(example_studies)

stopifnot(inherits(model, "rma"))
stopifnot(is.finite(as.numeric(coef(model))))

# Check that invalid sampling variances are rejected
invalid_studies <- data.frame(
  g = c(0.10, 0.20),
  var_g = c(0.05, -0.01)
)

invalid_result <- try(
  fit_meta_analysis(invalid_studies),
  silent = TRUE
)

stopifnot(inherits(invalid_result, "try-error"))

cat("Meta-analysis tests passed!\n")


# 2. Test the condition summaries -------------------------------

# Two example conditions, with two repetitions each
example_estimates <- data.frame(
  condition_id = c(1, 1, 2, 2),
  rep = c(1, 2, 1, 2),
  true_d = c(0, 0, 0.24, 0.24),
  tau = c(0, 0, 0.3, 0.3),
  prob_publish_nonsig = c(1, 1, 0.10, 0.10),
  estimate_published = c(0.10, -0.10, 0.40, 0.50),
  estimate_all = c(0.02, -0.02, 0.25, 0.23),
  estimate_trimfill = c(0.05, -0.05, 0.35, NA_real_),
  egger_detected = c(FALSE, FALSE, TRUE, FALSE)
)

# Calculate bias relative to the known true effect
example_estimates$bias_published <-
  example_estimates$estimate_published - example_estimates$true_d

example_estimates$bias_all <-
  example_estimates$estimate_all - example_estimates$true_d

example_estimates$bias_trimfill <-
  example_estimates$estimate_trimfill - example_estimates$true_d

example_summary <- summarise_conditions(example_estimates)

# Check the number of conditions and repetitions
stopifnot(nrow(example_summary) == 2)
stopifnot(all(example_summary$n_literatures == 2))

# Check that missing trim-and-fill results are handled correctly
stopifnot(sum(example_summary$n_trimfill_valid) == 3)

# In condition 2, one out of two Egger tests is significant
stopifnot(
  example_summary$egger_significant_rate[2] == 0.5
)

# Check the expected mean publication bias in condition 2
stopifnot(
  isTRUE(all.equal(example_summary$mean_bias_published[2], 0.21))
)

cat("Condition summary tests passed!\n")


# 3. Test the three plotting functions --------------------------

heatmap_test <- plot_bias_heatmap(example_summary)

comparison_test <- plot_estimate_comparison(example_summary)

egger_test_plot <- plot_egger_rates(example_summary)

stopifnot(inherits(heatmap_test, "ggplot"))
stopifnot(inherits(comparison_test, "ggplot"))
stopifnot(inherits(egger_test_plot, "ggplot"))

cat("Visualisation tests passed!\n")


# 4. Test one complete simulated literature ---------------------

# Use actual simulation data when the file is available.
# This checks the connection between simulation and analysis
# without rerunning all 6,000 literatures.

if (file.exists("data/simulated_literatures.rds")) {
  
  results <- readRDS("data/simulated_literatures.rds")
  
  one_literature <- subset(
    results,
    condition_id == 10 & rep == 1
  )
  
  one_result <- analyse_one_literature(one_literature)
  
  stopifnot(is.data.frame(one_result))
  stopifnot(nrow(one_result) == 1)
  stopifnot(one_result$n_published == 38)
  stopifnot(one_result$n_all == nrow(one_literature))
  stopifnot(one_result$n_all > one_result$n_published)
  stopifnot(is.finite(one_result$estimate_published))
  stopifnot(is.finite(one_result$estimate_all))
  
  cat("Single-literature integration test passed!\n")
  
} else {
  message("Integration test skipped: simulation data not available.")
}

cat("All completed analysis tests passed!\n")

