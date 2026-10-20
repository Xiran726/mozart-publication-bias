# =============================================================================
# test_simulation.R
#
# Checks that the functions in src/r/simulation.R do what they claim.
# Run from the project folder with:  source("tests/test_simulation.R")
# Every check prints "PASS"; a failing check stops with an error message.
# Takes about a minute. No packages needed.
# =============================================================================

source(file.path("src", "r", "simulation.R"))
set.seed(123)

passed <- 0
check <- function(description, condition) {
  if (!isTRUE(condition)) stop("FAILED: ", description, call. = FALSE)
  passed <<- passed + 1
  cat("PASS:", description, "\n")
}
gives_error <- function(expr) {
  inherits(tryCatch(expr, error = function(e) e), "error")
}


# --- 1. One experiment ------------------------------------------------------

data <- generate_study_data(n_per_group = 25, true_d = 0.5)
check("generate_study_data() creates 25 participants per group",
      nrow(data) == 50 && all(table(data$condition) == 25))

result <- summarise_study(data)
reference <- t.test(score ~ condition, data = data, var.equal = TRUE)
# t.test() computes control - mozart, so its t has the opposite sign
check("t and p match R's own t.test()",
      abs(result$t + reference$statistic) < 1e-10 &&
        abs(result$p_value - reference$p.value) < 1e-10)

# Without a real effect, about 5% of studies should be significant
null_p <- replicate(4000, simulate_mozart_study(n_per_group = 20,
                                                true_d = 0)$p_value)
check("false-positive rate is about 5% when there is no effect",
      abs(mean(null_p < 0.05) - 0.05) < 0.015)

# Hedges' g should be unbiased: on average it equals the true effect
g_values <- replicate(4000, simulate_mozart_study(n_per_group = 10,
                                                  true_d = 0.5)$g)
check("Hedges' g is on average equal to the true effect (0.5)",
      abs(mean(g_values) - 0.5) < 0.03)


# --- 2. Literature ------------------------------------------------------------

lit <- simulate_literature(n_studies = 200, true_d = 0.3)
check("simulate_literature() returns one row per study",
      nrow(lit) == 200)
check("sample sizes come from the real Mozart literature",
      all((lit$n_mozart + lit$n_control) %in%
            (2 * round(MOZART_TOTAL_N / 2))))
check("with tau = 0, every study has the same true effect",
      all(lit$study_true_d == 0.3))
lit_tau <- simulate_literature(n_studies = 2000, true_d = 0.3, tau = 0.2)
check("with tau = 0.2, true effects vary between studies with SD ~0.2",
      abs(sd(lit_tau$study_true_d) - 0.2) < 0.02)


# --- 3. Journal filter --------------------------------------------------------

check("prob_publish_nonsig = 1 publishes every study",
      all(apply_publication_bias(lit, prob_publish_nonsig = 1)$published))
only_sig <- apply_publication_bias(lit, prob_publish_nonsig = 0)
check("prob_publish_nonsig = 0 publishes only significant positive studies",
      identical(only_sig$published, only_sig$significant))
check("significant means p < .05 AND Mozart group higher",
      all(only_sig$p_value[only_sig$significant] < 0.05) &&
        all(only_sig$g[only_sig$significant] > 0))


# --- 4. Fixed number of published studies -------------------------------------

pub_lit <- simulate_published_literature(k_published = 38, true_d = 0,
                                         prob_publish_nonsig = 0.10)
check("exactly 38 studies are published",
      sum(pub_lit$published) == 38)
check("the last study run is the 38th published one (nothing extra kept)",
      pub_lit$published[nrow(pub_lit)])
check("file-drawer studies are kept (more than 38 studies were run)",
      nrow(pub_lit) > 38)
check("no publication bias: exactly 38 studies are run",
      nrow(simulate_published_literature(38, prob_publish_nonsig = 1)) == 38)


# --- 5./6. Full simulation ------------------------------------------------------

small_conditions <- make_conditions(true_d = c(0, 0.24),
                                    prob_publish_nonsig = c(1, 0.10),
                                    tau = 0)
run_a <- run_simulation(small_conditions, n_reps = 3, seed = 1,
                        verbose = FALSE)
run_b <- run_simulation(small_conditions, n_reps = 3, seed = 1,
                        verbose = FALSE)
check("make_conditions() crosses all factors (2 x 2 x 1 = 4 conditions)",
      nrow(small_conditions) == 4)
check("every condition and replication is present",
      all(table(unique(run_a[, c("condition_id", "rep")])$condition_id) == 3))
check("every simulated literature has exactly 38 published studies",
      all(tapply(run_a$published,
                 paste(run_a$condition_id, run_a$rep), sum) == 38))
check("the same seed gives exactly the same results (reproducible)",
      identical(run_a[names(run_a)], run_b[names(run_b)]))


# --- 7. Publication bias inflates published effects ----------------------------

summary_a <- describe_simulation(
  run_simulation(small_conditions, n_reps = 30, seed = 2, verbose = FALSE))
no_bias <- summary_a[summary_a$prob_publish_nonsig == 1, ]
strong  <- summary_a[summary_a$prob_publish_nonsig == 0.10, ]
check("without bias, published effects are close to the truth",
      all(abs(no_bias$mean_g_published - no_bias$true_d) < 0.05))
check("with strong bias, published effects are inflated",
      all(strong$mean_g_published > strong$true_d + 0.05))


# --- 8. Informative errors --------------------------------------------------

check("negative sample size gives an error",
      gives_error(generate_study_data(n_per_group = -5)))
check("text instead of a number gives an error",
      gives_error(simulate_literature(n_studies = 10, true_d = "big")))
check("probability above 1 gives an error",
      gives_error(apply_publication_bias(lit, prob_publish_nonsig = 1.5)))
check("negative heterogeneity gives an error",
      gives_error(simulate_literature(n_studies = 10, tau = -0.1)))
check("impossible publication target stops instead of running forever",
      gives_error(simulate_published_literature(
        k_published = 38, true_d = -1, prob_publish_nonsig = 0,
        max_studies = 200)))

cat("\nAll", passed, "checks passed.\n")