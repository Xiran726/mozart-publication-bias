# Tests for src/r/validation.R using small, manually constructed fixtures.
# Run from the repository root: source("tests/test_validation.R")
# Uses base R only; simulation.R and analysis scripts are not required.

local({
  source(file.path("src", "r", "validation.R"), local = TRUE)
  passed <- 0L
  
  check <- function(description, condition) {
    if (!isTRUE(condition)) {
      stop("FAILED: ", description, call. = FALSE)
    }
    passed <<- passed + 1L
    cat("PASS:", description, "\n")
  }
  
  # Verify the error message, not just the occurrence of any error.
  expect_error <- function(description, expr, message_fragment) {
    caught <- tryCatch(
      {
        force(expr)
        NULL
      },
      error = function(e) e
    )
    if (!inherits(caught, "error")) {
      stop("FAILED: ", description, " (no error was raised)", call. = FALSE)
    }
    actual_message <- conditionMessage(caught)
    if (!grepl(message_fragment, actual_message, fixed = TRUE)) {
      stop(
        "FAILED: ", description,
        "\nExpected message containing: ", message_fragment,
        "\nActual message: ", actual_message,
        call. = FALSE
      )
    }
    check(description, TRUE)
  }
  
  # A fixture is controlled example data, not a simulated research result.
  # Two conditions x two replications x two studies = eight rows.
  make_fixture <- function() {
    keys <- expand.grid(
      study_id = 1:2,
      rep = 1:2,
      condition_id = 1:2,
      KEEP.OUT.ATTRS = FALSE
    )
    n <- nrow(keys)
    data.frame(
      condition_id = keys$condition_id,
      true_d = ifelse(keys$condition_id == 1, 0, 0.24),
      tau = rep(0, n),
      prob_publish_nonsig = rep(1, n),
      rep = keys$rep,
      study_id = keys$study_id,
      study_true_d = ifelse(keys$condition_id == 1, 0, 0.24),
      n_mozart = rep(20, n),
      n_control = rep(20, n),
      g = rep(0, n),
      var_g = rep(0.1, n),
      se_g = rep(sqrt(0.1), n),
      t = rep(0, n),
      p_value = rep(1, n),
      significant_positive = rep(FALSE, n),
      published = rep(TRUE, n)
    )
  }
  
  valid <- make_fixture()
  
  # Valid input, non-default settings, and return behaviour.
  check("valid data passes with two published studies per literature",
        validate_simulation_results(valid, expected_k = 2))
  check("expected_k = NULL skips the exact-count requirement",
        validate_simulation_results(valid))
  returned <- withVisible(validate_simulation_results(valid, expected_k = 2))
  check("success returns TRUE invisibly",
        identical(returned$value, TRUE) && !returned$visible)
  before <- valid
  validate_simulation_results(valid, expected_k = 2)
  check("validation does not modify its input", identical(valid, before))
  shuffled <- valid[rev(seq_len(nrow(valid))), , drop = FALSE]
  check("row order does not affect validity",
        validate_simulation_results(shuffled, expected_k = 2))
  
  # Input structure.
  expect_error("non-data-frame input is rejected",
               validate_simulation_results(list()),
               "non-empty data frame")
  expect_error("empty data frame is rejected",
               validate_simulation_results(valid[FALSE, , drop = FALSE]),
               "non-empty data frame")
  bad <- valid
  bad$var_g <- NULL
  expect_error("missing variance column is identified",
               validate_simulation_results(bad),
               "Missing required columns: var_g")
  bad <- valid
  names(bad)[2] <- names(bad)[1]
  expect_error("duplicate column names are rejected",
               validate_simulation_results(bad),
               "duplicate column names")
  bad <- valid
  bad$g <- as.character(bad$g)
  expect_error("text in a numeric column is rejected",
               validate_simulation_results(bad),
               "Column `g` must be a numeric vector")
  
  # Missing, non-finite, and invalid numerical values.
  bad <- valid
  bad$g[1] <- NA_real_
  expect_error("missing effect size is identified",
               validate_simulation_results(bad),
               "Column `g` contains a missing or non-finite value at row 1")
  bad <- valid
  bad$g[1] <- Inf
  expect_error("infinite effect size is identified",
               validate_simulation_results(bad),
               "Column `g` contains a missing or non-finite value at row 1")
  bad <- valid
  bad$n_mozart[1] <- 1
  expect_error("group size below two is rejected",
               validate_simulation_results(bad),
               "Column `n_mozart` must contain whole numbers >= 2")
  bad <- valid
  bad$study_id[1] <- 1.5
  expect_error("fractional study identifier is rejected",
               validate_simulation_results(bad),
               "Column `study_id` must contain whole numbers >= 1")
  bad <- valid
  bad$p_value[1] <- 1.1
  expect_error("p-value above one is rejected",
               validate_simulation_results(bad),
               "Column `p_value` must be between 0 and 1")
  bad <- valid
  bad$prob_publish_nonsig[1] <- -0.1
  expect_error("negative publication probability is rejected",
               validate_simulation_results(bad),
               "Column `prob_publish_nonsig` must be between 0 and 1")
  bad <- valid
  bad$tau[1] <- -0.1
  expect_error("negative heterogeneity is rejected",
               validate_simulation_results(bad),
               "Column `tau` must be non-negative")
  bad <- valid
  bad$var_g[1] <- 0
  expect_error("zero variance is rejected",
               validate_simulation_results(bad),
               "Column `var_g` must be strictly positive")
  bad <- valid
  bad$var_g[1] <- -0.1
  expect_error("negative variance is rejected",
               validate_simulation_results(bad),
               "Column `var_g` must be strictly positive")
  bad <- valid
  bad$se_g[1] <- 0
  expect_error("zero standard error is rejected",
               validate_simulation_results(bad),
               "Column `se_g` must be strictly positive")
  bad <- valid
  bad$se_g[1] <- 0.9
  expect_error("standard error inconsistent with variance is rejected",
               validate_simulation_results(bad),
               "`se_g^2` does not match `var_g` at row 1")
  
  # Logical flags and grouping identifiers.
  bad <- valid
  bad$published <- as.integer(bad$published)
  expect_error("numeric publication flags are rejected",
               validate_simulation_results(bad),
               "Column `published` must be a logical vector without missing values")
  bad <- valid
  bad$published[1] <- NA
  expect_error("missing publication flags are rejected",
               validate_simulation_results(bad),
               "Column `published` must be a logical vector without missing values")
  bad <- valid
  bad$significant_positive[1] <- NA
  expect_error("missing significance flags are rejected",
               validate_simulation_results(bad),
               "Column `significant_positive` must be a logical vector without missing values")
  bad <- valid
  bad$true_d[1] <- 0.5
  expect_error("conflicting settings under the same condition ID are rejected",
               validate_simulation_results(bad),
               "exactly one set of condition settings")
  bad <- rbind(valid, valid[1, , drop = FALSE])
  expect_error("duplicate study within a literature is rejected",
               validate_simulation_results(bad),
               "Duplicate study identifier within a literature")
  
  # Publication-filter contracts. Use an intermediate probability where needed
  # so a different filter rule does not mask the intended test.
  selected <- valid
  selected$prob_publish_nonsig <- rep(0.33, nrow(selected))
  selected$published <- selected$study_id == 1
  check("retained file-drawer rows are accepted",
        validate_simulation_results(selected, expected_k = 1))
  bad <- selected
  bad$significant_positive[2] <- TRUE
  bad$g[2] <- 0.8
  bad$t[2] <- 3
  bad$p_value[2] <- 0.01
  expect_error("unpublished positive-significant study is rejected",
               validate_simulation_results(bad),
               "A study flagged `significant_positive` is unpublished")
  bad <- valid
  bad$published[1] <- FALSE
  expect_error("unpublished study under the no-bias filter is rejected",
               validate_simulation_results(bad),
               "With `prob_publish_nonsig = 1`, all studies must be published")
  only_positive <- selected
  only_positive$prob_publish_nonsig <- rep(0, nrow(only_positive))
  only_positive$significant_positive <- only_positive$published
  only_positive$g[only_positive$significant_positive] <- 0.8
  only_positive$t[only_positive$significant_positive] <- 3
  only_positive$p_value[only_positive$significant_positive] <- 0.01
  check("only-positive-significant publication policy is accepted",
        validate_simulation_results(only_positive, expected_k = 1))
  bad <- only_positive
  bad$published[2] <- TRUE
  expect_error("publication inconsistent with a zero-probability filter is rejected",
               validate_simulation_results(bad),
               "With `prob_publish_nonsig = 0`, publication must match positive significance")
  bad <- selected
  bad$published[bad$condition_id == 1 & bad$rep == 1] <- FALSE
  expect_error("literature with no published studies is identified",
               validate_simulation_results(bad),
               "No published studies in condition_id=1, rep=1")
  expect_error("wrong expected publication count identifies the literature",
               validate_simulation_results(valid, expected_k = 3),
               "Expected 3 published studies in condition_id=1, rep=1; found 2")
  
  # Validator arguments.
  expect_error("zero expected_k is rejected",
               validate_simulation_results(valid, expected_k = 0),
               "`expected_k` must be NULL or a single positive whole number")
  expect_error("fractional expected_k is rejected",
               validate_simulation_results(valid, expected_k = 1.5),
               "`expected_k` must be NULL or a single positive whole number")
  expect_error("zero tolerance is rejected",
               validate_simulation_results(valid, tolerance = 0),
               "`tolerance` must be a single positive finite number")
  expect_error("missing tolerance is rejected",
               validate_simulation_results(valid, tolerance = NA_real_),
               "`tolerance` must be a single positive finite number")
  
  cat("\nAll", passed, "validation checks passed.\n")
})

local({
  source(file.path("src", "r", "validation.R"), local = TRUE)
  
  # Minimal fixture for the reported missing-unpublished-studies case.
  published <- c(rep(FALSE, 106), rep(TRUE, 38))
  fixture <- data.frame(
    condition_id = rep(10, 144),
    true_d = rep(0, 144),
    tau = rep(0, 144),
    prob_publish_nonsig = rep(0.10, 144),
    rep = rep(1, 144),
    study_id = seq_len(144),
    study_true_d = rep(0, 144),
    n_mozart = rep(20, 144),
    n_control = rep(20, 144),
    g = rep(0, 144),
    var_g = rep(0.04, 144),
    se_g = rep(0.2, 144),
    t = rep(0, 144),
    p_value = rep(1, 144),
    significant_positive = rep(FALSE, 144),
    published = published
  )
  
  stopifnot(isTRUE(
    validate_simulation_results(fixture, expected_k = 38)
  ))
  
  # Row order should not affect validation.
  stopifnot(isTRUE(
    validate_simulation_results(fixture[144:1, ], expected_k = 38)
  ))
  
  published_only <- fixture[fixture$published, , drop = FALSE]
  stopifnot(sum(published_only$published) == 38)
  
  error <- tryCatch(
    {
      validate_simulation_results(published_only, expected_k = 38)
      NULL
    },
    error = function(e) e
  )
  stopifnot(
    inherits(error, "error"),
    grepl("study_id", conditionMessage(error), fixed = TRUE),
    grepl("condition_id=10, rep=1", conditionMessage(error), fixed = TRUE)
  )
  
  cat("PASS: Complete and reordered data pass validation.\n")
  cat("PASS: Removing 106 unpublished studies fails despite retaining 38 publications.\n")
})