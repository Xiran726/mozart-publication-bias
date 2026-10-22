# Tests for src/r/validation.R using manually constructed fixtures.
# Run from the repository root: source("tests/test_validation.R")
# Uses base R only. Does not require simulation.R or analysis scripts.

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
  
  expect_error <- function(description, expr, message_fragment) {
    caught <- tryCatch({ force(expr); NULL }, error = function(e) e)
    if (!inherits(caught, "error")) {
      stop("FAILED: ", description, " (no error was raised)", call. = FALSE)
    }
    actual <- conditionMessage(caught)
    if (!grepl(message_fragment, actual, fixed = TRUE)) {
      stop("FAILED: ", description,
           "\nExpected message containing: ", message_fragment,
           "\nActual message: ", actual, call. = FALSE)
    }
    check(description, TRUE)
  }
  
  # Two conditions x two replications x two studies = eight rows.
  # All studies initially have g = t = 0 and p_value = 1.
  make_fixture <- function() {
    keys <- expand.grid(study_id = 1:2, rep = 1:2, condition_id = 1:2,
                        KEEP.OUT.ATTRS = FALSE)
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
  
  # Construct a nonzero test result with mutually consistent g, t, p and flag.
  # g is derived from t using the usual approximate Hedges correction.
  # This helper does not run a statistical simulation.
  set_test_t <- function(x, rows, t, alpha = 0.05) {
    df <- x$n_mozart[rows] + x$n_control[rows] - 2
    correction <- 1 - 3 / (4 * df - 1)
    d <- t * sqrt(1 / x$n_mozart[rows] + 1 / x$n_control[rows])
    x$t[rows] <- t
    x$g[rows] <- correction * d
    x$var_g[rows] <-
      (x$n_mozart[rows] + x$n_control[rows]) /
      (x$n_mozart[rows] * x$n_control[rows]) + x$g[rows]^2 / (2 * df)
    x$se_g[rows] <- sqrt(x$var_g[rows])
    x$p_value[rows] <- 2 * pt(-abs(t), df = df)
    x$significant_positive[rows] <-
      x$p_value[rows] < alpha & x$g[rows] > 0
    x
  }
  
  valid <- make_fixture()
  
  # 1. Valid input and return behaviour.
  check("valid data passes", validate_simulation_results(valid, expected_k = 2))
  check("NULL expected_k is accepted", validate_simulation_results(valid))
  result <- withVisible(validate_simulation_results(valid, expected_k = 2))
  check("success returns TRUE invisibly",
        identical(result$value, TRUE) && !result$visible)
  before <- valid
  validate_simulation_results(valid, expected_k = 2)
  check("validation does not modify input", identical(valid, before))
  check("reversed row order is accepted",
        validate_simulation_results(valid[8:1, ], expected_k = 2))
  extra_column <- valid
  extra_column$note <- rep("example", nrow(valid))
  check("additional columns are accepted",
        validate_simulation_results(extra_column, expected_k = 2))
  
  # 2. Validator arguments.
  for (value in list(0, -1, NA_real_, Inf, "small", numeric(0), c(1, 2))) {
    expect_error("invalid tolerance is rejected",
                 validate_simulation_results(valid, tolerance = value),
                 "`tolerance` must be a single positive finite number.")
  }
  for (value in list(-0.1, 1.1, NA_real_, Inf, "0.05", numeric(0), c(0.01, 0.05))) {
    expect_error("invalid alpha is rejected",
                 validate_simulation_results(valid, alpha = value),
                 "`alpha` must be a single finite number between 0 and 1.")
  }
  for (value in list(0, -1, 1.5, NA_real_, Inf, "2", numeric(0), c(1, 2))) {
    expect_error("invalid expected_k is rejected",
                 validate_simulation_results(valid, expected_k = value),
                 "`expected_k` must be NULL or a single positive whole number.")
  }
  check("alpha zero is accepted",
        validate_simulation_results(valid, expected_k = 2, alpha = 0))
  check("alpha one is accepted",
        validate_simulation_results(valid, expected_k = 2, alpha = 1))
  
  # 3. Structure, required columns, types, and finite values.
  expect_error("non-data-frame input is rejected",
               validate_simulation_results(list()), "non-empty data frame")
  expect_error("empty data frame is rejected",
               validate_simulation_results(valid[FALSE, ]), "non-empty data frame")
  bad <- valid
  names(bad)[2] <- names(bad)[1]
  expect_error("duplicate column names are rejected",
               validate_simulation_results(bad), "duplicate column names")
  for (column in names(valid)) {
    bad <- valid
    bad[[column]] <- NULL
    expect_error(paste("missing required column", column),
                 validate_simulation_results(bad),
                 paste0("Missing required columns: ", column, "."))
  }
  numeric_columns <- setdiff(names(valid), c("significant_positive", "published"))
  for (column in numeric_columns) {
    bad <- valid
    bad[[column]] <- as.character(bad[[column]])
    expect_error(paste("nonnumeric column", column),
                 validate_simulation_results(bad),
                 sprintf("Column `%s` must be a numeric vector.", column))
    bad <- valid
    bad[[column]] <- matrix(bad[[column]], ncol = 1)
    expect_error(paste("matrix-valued numeric column", column),
                 validate_simulation_results(bad),
                 sprintf("Column `%s` must be a numeric vector.", column))
    for (value in list(NA_real_, NaN, Inf, -Inf)) {
      bad <- valid
      bad[[column]][1] <- value
      expect_error(paste("non-finite value in", column),
                   validate_simulation_results(bad),
                   sprintf("Column `%s` contains a missing or non-finite value at row 1.",
                           column))
    }
  }
  for (column in c("significant_positive", "published")) {
    for (kind in c("numeric", "missing", "matrix")) {
      bad <- valid
      if (kind == "numeric") bad[[column]] <- as.integer(bad[[column]])
      if (kind == "missing") bad[[column]][1] <- NA
      if (kind == "matrix") bad[[column]] <- matrix(bad[[column]], ncol = 1)
      expect_error(paste(kind, "logical column", column),
                   validate_simulation_results(bad),
                   sprintf("Column `%s` must be a logical vector without missing values.",
                           column))
    }
  }
  
  # 4. Integer constraints, ranges, and variance consistency.
  for (column in c("condition_id", "rep", "study_id", "n_mozart", "n_control")) {
    lower <- if (column %in% c("n_mozart", "n_control")) 2 else 1
    for (value in c(lower - 1, lower + 0.5)) {
      bad <- valid
      bad[[column]][1] <- value
      expect_error(paste("invalid whole number in", column),
                   validate_simulation_results(bad),
                   sprintf("Column `%s` must contain whole numbers >= %d; invalid row 1.",
                           column, lower))
    }
  }
  for (column in c("p_value", "prob_publish_nonsig")) {
    for (value in c(-0.1, 1.1)) {
      bad <- valid
      bad[[column]][1] <- value
      expect_error(paste("out-of-range value in", column),
                   validate_simulation_results(bad),
                   sprintf("Column `%s` must be between 0 and 1; invalid row 1.", column))
    }
  }
  bad <- valid
  bad$tau[1] <- -0.1
  expect_error("negative tau is rejected", validate_simulation_results(bad),
               "Column `tau` must be non-negative.")
  for (column in c("var_g", "se_g")) {
    for (value in c(0, -0.1)) {
      bad <- valid
      bad[[column]][1] <- value
      expect_error(paste("nonpositive", column), validate_simulation_results(bad),
                   sprintf("Column `%s` must be strictly positive.", column))
    }
  }
  bad <- valid
  bad$se_g[1] <- 0.9
  expect_error("inconsistent squared SE is rejected", validate_simulation_results(bad),
               "`se_g^2` does not match `var_g` at row 1.")
  bad <- valid
  bad$se_g[1] <- .Machine$double.xmax
  expect_error("overflow in squared SE is rejected", validate_simulation_results(bad),
               "`se_g^2` does not match `var_g` at row 1.")
  near <- valid
  near$se_g[1] <- sqrt(near$var_g[1] + 1e-10)
  check("SE discrepancy within tolerance is accepted",
        validate_simulation_results(near, expected_k = 2))
  
  # 5. Review 2: no heterogeneity means equal true effects.
  bad <- valid
  bad$study_true_d[1] <- 0.5
  expect_error("different true effect with tau zero is rejected",
               validate_simulation_results(bad),
               "`study_true_d` does not match `true_d` when `tau = 0` at row 1.")
  near <- valid
  near$study_true_d[1] <- 1e-10
  check("true-effect discrepancy within tolerance is accepted",
        validate_simulation_results(near, expected_k = 2))
  heterogeneous <- valid
  heterogeneous$tau <- rep(0.3, nrow(heterogeneous))
  heterogeneous$study_true_d[1] <- 0.5
  check("different true effect with positive tau is accepted",
        validate_simulation_results(heterogeneous, expected_k = 2))
  
  # 6. Review 3: t direction and two-sided p-value.
  positive <- set_test_t(valid, rows = 1, t = 3)
  check("consistent positive t result is accepted",
        validate_simulation_results(positive, expected_k = 2))
  negative <- set_test_t(valid, rows = 1, t = -3)
  check("consistent negative t result is accepted",
        validate_simulation_results(negative, expected_k = 2))
  bad <- positive
  bad$t[1] <- -bad$t[1]
  expect_error("flipped t sign is rejected", validate_simulation_results(bad),
               "`t` and `g` do not have the same sign at row 1.")
  bad <- positive
  bad$g[1] <- -bad$g[1]
  expect_error("flipped g sign is rejected", validate_simulation_results(bad),
               "`t` and `g` do not have the same sign at row 1.")
  for (column in c("g", "t")) {
    bad <- positive
    bad[[column]][1] <- 0
    expect_error(paste("zero", column, "with nonzero counterpart is rejected"),
                 validate_simulation_results(bad),
                 "`t` and `g` do not have the same sign at row 1.")
  }
  bad <- positive
  bad$p_value[1] <- 0.02
  expect_error("p-value inconsistent with t is rejected", validate_simulation_results(bad),
               "`p_value` does not match the two-sided t test at row 1.")
  near <- positive
  near$p_value[1] <- near$p_value[1] + 1e-10
  check("p-value discrepancy within tolerance is accepted",
        validate_simulation_results(near, expected_k = 2))
  unequal <- valid
  unequal$n_mozart[1] <- 21
  unequal <- set_test_t(unequal, rows = 1, t = 3)
  check("unequal group sizes use the correct degrees of freedom",
        validate_simulation_results(unequal, expected_k = 2))
  bad <- valid
  bad$n_mozart[1] <- .Machine$double.xmax
  bad$n_control[1] <- .Machine$double.xmax
  expect_error("overflow in degrees of freedom is rejected",
               validate_simulation_results(bad),
               "Non-finite t-test degrees of freedom at row 1.")
  
  # 7. Condition settings and study identifiers.
  for (column in c("true_d", "tau", "prob_publish_nonsig")) {
    bad <- valid
    bad[[column]][1] <- if (column == "prob_publish_nonsig") 0.33 else 0.5
    if (column == "true_d") bad$study_true_d[1] <- bad$true_d[1]
    expect_error(paste("conflicting condition setting", column),
                 validate_simulation_results(bad),
                 "Each `condition_id` must correspond to exactly one set of condition settings.")
  }
  bad <- rbind(valid, valid[1, ])
  expect_error("duplicate study key is rejected", validate_simulation_results(bad),
               "Duplicate study identifier within a literature at row 9.")
  bad <- valid
  bad$study_id[2] <- 3
  expect_error("gap in study IDs is rejected", validate_simulation_results(bad),
               "Missing or non-consecutive `study_id` values in condition_id=1, rep=1.")
  bad <- valid
  bad$study_id <- bad$study_id + 1
  expect_error("study IDs must start at one", validate_simulation_results(bad),
               "Expected IDs 1 through 2 with no gaps")
  
  # 8. Review 1: check the significance flag itself.
  bad <- valid
  bad$significant_positive[1] <- TRUE
  expect_error("incorrect TRUE significance flag is rejected",
               validate_simulation_results(bad),
               "`significant_positive` does not match `p_value < alpha & g > 0` at row 1")
  bad <- positive
  bad$significant_positive[1] <- FALSE
  expect_error("incorrect FALSE significance flag is rejected",
               validate_simulation_results(bad),
               "`significant_positive` does not match `p_value < alpha & g > 0` at row 1")
  bad <- negative
  bad$significant_positive[1] <- TRUE
  expect_error("negative significant effect cannot be flagged positive",
               validate_simulation_results(bad),
               "`significant_positive` does not match")
  check("zero effect, zero t and p one are accepted",
        validate_simulation_results(valid, expected_k = 2))
  
  boundary <- set_test_t(valid, rows = 1,
                         t = qt(0.05 / 2, df = 38, lower.tail = FALSE))
  # Preserve the exact threshold under test; round-trip t-to-p error is tiny.
  boundary$p_value[1] <- 0.05
  boundary$significant_positive[1] <- FALSE
  check("p equal to alpha is not significant",
        validate_simulation_results(boundary, expected_k = 2))
  bad <- boundary
  bad$significant_positive[1] <- TRUE
  expect_error("TRUE flag at exact alpha boundary is rejected",
               validate_simulation_results(bad),
               "`significant_positive` does not match")
  custom <- set_test_t(valid, rows = 1,
                       t = qt(0.08 / 2, df = 38, lower.tail = FALSE), alpha = 0.10)
  custom$p_value[1] <- 0.08
  check("custom alpha is respected",
        validate_simulation_results(custom, expected_k = 2, alpha = 0.10))
  expect_error("custom-alpha flag fails at default alpha",
               validate_simulation_results(custom, expected_k = 2),
               "`significant_positive` does not match")
  
  # 9. Publication rules and counts.
  # Every literature has an unpublished ID 1 and a published ID 2.
  selected <- valid
  selected$prob_publish_nonsig <- rep(0.33, nrow(selected))
  selected$published <- selected$study_id == 2
  check("retained unpublished rows are accepted",
        validate_simulation_results(selected, expected_k = 1))
  check("NULL expected_k permits another publication count",
        validate_simulation_results(selected))
  bad <- set_test_t(selected, rows = 1, t = 3)
  expect_error("unpublished positive-significant study is rejected",
               validate_simulation_results(bad),
               "A study flagged `significant_positive` is unpublished, contrary to the simulation filter.")
  bad <- valid
  bad$published[1] <- FALSE
  expect_error("probability one requires every study to be published",
               validate_simulation_results(bad),
               "With `prob_publish_nonsig = 1`, all studies must be published.")
  only_positive <- selected
  only_positive$prob_publish_nonsig <- rep(0, nrow(only_positive))
  only_positive <- set_test_t(only_positive, which(only_positive$published), t = 3)
  check("zero-probability publication policy is accepted",
        validate_simulation_results(only_positive, expected_k = 1))
  bad <- only_positive
  bad$published[1] <- TRUE
  expect_error("nonsignificant publication with probability zero is rejected",
               validate_simulation_results(bad),
               "With `prob_publish_nonsig = 0`, publication must match positive significance.")
  bad <- selected
  bad$published[bad$condition_id == 1 & bad$rep == 1] <- FALSE
  expect_error("literature with no publications is identified",
               validate_simulation_results(bad),
               "No published studies in condition_id=1, rep=1.")
  expect_error("too few publications relative to expected_k are identified",
               validate_simulation_results(valid, expected_k = 3),
               "Expected 3 published studies in condition_id=1, rep=1; found 2.")
  expect_error("too many publications relative to expected_k are identified",
               validate_simulation_results(valid, expected_k = 1),
               "Expected 1 published studies in condition_id=1, rep=1; found 2.")
  
  # 10. Review 4: the last study must be published, independently of row order.
  extra <- selected[1, ]
  extra$study_id <- 3
  bad <- rbind(selected, extra)
  expect_error("extra unpublished study after stopping point is rejected",
               validate_simulation_results(bad, expected_k = 1),
               "The last study is unpublished in condition_id=1, rep=1.")
  expect_error("stopping rule also applies with NULL expected_k",
               validate_simulation_results(bad),
               "The last study is unpublished in condition_id=1, rep=1.")
  check("last-study check uses IDs, not row order",
        validate_simulation_results(selected[8:1, ], expected_k = 1))
  
  # 11. Regression: deleting 106 unpublished rows leaves 38 publications,
  # but the original study IDs reveal missing studies.
  fixture <- valid[rep(1L, 144), , drop = FALSE]
  rownames(fixture) <- NULL
  fixture$condition_id <- rep(10, 144)
  fixture$rep <- rep(1, 144)
  fixture$study_id <- seq_len(144)
  fixture$prob_publish_nonsig <- rep(0.10, 144)
  fixture$published <- c(rep(FALSE, 106), rep(TRUE, 38))
  check("complete 144-study file-drawer fixture is accepted",
        validate_simulation_results(fixture, expected_k = 38))
  check("reordered file-drawer fixture is accepted",
        validate_simulation_results(fixture[144:1, ], expected_k = 38))
  published_only <- fixture[fixture$published, , drop = FALSE]
  check("regression fixture retains exactly 38 publications",
        nrow(published_only) == 38L && sum(published_only$published) == 38)
  expect_error("removing 106 unpublished studies fails despite 38 publications",
               validate_simulation_results(published_only, expected_k = 38),
               "Missing or non-consecutive `study_id` values in condition_id=10, rep=1.")
  
  cat("\nAll", passed, "validation checks passed.\n")
})
