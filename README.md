# Publication Bias in Psychological Research

This project uses simulated Mozart-effect experiments to examine how selective publication changes meta-analytic estimates. It was developed for the Programming
in Psychological Science (PIPS) course in the Research Master Psychology at the University of Amsterdam.

Each experiment compares spatial reasoning scores between a Mozart-listening group and a control group. The analysis compares all-studies and published-only random-effects meta-analyses, applies trim-and-fill, and examines how often Egger's test detects funnel-plot asymmetry.

All participant scores and study results are simulated. The Mozart-effect literature motivates the design; this project does not establish the cause or size of the historical Mozart effect.

## Getting started

1. Clone this repository using its GitHub **Code** menu, or download and extract
the ZIP file.
2. Open the repository-root `main.Rmd` in RStudio.
3. Install the required packages in the R Console:

```r
install.packages(c("knitr", "rmarkdown", "metafor", "ggplot2"))
```

4. Click **Knit** to run the full workflow and generate `main.html`.
Alternatively, run this command in R from the repository root:

```r
rmarkdown::render("main.Rmd")
```

R is required; RStudio is recommended. HTML rendering also requires Pandoc, which is normally supplied with RStudio. No external dataset needs to be downloaded.

The document checks for missing packages and scripts. It does not install packages automatically. Keep `main.Rmd` in the repository root and preserve the folder structure below.

## Repository structure

```text
.
├── README.md
├── main.Rmd                     # Main report and workflow
├── src/
│   └── r/
│       ├── simulation.R         # Data generation and publication selection
│       ├── validation.R         # Study-level data checks
│       └── analysis.R           # Meta-analysis, summaries, and plots
├── tests/
|   ├── test_simulation.R        # Simulation checks
|   ├── test_validation.R        # Validation checks
|   └── test_analysis.R          # Analysis checks
├── figures/
|   ├── effect_estimate_comparison.png        # meta-analytic estimates across simulation conditions
|   ├── egger_test_rates.png                  # Egger's test across simulation consitions
|   └── publication_bias_heatmap.png          # Publication selection and estimation bias
```

## What the main document does

Running `main.Rmd` from top to bottom:

1. Generates simulated studies under 12 conditions.
2. Applies publication selection and retains unpublished studies.
3. Validates the study-level output.
4. Fits separate random-effects meta-analyses to all studies and published studies within each simulated literature.
5. Runs Egger's test and trim-and-fill on the published-study models.
6. Reports analysis diagnostics, valid result counts, effect estimates, and estimation errors.
7. Produces three plots: an estimate comparison, Egger rejection rates, and a published-only bias heatmap.

Each render generates fresh data using the specified seed. 
Simulation data and analysis results are used in memory; the main document does not automatically save them as separate data files. 
The HTML report contains the tables, plots, explanations, and R session information.

## Simulation settings

The default design crosses these values:

| Setting | Values |
|---|---|
| Mean population effect (`true_d`) | 0, 0.24 |
| Between-study standard deviation (`tau`) | 0, 0.3 |
| Publication probability for other results (`prob_publish_nonsig`) | 1, 0.33, 0.10 |

Positive effects with a two-sided p-value below 0.05 are always published. All other results, including significant negative effects, are published with the specified probability.

Each simulated literature contains 38 published studies, plus the unpublished studies generated before the 38th publication. Stronger selection usually requires more studies to reach this target. A literature is identified by `condition_id` and `rep`; the analysis treats each literature separately.

The YAML header in `main.Rmd` controls the run size:

```yaml
params:
  n_reps: 20
  k_published: 38
  seed: 2026
```

- `n_reps`: simulated literatures per condition.
- `k_published`: target number of published studies per literature.
- `seed`: random seed for data generation.

The default run requests 240 literatures in total. Twenty repetitions per condition are intended for development, not as evidence of adequate simulation precision. Larger runs take longer. To request 500 repetitions without editing the file:

```r
rmarkdown::render("main.Rmd", params = list(n_reps = 500))
```

Condition values can be edited in the `simulation-settings` chunk. The supplied plotting code uses labels for the default design; check that label handling supports any new values before changing the conditions.

## Dependencies

| Dependency | Purpose |
|---|---|
| R | Runs the simulation and analysis |
| `knitr`, `rmarkdown` | Execute and render `main.Rmd` |
| `metafor` | Random-effects meta-analysis, Egger's test, and trim-and-fill |
| `ggplot2` | Figures |
| Pandoc | HTML rendering |

Package dependencies are installed by `install.packages()`. The report records the R and package versions used for each render through `sessionInfo()`.

## Tests and checks

From the repository root, run these commands in R:

```r
source("tests/test_simulation.R")
source("tests/test_validation.R")
source("tests/test_analysis.R")
```

These tests are separate from rendering. During rendering, validation stops the workflow if study-level checks fail. Analysis warnings and errors are recorded in diagnostic tables; an HTML file can still be produced when some literature-level analyses fail. Check the diagnostics and valid counts before interpreting results.

For a quick integration check, render with two repetitions per condition:

```r
rmarkdown::render("main.Rmd", params = list(n_reps = 2))
```

## Interpreting the results

Mean bias is the average estimate minus the generating population effect. Mean absolute error measures distance from that effect without cancelling positive and negative errors.

Egger's test examines asymmetry, not publication bias alone. A non-significant result does not establish that publication bias is absent. Trim-and-fill is not guaranteed to recover the true effect. If methods fail, their summaries may use different sets of literatures; consult the valid counts.

The plots report point summaries without Monte Carlo uncertainty intervals. The model also assumes normal scores, equal group variances, equal group sizes, and a simplified publication rule. These limits should guide the interpretation of the simulation.

## References

Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). *Introduction to Meta-Analysis*. Wiley.

Duval, S., & Tweedie, R. (2000). Trim and fill: A simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. *Biometrics, 56*(2), 455–463.

Egger, M., Davey Smith, G., Schneider, M., & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. *BMJ, 315*, 629–634.

Franco, A., Malhotra, N., & Simonovits, G. (2014). Publication bias in the social sciences: Unlocking the file drawer. *Science, 345*, 1502–1505.

Pietschnig, J., Voracek, M., & Formann, A. K. (2010). Mozart effect–Shmozart effect: A meta-analysis. *Intelligence, 38*, 314–323.

Rauscher, F. H., Shaw, G. L., & Ky, K. N. (1993). Music and spatial task performance. *Nature, 365*, 611.

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. *Journal of Statistical Software, 36*(3), 1–48. https://doi.org/10.18637/jss.v036.i03