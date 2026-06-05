# Scripts

This directory contains all code for data preparation, model fitting, diagnostics, and results generation. Scripts are the authoritative record of every analytical decision — if something was done to the data or model, it must be in a script.

## Script Execution Order

Scripts should be numbered and run in sequence. Document that order here so a reader can reproduce the full analysis without guessing:

| Order | Script | Description |
|-------|--------|-------------|
| 01 | `01_data_prep.R` | Load raw data, clean, filter, and write to `data/processed/` |
| 02 | `02_fit_model.R` | Fit the statistical or process model; save model object to `results/model-fits/` |
| 03 | `03_diagnostics.R` | Run convergence checks, posterior predictive checks, or cross-validation |
| 04 | `04_results.R` | Generate all figures and tables; write to `results/figures/` and `results/tables/` |

A top-level `run_all.R` or `Makefile` that runs scripts in order is strongly encouraged so the full analysis can be reproduced with a single command. Note that we intend to transition this to a `targets` workflow.

## Best Practices

### Reproducibility

- **Set a random seed** at the top of every script that involves any stochastic process (MCMC sampling, bootstrapping, simulation, train/test splits). Use a consistent seed across the project.
  ```r
  set.seed(42)  # R
  ```
  ```python
  import random, numpy as np
  random.seed(42)
  np.random.seed(42)
  ```
- **Use relative paths** rooted at the project directory. Never hardcode absolute paths. In R, use the `here` package:
  ```r
  library(here)
  data <- read.csv(here("data", "processed", "model-data.csv"))
  ```
- **Run scripts from the project root**, not from the `scripts/` subdirectory. This ensures relative paths resolve correctly.
- **Save a record of the software environment** at the end of a complete run:
  ```r
  sink(here("scripts", "session_info.txt"))
  sessionInfo()
  sink()
  ```
  ```python
  # run: pip freeze > scripts/requirements.txt
  ```

### Code Organization and Style

- One script, one purpose. Each script should do one logical thing (data prep, model fitting, results). Avoid monolithic scripts that do everything.
- Scripts should be self-contained: any object a script needs that was created by a prior script should be loaded from a saved file, not passed through a global environment. This makes it easy to re-run individual steps without re-running the entire pipeline.
- Use consistent naming conventions throughout: `snake_case` for variables and functions is conventional in R and Python.
- Keep lines to a reasonable length (≤ 100 characters) for readability in diff views and printed manuscripts.
- Delete commented-out code before publishing. Use git history to recover old approaches if needed.

### Documentation Within Scripts

Each script should begin with a header block that describes its purpose, inputs, and outputs:

```r
# =============================================================================
# Script: 02_fit_model.R
# Purpose: Fit the Bayesian hierarchical model to processed survey data
# Inputs:  data/processed/model-data.csv
# Outputs: results/model-fits/model_fit.rds
# Author:  [Name]
# Date:    [YYYY-MM-DD]
# =============================================================================
```

- Comment *why*, not *what*. Code explains what it does; comments explain why a choice was made.
  ```r
  # Exclude years prior to 1995 — monitoring protocol changed, data are not comparable
  df <- df[df$year >= 1995, ]
  ```
- Document any hard-coded constants (prior parameters, cutoff values, thresholds) with a justification or citation.

### Functions and Modularity

- If a piece of logic is used more than once, put it in a function. Store shared utility functions in a dedicated file (e.g., `scripts/utils.R` or `scripts/helpers.py`) and source/import it at the top of each script that needs it.
- Write functions that accept arguments rather than reading from the global environment. This makes functions testable and reusable.
- Consider writing unit tests for key functions, especially any custom likelihood, data transformation, or summary calculation.

### Version Control

- Commit scripts to git frequently and with informative commit messages:
  ```
  Add prior sensitivity analysis for detection probability
  ```
  Not:
  ```
  updates
  ```
- Use a `.gitignore` to exclude system files, large outputs, and environment directories (`.DS_Store`, `__pycache__/`, `.Rproj.user/`, `renv/library/`).
- Do not commit rendered outputs (HTML notebooks, PDFs) to git unless they are the final manuscript supplement. These are large binary files and inflate repository size quickly.

### Dependency Management

- **R:** Use [`renv`](https://rstudio.github.io/renv/) to snapshot and restore the package library. Commit `renv.lock` to git.
  ```r
  renv::init()    # initialize
  renv::snapshot() # update lockfile after adding packages
  renv::restore()  # recreate library from lockfile
  ```
- **Python:** Pin dependencies in `requirements.txt` (for pip) or `environment.yml` (for conda). Include both direct and transitive dependencies.
- Document any system-level dependencies (e.g., JAGS, Stan/CmdStan, GDAL) and the versions used, since these are not captured by language-level package managers.

### Long-Running Jobs

- For models that take hours to run, save intermediate checkpoints so the job can be resumed without starting over.
- Save fitted model objects to `results/model-fits/` immediately after fitting. Downstream scripts should load the saved object rather than re-fitting.
- For MCMC models, save the full posterior samples, not just posterior summaries — summaries can always be recomputed, but re-running a long chain is costly.
- Consider adding a `--dry-run` flag or short test mode (small data subset, fewer iterations) that lets a reader verify the pipeline runs end-to-end quickly before committing to a full run.
