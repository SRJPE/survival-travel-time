# Joint CJS + CWT simulation (draft, 30 Sep 2026)

## Question
Can the acoustic telemetry (CJS) and Knights Landing (KL) CWT data be fit as **one**
survival and travel-time process with two observation models, instead of separate
models linked by soft-shared slopes (Flora's 09)? And does that fix the CJS vs CWT flow
disagreement if KL trap efficiency falls with flow?

## Files
| File | What it does |
|---|---|
| `joint_cjs_cwt.stan` | One model, switches select the four fits below |
| `sim_joint_functions.R` | Design, truth, flows, simulation, Stan data, fitting (cmdstanr or rstan) |
| `30_run_sim_joint.R` | Scenario x replicate x model runner; restart-safe, one .rds per fit |
| `31_summarise_sim_joint.R` | Bias, coverage, interval width, sampler summary, Okabe-Ito figures |

Outputs go to `results/sims/joint_cjs_cwt/run2/` (summaries in `summary/`). `run1/` holds the 30 Sep draft test fits (rstan, 2 chains x 400/400, 5 replicates); set `run_dir` to it in 31 to re-summarise them.

Run: `source(here::here("scripts", "sims_joint", "30_run_sim_joint.R"))`, then
`source(here::here("scripts", "sims_joint", "31_summarise_sim_joint.R"))`.
Set `cfg$backend = "cmdstanr"` (default, works on the laptop). Draft run here used rstan,
2 chains x 400/400 because this test machine has 2 cores.

## Process (shared by all fish)
- Reaches: 1 release-Woodson, 2 Woodson-Butte, 3 Butte-Sacramento (KL is inside), 4 Sacramento-Delta.
- Survival per 100 km: `logit s_j = b_j + bQ*Q + bL*L + group RE (+ delta_AT for acoustic fish)`;
  over a reach, `s_j^(km/100)`.
- Travel time: `log days per 100 km = t_j + gQ*Q + gL*L + group RE`; cumulative median to a
  station is the sum over reaches; individual times lognormal around it. Computed on the log
  scale with log_sum_exp (no exp-then-log).
- Q = group flow covariate (mean standardized log flow over 30 days after release);
  L = group mean fork length, standardized on one scale for both sources (80 mm, SD 15).
- Non-centered group REs throughout (stage 1 result).

## Observation models
- **Telemetry:** CJS over Woodson, Butte, Sacramento, Delta; detection by station x year.
  Histories are grouped (16 possible per group) and travel times enter as sufficient statistics,
  so the CJS part is fast.
- **CWT (new):** recoveries at KL in 5-day bins after release:
  `y[g,d] ~ Poisson(N_g x S_KL,g x P(arrive in bin d) x p_KL(d))`.
  Arrival probabilities come from the shared travel-time model. This is the time-stratified
  idea of Hance et al. (2021) at group level (CWT groups are tens of thousands of fish, so
  individual latent states are neither needed nor feasible).
- **KL trap:** `logit p_KL(d) = a0 + year RE + bp x KL flow in bin d`.
- **Efficiency trials:** marked releases just above the trap, `r ~ Binomial(m, p_KL)`, 5 per
  year x 1,500 fish (stand-in for the BT-SPAS-X KL trials).
- **Tag effect:** acoustic fish have survival offset `delta_AT` (tag burden, and by extension
  premature tag failure, Holbrook et al. 2013); CWT fish do not.

## Models compared (same simulated data each replicate)
| Model | Data | KL pCap | Notes |
|---|---|---|---|
| A CWT only | CWT | one constant, prior (like 06) | reaches 1-3 pooled into one rate (not separable from CWT alone) |
| B Joint, constant pCap | CJS + CWT + trials | year RE, not flow-dependent | delta_AT estimated |
| C Joint, flow pCap | CJS + CWT + trials | year RE + flow slope | the proposed model |
| D CJS only | CJS | n/a | like 03; delta_AT not estimable, so S_KL is tagged-fish survival |

## Scenarios
- **S1:** KL efficiency constant in flow (bp = 0).
- **S2:** KL efficiency falls with flow (bp = -0.6 per SD of log flow).

## Truth (draft values)
b = (0.3, 1.0, 0.8, 1.0); bQ = 0.5; bL = 0.3; delta_AT = -0.3; sd_S = 0.5;
t = (1.2, 1.0, 1.0, 1.1); gQ = gL = -0.15; KL efficiency 0.3% at mean flow, year SD 0.3.
Estimand S_KL: untagged 75 mm fish from Battle to KL, typical group, mean flow (truth 0.244).

## Design (reduced for speed)
25 telemetry groups x 150 fish over 2013-2024 (60% Battle, 40% RBDD, 70-130 mm);
300 CWT groups over 22 years (94% Battle 61-87 mm, 6% RBDD 43-60 mm, confounded as in the
real data), about 60,000 released per group; 60-day recovery window.

## Simplifications to remember
- Placeholder river km (Battle 522, RBDD 462, Woodson 395, Butte 318, KL 222, Sacramento 172,
  Delta 71). Replace with 02's ReachKM before any quantitative use.
- Group-level flow covariate (as in Flora's models), not Perry et al. (2018) reach-entry
  time-varying covariates.
- Data simulated from the fitted model's own structure: tests estimability, not robustness.
- The pCap prior is centred on the truth (SD 1.39, the BT-SPAS-X prior SD). In the real data
  the CWT posterior sits about 4x below the prior mean.
- No CWT processing fraction, no fry rearing delay, tag effect constant in size.
- Telemetry travel times: one individual multiplier across stations, so the model's
  per-station likelihood has the same pseudo-replication as Flora's (open item).

## Next scenarios to add
1. No efficiency trials (tests the delta_AT vs p_KL trade-off directly).
2. Unknown CWT processing fraction by year (the 2023 problem).
3. Genuinely different flow slopes for CWT and acoustic fish: does the joint model reveal the
   conflict (source-specific slope, posterior predictive checks)?
4. RBDD fry rearing delay (the 07 travel-time misfit).
5. Mis-centred pCap prior (like the real 4x gap).
6. Flora's soft-sharing combined model as a fifth comparison.
