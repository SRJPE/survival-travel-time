# Joint CJS + CWT simulation: first results (draft, 30 Sep 2026)

5 replicates per scenario, 4 models per replicate, 40 fits. rstan, 2 chains x 400 warmup +
400 sampling (test machine has 2 cores). 0 divergences, 0 treedepth hits in all 40 fits.
Max Rhat up to 1.09 (model A) and about 1.03 to 1.06 elsewhere: expected with 800 draws;
the laptop settings (4 x 750 + 750) should clear it. Treat these as a first look, not
evidence: 5 replicates, reduced design, data simulated from the fitted model's structure.

Models: A CWT only, constant pCap (like 06); B joint, constant pCap + trials; C joint,
flow-dependent pCap + trials; D CJS only (like 03).
Scenarios: S1 KL efficiency constant in flow; S2 KL efficiency falls with flow (bp = -0.6).

## Flow effect on survival (bQ, truth 0.5)
| Scenario | A CWT only | B joint, const pCap | C joint, flow pCap | D CJS only |
|---|---|---|---|---|
| S1 mean (coverage, width) | 0.40 (0.4, 0.17) | 0.48 (1.0, 0.18) | 0.49 (1.0, 0.18) | 0.42 (1.0, 0.49) |
| S2 mean (coverage, width) | **-0.03** (0.0, 0.11) | **0.24** (0.0, 0.18) | **0.48** (0.8, 0.19) | 0.57 (1.0, 0.53) |

## Survival Battle to KL, untagged 75 mm fish (S_KL, truth 0.244)
| Scenario | A | B | C | D |
|---|---|---|---|---|
| S1 mean (coverage, width) | 0.18 (0.2, 0.11) | 0.24 (0.8, 0.06) | 0.24 (0.8, 0.06) | 0.21 (0.8, 0.17) |
| S2 mean (coverage, width) | 0.08 (0.0, 0.09) | 0.23 (0.6, 0.06) | 0.25 (1.0, 0.06) | 0.19 (0.8, 0.18) |

D estimates tagged-fish survival (it cannot separate the tag effect), so it runs low for untagged fish.

## Other parameters
- KL efficiency flow slope bp (C): S1 0.005 (truth 0); S2 -0.58 (truth -0.6, coverage 0.8, width 0.11).
- Acoustic tag offset delta_AT (truth -0.3): C -0.25 (S1), -0.28 (S2), unbiased but wide (width
  about 0.5). B in S2: -0.16.
- Size effect bL (truth 0.3): B and C 0.27 to 0.31 (width 0.18 to 0.21); D 0.20 to 0.22 (width
  0.58 to 0.72); A 0.21 to 0.26.
- Travel time (TT_KL, gQ): recovered by all models with CWT data; D about 6x wider on TT_KL.

## Reading
1. **A flow-dependent KL trap reproduces the real CWT result.** With trap efficiency falling with
   flow (S2), the CWT-only model (like 06) estimates a flow effect near zero (-0.03) when the
   truth is 0.5, with a narrow interval that confidently excludes the truth. The real 06 gives
   0.09. This doesn't prove that is what's happening in the real data, but it shows a mundane,
   testable mechanism that is sufficient to produce it.
2. **Joining the data with the wrong observation model makes things worse, not better.** B
   (joint, constant pCap) pulls the shared flow slope to 0.24 in S2: worse than the CJS alone
   (0.57), with a narrow interval that excludes the truth in every replicate. So soft or hard
   sharing without modelling the trap is a trap.
3. **With a flow-dependent pCap and efficiency trials, the joint model recovers the flow effect
   (0.48) and costs nothing when efficiency is actually constant (S1: 0.49, bp estimated at 0).**
4. **Precision gain from joining:** flow-effect interval about 2.8x narrower than CJS only
   (0.19 vs 0.53); survival-to-KL interval about 3x narrower (0.06 vs 0.17); size effect about
   3x narrower.
5. **The acoustic tag effect is only weakly identified** (interval width about 0.5 on the logit
   scale), even with 5 efficiency trials per year for 22 years. This is the delta_AT vs p_KL
   trade-off: the trials are what separate them. Next scenario: no trials.
6. Model A misses the survival level even in S1 (0.18 vs 0.24) with the pCap prior centred on
   the truth: CWT alone cannot pin the level, consistent with the last-interval point.

## Next
- Laptop run: 20 replicates, cmdstanr, 4 x 750/750 (about 4 fits per replicate x 2 scenarios;
  roughly 10 hours, so overnight on a free machine, not the night 09 runs).
- Add scenarios from the README list, first "no efficiency trials" and "different true flow
  slopes for CWT and acoustic fish" (can the joint model reveal a genuine conflict?).
- Replace placeholder river km with 02's values; use the real release-group design.
