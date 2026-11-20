# Simple joint acoustic telemetry + CWT simulation

This is the stripped-down version of the joint model.

## Biological process

Both tag types are assumed to sample the same underlying fish population:

- Survival from release to the endpoint is a single probability, `S`.
- Travel time to the endpoint is lognormal:

  `log(T) ~ Normal(mu_logT, sigma_logT^2)`.

No reaches, routing, flow, size, year effects, random effects, or time-varying parameters are included.

## Observation models

For acoustic tags:

`n_AT ~ Binomial(N_AT, S * p_AT)`

For CWT fish:

`n_CWT ~ Binomial(N_CWT, S * p_CWT)`

Observed travel times from detected acoustic fish and recovered CWT fish are then modeled as draws from the same lognormal travel-time distribution.

The key assumption is that `p_AT` and `p_CWT` do not vary with travel time in this first model.

## Why p_AT and p_CWT are known here

With one endpoint, the acoustic data identify `S * p_AT`, not `S` and `p_AT` separately. The CWT data similarly identify `S * p_CWT`.

So this first simulation treats both observation probabilities as known. This is intentional. It lets us ask the simplest question first:

> If both observation systems truly sample the same survival and travel-time process, does a joint likelihood recover the known survival and travel-time parameters?

Later versions can estimate acoustic detection using a downstream receiver and CWT capture efficiency using efficiency trials or an informative external estimate.

## Relation to the uploaded papers

This is a deliberate simplification, not a direct reproduction of any one published model.

- Perry et al. developed the acoustic mark-recapture framework and later jointly modeled survival and lognormal travel times.
- Holbrook et al. wrote acoustic apparent survival as the product of fish survival and tag survival, and used CWT travel-time distributions to inform transmitter-failure correction.
- Newman & Brandes treated Chipps Island CWT recoveries as release-recovery information whose expectation depends on survival and capture/recovery probability.
- Hance et al. extended the telemetry framework to integrate over unobserved travel times. That machinery is unnecessary in this first endpoint-only model.

## Files

- `simple_joint_AT_CWT_JAGS.R`: self-contained simulation, JAGS fit, summary, and figures.
- `simple_joint_AT_CWT.jags`: the short JAGS model by itself.

## Outputs created by the R script

- `simple_joint_parameter_recovery.csv`
- `figure1_travel_time.png`
- `figure2_posterior_recovery.png`

## Requirements

In R:

```r
install.packages(c("rjags", "coda"))
```

JAGS itself must also be installed on the computer.

## Next logical extension

Do not add routes or covariates yet. The next useful extension is tag failure:

`Pr(acoustic observation at time T) = S * p_AT * q(T)`

where `q(T)` is the probability the transmitter is still operating at travel time `T`. This follows the core Holbrook et al. idea and would let the CWT travel-time distribution help correct telemetry survival when slow acoustic migrants disappear because tags expire.
