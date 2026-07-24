# Convergence diagnostics for CJS model.

rm(list = ls())

library(rstan)
library(bayesplot)
library(dplyr)
library(ggplot2)
library(here)
library(tibble)

dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)

rhat_threshold <- 1.02
ess_threshold <- 400
core_trace_pars <- c(
  "S_bReach[1]", "S_bReach[2]", "S_bReach[3]", "S_bReach[4]",
  "S_bTrib[1]", "S_bTrib[2]",
  "S_bCov", "S_bCovT", "S_bSz",
  "RE_sd[1]", "RE_sd[2]", "RE_sd[3]", "RE_sd[4]",
  "RE_sdT[1]", "RE_sdT[2]",
  "TT_bCov", "T_bReach", "T_bSz",
  "Pro_sd", "Pro_sdT"
)

# Change this one block to diagnose a different CJS fit.
fit_config <- list(
  fit_file = here("results", "fit_CovIndCont_MaxFlow_FL.Rdata"),
  output_tag = "MaxFlowFL",
  model_label = "MaxFlow + FL",
  core_trace_pars = core_trace_pars,
  # Keep the default diagnostics fast; trace/pairs plots can be enabled for closer inspection.
  write_posterior_plots = TRUE,
  write_trace_plots = FALSE,
  write_pairs_plots = FALSE
)

# Example alternative:
# fit_config <- list(
#   fit_file = here("outputs", "fit_CovIndCont_syre_MaxFlow_FL.Rdata"),
#   output_tag = "Syre_MaxFlowFL",
#   model_label = "Syre MaxFlow + FL",
#   core_trace_pars = core_trace_pars_syre
# )

run_diagnostics <- function(config) {
  fit_file <- config$fit_file
  output_tag <- config$output_tag
  model_label <- config$model_label
  write_posterior_plots <- is.null(config$write_posterior_plots) || isTRUE(config$write_posterior_plots)
  write_trace_plots <- is.null(config$write_trace_plots) || isTRUE(config$write_trace_plots)
  write_pairs_plots <- is.null(config$write_pairs_plots) || isTRUE(config$write_pairs_plots)
  
  if (!file.exists(fit_file)) {
    warning("Skipping ", fit_file, ": file does not exist.")
    return(invisible(NULL))
  }
  
  load(fit_file)
  if (!exists("fit")) {
    stop("Expected object named `fit` in ", fit_file)
  }
  
  available_pars <- fit@sim$fnames_oi
  n_warmup <- fit@sim$warmup
  if (length(n_warmup) > 1) {
    n_warmup <- unique(n_warmup)
  }
  if (length(n_warmup) != 1 || !is.finite(n_warmup)) {
    warning("Could not determine warmup from fit; using 0 for trace plots.")
    n_warmup <- 0
  }
  
  keep_existing_pars <- function(pars) {
    pars[pars %in% available_pars]
  }

  keep_existing_requested_pars <- function(pars) {
    pars[vapply(
      pars,
      function(par) par %in% available_pars || any(startsWith(available_pars, paste0(par, "["))),
      logical(1)
    )]
  }
  
  keep_existing_prefixes <- function(prefixes) {
    prefixes[vapply(
      prefixes,
      function(prefix) any(startsWith(available_pars, paste0(prefix, "[")) | available_pars == prefix),
      logical(1)
    )]
  }

  first_existing_by_prefix <- function(prefixes, n_per_prefix = 1) {
    unique(unlist(lapply(prefixes, function(prefix) {
      matches <- available_pars[
        startsWith(available_pars, paste0(prefix, "[")) |
          available_pars == prefix
      ]
      head(matches, n_per_prefix)
    }), use.names = FALSE))
  }
  
  model_parameter_prefixes <- c(
    "P_b", "muPb", "sdPb",
    "S_bReach", "S_bCov", "S_bCovT", "S_bSz", "S_bTrib",
    "S_RE", "S_REt", "S_REy", "S_REs",
    "RE_sd", "RE_sdT", "REy_sd", "REs_sd",
    "T_bReach", "T_bCov", "TT_bCov", "TT_bCovT", "T_bSz", "T_bTrib",
    "TT_RE", "TT_RET", "TT_REy", "TT_REs",
    "TTRE_sd", "TTRE_sdT", "TTREy_sd", "TTREs_sd",
    "TT", "TTT", "lgB0", "lgB0T", "Pro_sd", "Pro_sdT"
  )
  diagnostic_pars <- keep_existing_prefixes(model_parameter_prefixes)
  core_trace_pars <- keep_existing_pars(config$core_trace_pars)
  diagnostic_pars <- unique(c(core_trace_pars, diagnostic_pars))
  
  if (length(diagnostic_pars) == 0) {
    warning("Skipping ", fit_file, ": no diagnostic parameters found.")
    return(invisible(NULL))
  }
  
  model_summary <- summary(fit, pars = diagnostic_pars, probs = c(0.025, 0.5, 0.975))$summary %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`) %>%
    mutate(
      low_ess = n_eff < ess_threshold,
      high_rhat = is.na(Rhat) | Rhat > rhat_threshold
    ) %>%
    arrange(desc(high_rhat), desc(low_ess), Rhat, n_eff)
  
  flagged_summary <- model_summary %>%
    filter(high_rhat | low_ess)
  
  core_model_summary <- model_summary %>%
    filter(parameter %in% core_trace_pars)
  
  core_flagged_summary <- core_model_summary %>%
    filter(high_rhat | low_ess)

  write.csv(
    model_summary,
    here("results/tables", paste0("cjs_convergence_diagnostics_", output_tag, ".csv")),
    row.names = FALSE
  )
  write.csv(
    flagged_summary,
    here("results/tables", paste0("cjs_convergence_diagnostics_flagged_", output_tag, ".csv")),
    row.names = FALSE
  )
  write.csv(
    core_model_summary,
    here("results/tables", paste0("cjs_core_convergence_diagnostics_", output_tag, ".csv")),
    row.names = FALSE
  )
  write.csv(
    core_flagged_summary,
    here("results/tables", paste0("cjs_core_convergence_diagnostics_flagged_", output_tag, ".csv")),
    row.names = FALSE
  )
  
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  n_divergent <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
  n_transitions <- sum(vapply(sampler_params, nrow, integer(1)))
  n_treedepth <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 10), numeric(1)))
  
  sampler_diagnostics <- tibble(
    diagnostic = c(
      "post_warmup_transitions",
      "divergent_transitions",
      "divergent_transition_rate",
      "max_treedepth_hits",
      "max_treedepth_hit_rate"
    ),
    value = c(
      n_transitions,
      n_divergent,
      n_divergent / n_transitions,
      n_treedepth,
      n_treedepth / n_transitions
    )
  )
  
  sampler_parameter_summary <- bind_rows(
    lapply(seq_along(sampler_params), function(chain_id) {
      as.data.frame(sampler_params[[chain_id]]) %>%
        summarise(
          across(
            everything(),
            list(mean = mean, min = min, max = max),
            .names = "{.col}_{.fn}"
          )
        ) %>%
        mutate(chain = chain_id, .before = 1)
    })
  )
  
  write.csv(
    sampler_diagnostics,
    here("results/tables", paste0("cjs_sampler_diagnostics_", output_tag, ".csv")),
    row.names = FALSE
  )
  write.csv(
    sampler_parameter_summary,
    here("results/tables", paste0("cjs_sampler_parameter_summary_", output_tag, ".csv")),
    row.names = FALSE
  )
  
  rhat_plot <- mcmc_rhat(model_summary$Rhat) +
    ggtitle(paste("Rhat diagnostics:", model_label))
  ggsave(
    here("results/figures", paste0("rhat_", output_tag, ".png")),
    plot = rhat_plot,
    width = 8,
    height = 5,
    dpi = 350
  )
  
  neff_plot <- ggplot(model_summary, aes(x = n_eff)) +
    geom_histogram(bins = 40, fill = "grey70", color = "white") +
    geom_vline(xintercept = ess_threshold, color = "red", linetype = "dashed") +
    theme_bw() +
    labs(
      x = "Effective sample size",
      y = "Parameter count",
      title = paste("Effective sample size diagnostics:", model_label)
    )
  ggsave(
    here("results/figures", paste0("n_eff_", output_tag, ".png")),
    plot = neff_plot,
    width = 8,
    height = 5,
    dpi = 350
  )
  
  np <- nuts_params(fit)
  lp <- log_posterior(fit)
  
  divergence_plot <- mcmc_nuts_divergence(np, lp)
  ggsave(
    here("results/figures", paste0("divergence_", output_tag, ".png")),
    plot = divergence_plot,
    width = 8,
    height = 5,
    dpi = 350
  )
  
  energy_plot <- mcmc_nuts_energy(np)
  ggsave(
    here("results/figures", paste0("energy_", output_tag, ".png")),
    plot = energy_plot,
    width = 8,
    height = 5,
    dpi = 350
  )
  
  write_trace_pdf <- function(filename, pars, title, nrow = NULL) {
    pars <- keep_existing_pars(pars)
    if (length(pars) == 0) {
      message("Skipping ", filename, ": no requested parameters found in fit.")
      return(invisible(NULL))
    }
    
    pdf(here("results/figures", filename), onefile = TRUE, width = 11, height = 8.5)
    trace_args <- list(fit, pars = pars, n_warmup = n_warmup)
    if (!is.null(nrow)) {
      trace_args$facet_args <- list(nrow = nrow)
    }
    trace_plot <- do.call(mcmc_trace, trace_args) +
      ggtitle(title)
    print(trace_plot)
    dev.off()
  }
  
  write_pairs_pdf <- function(filename, pars, title) {
    pars <- keep_existing_pars(pars)
    if (length(pars) < 2) {
      message("Skipping ", filename, ": fewer than two requested parameters found in fit.")
      return(invisible(NULL))
    }
    
    posterior_array <- rstan::extract(
      fit,
      pars = pars,
      permuted = FALSE,
      inc_warmup = FALSE
    )
    pdf(here("results/figures", filename), onefile = TRUE, width = 11, height = 8.5)
    pairs_plot <- mcmc_pairs(
      posterior_array,
      np = np,
      pars = pars,
      off_diag_args = list(size = 0.75),
      grid_args = list(top = title)
    )
    print(pairs_plot)
    dev.off()
  }

  write_posterior_distributions <- function(filename, pars, title) {
    pars <- keep_existing_requested_pars(pars)
    if (length(pars) == 0) {
      message("Skipping ", filename, ": no requested parameters found in fit.")
      return(invisible(NULL))
    }

    posterior_array <- rstan::extract(
      fit,
      pars = pars,
      permuted = FALSE,
      inc_warmup = FALSE
    )
    posterior_plot <- mcmc_dens(
      posterior_array,
      pars = pars,
      facet_args = list(scales = "free")
    ) +
      ggtitle(title)
    ggsave(
      here("results/figures", filename),
      plot = posterior_plot,
      width = 12,
      height = max(8, 0.45 * length(pars)),
      dpi = 350
    )
  }
  
  core_pairs_pars <- core_trace_pars[seq_len(min(length(core_trace_pars), 12))]

  if (write_posterior_plots) {
    write_posterior_distributions(
      paste0("posterior_core_", output_tag, ".png"),
      core_trace_pars,
      paste("Core parameter posterior distributions:", model_label)
    )
  }
  
  if (write_pairs_plots && length(core_pairs_pars) >= 2) {
    posterior_array <- rstan::extract(
      fit,
      pars = core_pairs_pars,
      permuted = FALSE,
      inc_warmup = FALSE
    )
    png(
      filename = here("results/figures", paste0("cjs_pairs_diagnostics_", output_tag, ".png")),
      width = 1200,
      height = 900,
      res = 130
    )
    pairs_plot <- mcmc_pairs(
      posterior_array,
      np = np,
      pars = core_pairs_pars,
      off_diag_args = list(size = 0.75),
      grid_args = list(top = paste("Core parameter pairs:", model_label))
    )
    print(pairs_plot)
    dev.off()
  } else if (write_pairs_plots) {
    message("Skipping cjs_pairs_diagnostics_", output_tag, ".png: fewer than two core parameters found.")
  }
  
  if (write_trace_plots) {
    write_trace_pdf(
      paste0("trace_core_", output_tag, ".pdf"),
      core_trace_pars,
      paste("Core parameter trace plots:", model_label)
    )
    
    write_trace_pdf(
      paste0("trace_pcap_", output_tag, ".pdf"),
      c("P_b[1,1]", "P_b[1,2]", "P_b[1,3]", "P_b[1,4]"),
      paste("Detection probability trace plots:", model_label),
      nrow = 4
    )
    
    write_trace_pdf(
      paste0("trace_survival_", output_tag, ".pdf"),
      first_existing_by_prefix(
        c("S_bReach", "S_bTrib", "S_bCov", "S_bCovT", "S_bSz", "S_RE", "S_REt"),
        n_per_prefix = 2
      ),
      paste("Survival parameter trace plots:", model_label),
      nrow = 4
    )
    
    write_trace_pdf(
      paste0("trace_travel_time_", output_tag, ".pdf"),
      first_existing_by_prefix(
        c("T_bReach", "T_bTrib", "TT_bCov", "TT_bCovT", "T_bSz", "TT_RE", "TT_RET"),
        n_per_prefix = 2
      ),
      paste("Travel-time parameter trace plots:", model_label),
      nrow = 4
    )
  }
  
  if (write_pairs_plots) {
    write_pairs_pdf(
      paste0("pairs_core_", output_tag, ".pdf"),
      core_pairs_pars,
      paste("Core parameter pairs:", model_label)
    )
  }
  
  message("Diagnostics complete for ", fit_file)
  message("Flagged parameters: ", nrow(flagged_summary))
  message("Divergences: ", n_divergent, " of ", n_transitions, " post-warmup transitions")
  
  invisible(
    list(
      fit_file = fit_file,
      output_tag = output_tag,
      n_flagged = nrow(flagged_summary),
      n_divergent = n_divergent,
      n_transitions = n_transitions
    )
  )
}

diagnostic_results <- run_diagnostics(fit_config)

