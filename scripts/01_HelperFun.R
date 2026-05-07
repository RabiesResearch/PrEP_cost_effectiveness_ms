

# =============================================================================
# Helper functions — optimised version
# -----------------------------------------------------------------------------
# ### OPT: Pre-load model files and posterior samples ONCE at package/script
#          load time, not inside predict_cases on every call.
#          Call load_rabies_models() once before running any scenarios.
# -----------------------------------------------------------------------------
.rabies_model_cache <- new.env(parent = emptyenv())

load_rabies_models <- function(
    vax_model_path      = "./data/dog_incidence_model/cases_from_vax_model.rds",
    vax_case_model_path = "./data/dog_incidence_model/cases_from_vax+cases_model.rds"
) {
  vax_model      <- readRDS(vax_model_path)
  vax_case_model <- readRDS(vax_case_model_path)
  .rabies_model_cache$vax_samples      <- posterior_samples(vax_model)[, 1:3]
  .rabies_model_cache$vax_case_samples <- posterior_samples(vax_case_model)[, 1:4]
  invisible(NULL)
}


# Estimate dog population
estimate_dog_population <- function(N, pop, HDR, horizon) {
  HDR_draws <- runif(n = N, min = HDR[1], max = HDR[2])
  # ### OPT: vectorised outer product instead of year-loop
  matrix(pop / HDR_draws, nrow = N, ncol = horizon)
}


# Vaccination coverage ramp
calc_vax_coverage <- function(base_vax_cov,
                              target_vax_cov    = NULL,
                              mdv_campaign_budget = NULL,
                              vaccinate_dog_cost  = NULL,
                              dog_pop             = NULL,
                              horizon,
                              discount,
                              years_to_target,
                              dog_burnin = 0) {
  
  if (is.null(target_vax_cov) && is.null(mdv_campaign_budget))
    stop("At least one of 'target_vax_cov' or 'mdv_campaign_budget' must be provided.")
  if (dog_burnin < 0 || dog_burnin >= horizon)
    stop("'dog_burnin' must be >= 0 and < horizon.")
  
  cost_per_dog <- if (!is.null(vaccinate_dog_cost) && length(vaccinate_dog_cost) > 1)
    mean(vaccinate_dog_cost) else vaccinate_dog_cost
  
  if (!is.null(target_vax_cov))
    target_vax_cov <- min(target_vax_cov, 0.8)
  
  generate_vax_coverage <- function(final_coverage, years_to_target) {
    active_horizon <- horizon - dog_burnin
    if (active_horizon <= 0) return(rep(base_vax_cov, horizon))
    
    # ### OPT: vectorised ramp instead of for-loop
    t      <- seq_len(active_horizon)
    ramp_t <- pmin(t, years_to_target)
    # Geometric-ish interpolation used in original: iterative, so we replicate
    # with Reduce to stay faithful to original logic
    vax_active <- Reduce(
      function(prev, yr) {
        if (yr <= years_to_target)
          prev + (final_coverage - prev) * (yr / years_to_target)
        else
          final_coverage
      },
      seq_len(active_horizon),
      accumulate = TRUE,
      init = base_vax_cov
    )[-1]  # drop the init value
    
    vax_active <- c(base_vax_cov, vax_active[-active_horizon])
    c(rep(base_vax_cov, dog_burnin), vax_active)[seq_len(horizon)]
  }
  
  if (!is.null(target_vax_cov) && is.null(mdv_campaign_budget))
    return(generate_vax_coverage(target_vax_cov, years_to_target))
  
  if (!is.null(mdv_campaign_budget) && is.null(target_vax_cov)) {
    if (is.null(cost_per_dog) || is.null(dog_pop))
      stop("Both 'vaccinate_dog_cost' and 'dog_pop' must be provided.")
    max_achievable <- min(floor(mdv_campaign_budget / cost_per_dog) / dog_pop, 0.8)
    return(generate_vax_coverage(max_achievable, years_to_target))
  }
  
  # Both provided
  if (is.null(cost_per_dog) || is.null(dog_pop))
    stop("Both 'vaccinate_dog_cost' and 'dog_pop' must be provided.")
  max_achievable <- min(floor(mdv_campaign_budget / cost_per_dog) / dog_pop, 0.8)
  if (max_achievable < target_vax_cov)
    message("Budget is insufficient to achieve target coverage; using budget-based maximum coverage.")
  generate_vax_coverage(min(target_vax_cov, max_achievable), years_to_target)
}


# ### OPT: Vectorised — no year-loop; single matrix multiply per output
calculate_vaccinated_and_susceptible <- function(N, horizon, dog_pop, vax_cov) {
  cov_mat            <- matrix(unlist(vax_cov), nrow = N, ncol = horizon, byrow = FALSE)
  ts_dogs_vaccinated <- dog_pop * cov_mat
  list(ts_dogs_vaccinated = ts_dogs_vaccinated,
       sus_dogs           = dog_pop - ts_dogs_vaccinated)
}


# ### OPT: Vectorised — single runif draw of N*horizon values instead of year-loop
calculate_campaign_cost <- function(N, horizon, mdv_campaign_budget,
                                    ts_dogs_vaccinated, vaccinate_dog_cost) {
  if (!is.null(mdv_campaign_budget)) {
    matrix(mdv_campaign_budget, nrow = N, ncol = horizon)
  } else {
    cost_draws <- runif(n = N * horizon,
                        min = vaccinate_dog_cost[1],
                        max = vaccinate_dog_cost[2])
    ts_dogs_vaccinated * matrix(cost_draws, nrow = N, ncol = horizon)
  }
}


# ### OPT: No longer reads RDS or calls posterior_samples — uses pre-loaded cache.
#          The rep-level loop is unavoidable (Markov chain), but the year-level
#          inner loop is kept tight.
predict_cases <- function(nreps, vax_cov, horizon, dog_pop,
                          rabies_inc, seed = NULL) {
  
  if (!exists("vax_samples", envir = .rabies_model_cache))
    stop("Call load_rabies_models() before running the model.")
  
  if (!is.null(seed)) set.seed(seed)
  
  vax_samples      <- .rabies_model_cache$vax_samples
  vax_case_samples <- .rabies_model_cache$vax_case_samples
  incidence_adjust <- 0.0001394842
  
  cases_mat <- matrix(NA_real_, nrow = nreps, ncol = horizon)
  
  for (rep in seq_len(nreps)) {
    # Year 1 — vax-only model
    ps1    <- vax_samples[sample.int(nrow(vax_samples), 1L), ]
    mu_y1  <- exp(ps1[[1]] + ps1[[2]] * vax_cov[[1]] + log(dog_pop[rep]))
    cases_mat[rep, 1] <- min(
      rnbinom(1L, mu = mu_y1, size = ps1[[3]]),
      rabies_inc[2] * dog_pop[rep]
    )
    
    # Years 2+ — vax + lagged-cases model
    ps2 <- vax_case_samples[sample.int(nrow(vax_case_samples), 1L), ]
    for (year in 2:horizon) {
      log_prev <- log(cases_mat[rep, year - 1] / dog_pop[rep] + incidence_adjust)
      mu_yr    <- exp(ps2[[1]] + ps2[[2]] * vax_cov[[year]] + ps2[[3]] * log_prev +
                        log(dog_pop[rep]))
      cases_mat[rep, year] <- min(
        rnbinom(1L, mu = mu_yr, size = ps2[[4]]),
        rabies_inc[2] * dog_pop[rep]
      )
    }
  }
  cases_mat
}


nBitesBiters <- function(n_rabid, pBite, pBiteK) {
  bites_by_dog <- rnbinom(n = n_rabid, mu = pBite, size = pBiteK)
  list(
    nBites      = sum(bites_by_dog),
    nBiters     = sum(bites_by_dog > 0L),
    bites_by_dog = bites_by_dog[bites_by_dog > 0L]
  )
}


predict_dograbies_n_humanexposures <- function(N, horizon, vax_cov, dog_pop,
                                               rabies_inc, mu, k, seed = NULL) {
  ts_rabid_dogs <- predict_cases(nreps = N, vax_cov = vax_cov, horizon = horizon,
                                 dog_pop = dog_pop, rabies_inc = rabies_inc, seed = seed)
  
  output               <- sapply(ts_rabid_dogs, FUN = nBitesBiters, pBite = mu, pBiteK = k)
  ts_exposures         <- matrix(unlist(output[1, ]), N, horizon)
  ts_rabid_biting_dogs <- matrix(unlist(output[2, ]), N, horizon)
  bites_by_dog         <- output[3, ]
  
  list(ts_rabid_dogs        = ts_rabid_dogs,
       ts_exposures         = ts_exposures,
       ts_rabid_biting_dogs = ts_rabid_biting_dogs,
       bites_by_dog         = bites_by_dog)
}


split_dog_population <- function(dog_pop, n_splits) {
  base <- dog_pop %/% n_splits
  rem  <- dog_pop %% n_splits
  lapply(seq_len(n_splits), function(i) base + (rem >= i))
}


combine_rabies_results <- function(results_list) {
  list(
    ts_rabid_dogs        = Reduce(`+`, lapply(results_list, `[[`, "ts_rabid_dogs")),
    ts_exposures         = Reduce(`+`, lapply(results_list, `[[`, "ts_exposures")),
    ts_rabid_biting_dogs = Reduce(`+`, lapply(results_list, `[[`, "ts_rabid_biting_dogs")),
    bites_by_dog         = Reduce(function(a, b) Map(c, a, b),
                                  lapply(results_list, `[[`, "bites_by_dog"))
  )
}


predict_dograbies_split <- function(N, horizon, vax_cov, dog_pop, rabies_inc,
                                    mu, k, seed,
                                    pop_serengeti = 1e5,
                                    split_by = c("max", "first", "mean")) {
  split_by <- match.arg(split_by)
  ref_dogs <- switch(split_by,
                     max   = max(dog_pop,  na.rm = TRUE),
                     first = dog_pop[1, 1],
                     mean  = mean(dog_pop, na.rm = TRUE))
  
  n_splits       <- max(100L, ceiling(ref_dogs / pop_serengeti))
  dog_pop_splits <- split_dog_population(dog_pop, n_splits)
  
  split_results <- lapply(seq_len(n_splits), function(i) {
    predict_dograbies_n_humanexposures(
      N = N, horizon = horizon, vax_cov = vax_cov,
      dog_pop = dog_pop_splits[[i]],
      rabies_inc = rabies_inc, mu = mu, k = k,
      seed = seed + i - 1L
    )
  })
  
  combine_rabies_results(split_results)
}




###########
# Post model helper functions


summarise_stochasticity <- function(mat=status_quo$ts_deaths_averted_PEP, scenario="status quo") {
  # Calculate summary statistics per column
  result <- apply(mat, 2, function(x) {
    c(
      LL = quantile(x, 0.025, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      UL = quantile(x, 0.975, na.rm = TRUE)
    )
  })
  
  # Convert result to a data frame
  result_df <- as.data.frame(t(result)) %>%
    dplyr::rename(
      LL = `LL.2.5%`, 
      UL = `UL.97.5%`
    )
  
  # Add 'year' and 'scenario' columns
  result_df$year <- seq_len(ncol(mat))  
  result_df$scenario <- scenario
  
  return(result_df)
}

summarise_across_horizon <- function(
    mat,
    probs = c(0.025, 0.5, 0.95),
    na.rm = TRUE
) {
  
  totals <- rowSums(mat)
  
  qs <- quantile(totals, probs = probs, na.rm = na.rm)
  
  tibble::tibble(
    LL     = as.numeric(qs[1]),
    Median = as.numeric(qs[2]),
    UL     = as.numeric(qs[3])
  )
}

summarise_across_scenarios <- function(country_results, matrix_names) {
  
  output_list <- list()
  
  for (mat_name in matrix_names) {
    
    scen_summaries <- lapply(names(country_results), function(scenario) {
      
      mat <- country_results[[scenario]][[mat_name]]
      
      df <- summarise_stochasticity(mat, scenario = scenario)
      df$scenario <- scenario
      df$matrix   <- mat_name    # optional label
      
      df
    })
    
    # bind scenarios for this matrix
    combined <- dplyr::bind_rows(scen_summaries)
    
    # store under matrix name
    output_list[[mat_name]] <- combined
  }
  
  return(output_list)
}


summarise_variables_across_scenarios <- function(results_list,
                                                 variables,
                                                 scenario_names = names(results_list)) {
  
  if (length(results_list) != length(scenario_names)) {
    stop("The number of scenarios in results_list must match scenario_names.")
  }
  
  per_var_tables <- map(variables, function(v) {
    mats <- map(results_list, ~ .x[[v]])
    
    if (any(map_lgl(mats, is.null))) {
      missing_in <- scenario_names[map_lgl(mats, is.null)]
      stop(sprintf("Variable '%s' missing in scenarios: %s", v, paste(missing_in, collapse = ", ")))
    }
    
    create_summary_across_horizon(
      variable_name  = v,
      matrices       = mats,
      scenario_names = scenario_names
    ) %>%
      select(scenario, LL, Median, UL) %>%
      rename(
        !!paste0(v, "_LL")     := LL,
        !!paste0(v, "_Median") := Median,
        !!paste0(v, "_UL")     := UL
      )
  })
  
  reduce(per_var_tables, full_join, by = "scenario")
}

summarise_multiple_vars <- function(variables, scenario_list, scenario_names) {
  
  expand_grid(
    scenario = scenario_names,
    variable = variables
  ) %>%
    mutate(
      data = purrr::map2(scenario, variable, ~ {
        mat <- scenario_list[[.x]][[.y]]
        
        if (is.null(mat) || is.null(dim(mat))) {
          stop("Missing matrix: ", .x, " / ", .y)
        }
        
        summarise_stochasticity(mat) %>%
          dplyr::select(-scenario)
      })
    ) %>%
    tidyr::unnest(data)
}


# How to use
# long_all <- summarise_multiple_vars(
#   variables = c("ts_exposures", "ts_deaths"),
#   scenario_list = Kerala_general,
#   scenario_names = scenarios
# )




# Summarise and visualize
create_summary_across_horizon <- function(variable_name, matrices, scenario_names) {
  
  map2_dfr(matrices, scenario_names, function(mat, scen_name) {
    
    # Sum across years WITHIN each simulation (rowSums if rows = simulations)
    # then take quantiles across simulations
    sim_totals <- rowSums(mat)
    
    tibble(
      scenario = scen_name,
      LL       = quantile(sim_totals, 0.025),
      Median   = median(sim_totals),
      UL       = quantile(sim_totals, 0.975)
    )
  })
}



plot_ribbon <- function(mydata, x_axis, y_axis, xlab, ylab, palette = NULL) {
  
  #  color palette if none is provided
  if (is.null(palette)) {
    palette <- c("#8195b2", "#648c67", "orange2", "turquoise1") 
  }
  
  ggplot(mydata, aes(x = {{ x_axis }}, y = {{ y_axis }}, group = scenario, color = scenario, fill = scenario)) +  
    geom_line(size = 1) +  
    geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.4, color = NA) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    scale_y_continuous(labels = scales::comma) +
    scale_color_manual(values = palette) +  # Apply custom colors
    scale_fill_manual(values = palette) +   # Apply custom fill colors
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
}


# --- Function to derive the ICER frontier (not Pareto) ---
get_icer_frontier <- function(df, effect_col = "x", cost_col = "inc_cost_musd") {
  
  out <- df %>%
    filter(!is.na(.data[[effect_col]]), !is.na(.data[[cost_col]])) %>%
    distinct(scenario, .keep_all = TRUE) %>%
    arrange(.data[[effect_col]], .data[[cost_col]])
  
  # 1) Remove strictly dominated strategies
  changed <- TRUE
  while (changed) {
    n_before <- nrow(out)
    keep <- rep(TRUE, nrow(out))
    
    for (i in seq_len(nrow(out))) {
      for (j in seq_len(nrow(out))) {
        if (i == j) next
        
        eff_i  <- out[[effect_col]][i]
        cost_i <- out[[cost_col]][i]
        eff_j  <- out[[effect_col]][j]
        cost_j <- out[[cost_col]][j]
        
        if (eff_j >= eff_i &&
            cost_j <= cost_i &&
            (eff_j > eff_i || cost_j < cost_i)) {
          keep[i] <- FALSE
          break
        }
      }
    }
    
    out <- out[keep, , drop = FALSE]
    changed <- nrow(out) < n_before
  }
  
  # 2) Remove extendedly dominated strategies
  changed <- TRUE
  while (changed && nrow(out) >= 3) {
    changed <- FALSE
    out <- out %>% arrange(.data[[effect_col]], .data[[cost_col]])
    
    d_cost <- diff(out[[cost_col]])
    d_eff  <- diff(out[[effect_col]])
    icers  <- d_cost / d_eff
    
    # ICERs on the frontier should increase monotonically
    bad <- which(diff(icers) < 0)
    
    if (length(bad) > 0) {
      # remove the middle strategy causing the violation
      remove_idx <- bad[1] + 1
      out <- out[-remove_idx, , drop = FALSE]
      changed <- TRUE
    }
  }
  
  out %>%
    arrange(.data[[effect_col]], .data[[cost_col]]) %>%
    mutate(
      frontier_inc_cost   = c(NA, diff(.data[[cost_col]])),
      frontier_inc_effect = c(NA, diff(.data[[effect_col]])),
      frontier_icer       = frontier_inc_cost / frontier_inc_effect
    )
}





