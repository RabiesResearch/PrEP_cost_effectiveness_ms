

require(pacman)
pacman::p_load(tidyverse,   # cleaning, wrangling
               brms,      # models
               rlang,
               lubridate, # dates
               matrixStats,
               scales,
               readxl,
               ggrepel,
               cowplot,
               patchwork,
               purrr,
               furrr
)



# Decision tree model that can be applied to create different scenarios 

# source helper functions  
source("./scripts/01_HelperFun.R")


decision_tree_PrEP <- function(N = 10, pop = 35e6, HDR = c(16,17), unowned_prop = 0.633,
                               horizon = 5, discount = 0.03, mu = 0.38, k = 0.72,
                               bpi = 15.3, pCompliance_healthy = 0.9654, pStart_healthy = 0.9684211,
                               pSeek_exposure = 0.9, pStart_exposure = 0.9684211,
                               pCompliance_exp = 0.9654, pDeath = 0.17,
                               pPrevent_complete = 0.999, pPrevent_incomplete = 0.986,
                               rabies_inc = c(0.0075, 0.0125),
                               mdv_unowned_budget = NULL, mdv_owned_budget = NULL,
                               vaccinate_owned_dog_cost = c(0.5, 1),
                               vaccinate_unowned_dog_cost = c(3.5, 4.5),
                               base_vax_cov_owned = 0.5,   target_vax_cov_owned = 0.5,
                               base_vax_cov_unowned = 0.03, target_vax_cov_unowned = 0.03,
                               years_to_target = 3,
                               pInvestigate = 0.5, pFound = 0.4, pTestable = 0.2,
                               pFS = 0.05, RIG_cov = 0.38,
                               init_PrEP_cov = 0.15, PrEP_effectiveness = 0,
                               PrEP_vials_per_pt = 0.44, PEP_vials_per_pt = 0.66,
                               human_vaccine_cost_per_vial_PEP = 5, RIG_cost = 12,
                               human_vaccine_cost_per_vial_PrEP=5, 
                               birth_rate = 0.012, routine_vax_cov = 0.78,
                               seed = 123, ibcm = "no", dog_burnin = 1) {
  
  total_horizon <- horizon + dog_burnin
  
  # ---------------------------------------------------------------------------#
  # 1: Dog population & vaccination   #######
  # ---------------------------------------------------------------------------#
  dog_pop      <- estimate_dog_population(N, pop, HDR, total_horizon)
  unowned_dogs <- matrix(rbinom(N * total_horizon, as.integer(dog_pop), unowned_prop),
                         nrow = N, ncol = total_horizon)
  owned_dogs   <- dog_pop - unowned_dogs
  owned_prop   <- 1 - unowned_prop
  
  vax_cov_owned <- calc_vax_coverage(
    base_vax_cov = base_vax_cov_owned, target_vax_cov = target_vax_cov_owned,
    mdv_campaign_budget = mdv_owned_budget, vaccinate_dog_cost = vaccinate_owned_dog_cost,
    dog_pop = owned_dogs, horizon = total_horizon, discount = discount,
    years_to_target = years_to_target, dog_burnin = dog_burnin)
  
  vax_cov_unowned <- calc_vax_coverage(
    base_vax_cov = base_vax_cov_unowned, target_vax_cov = target_vax_cov_unowned,
    mdv_campaign_budget = mdv_unowned_budget, vaccinate_dog_cost = vaccinate_unowned_dog_cost,
    dog_pop = unowned_dogs, horizon = total_horizon, discount = discount,
    years_to_target = years_to_target, dog_burnin = dog_burnin)
  
  dog_vax_cov <- (vax_cov_unowned * unowned_prop) + (vax_cov_owned * owned_prop)
  
  vax_results        <- calculate_vaccinated_and_susceptible(N, total_horizon, dog_pop, dog_vax_cov)
  vaccinated_unowned <- calculate_vaccinated_and_susceptible(N, total_horizon, unowned_dogs, vax_cov_unowned)$ts_dogs_vaccinated
  vaccinated_owned   <- calculate_vaccinated_and_susceptible(N, total_horizon, owned_dogs,   vax_cov_owned)$ts_dogs_vaccinated
  
  MDV_campaign_cost <-
    calculate_campaign_cost(N, total_horizon, mdv_unowned_budget, vaccinated_unowned, vaccinate_unowned_dog_cost) +
    calculate_campaign_cost(N, total_horizon, mdv_owned_budget,   vaccinated_owned,   vaccinate_owned_dog_cost)
  
  # ---------------------------------------------------------------------------#
  # 2: Dogs and human bites: exposures & healthy
  # 2A: Dog rabies & human exposures — strip burn-in immediately   #######
  # ---------------------------------------------------------------------------#
  rabies_results_full <- predict_dograbies_split(
    N = N, horizon = total_horizon, vax_cov = dog_vax_cov, dog_pop = dog_pop,
    rabies_inc = rabies_inc, mu = mu, k = k, seed = seed,
    pop_serengeti = 1e5, split_by = "mean")
  
  keep_cols   <- (dog_burnin + 1):total_horizon
  strip_burnin <- function(mat) mat[, keep_cols, drop = FALSE]
  
  # ### OPT: bites_by_dog flat-list stripping (column-major: year = ceiling(idx/N))
  bites_keep_idx <- which(ceiling(seq_along(rabies_results_full$bites_by_dog) / N) > dog_burnin)
  
  rabies_results <- list(
    ts_rabid_dogs        = strip_burnin(rabies_results_full$ts_rabid_dogs),
    ts_exposures         = strip_burnin(rabies_results_full$ts_exposures),
    ts_rabid_biting_dogs = strip_burnin(rabies_results_full$ts_rabid_biting_dogs),
    bites_by_dog         = rabies_results_full$bites_by_dog[bites_keep_idx])
  
  vax_results$ts_dogs_vaccinated <- strip_burnin(vax_results$ts_dogs_vaccinated)
  vax_results$sus_dogs           <- strip_burnin(vax_results$sus_dogs)
  MDV_campaign_cost              <- strip_burnin(MDV_campaign_cost)
  
  # ---------------------------------------------------------------------------#
  #  2B: Healthy bites (human horizon only)   #######
  # ---------------------------------------------------------------------------#
  ts_recorded_bites <- matrix(rbinom(N * horizon, pop, bpi / 1000), N, horizon)
  
  # ---------------------------------------------------------------------------#
  # 3: PrEP coverage (human horizon only)
  # ---------------------------------------------------------------------------#
  
  annual_PrEP_cov <- rep(0, horizon)
  if (init_PrEP_cov > 0) {
    annual_PrEP_cov[1] <- init_PrEP_cov
    if (horizon > 1) annual_PrEP_cov[-1] <- birth_rate * routine_vax_cov
  }
  cum_PrEP_cov <- cumsum(annual_PrEP_cov)
  
  ts_pop_PrEP     <- matrix(rbinom(N * horizon, pop, rep(annual_PrEP_cov, each = N)), N, horizon)
  ts_cum_pop_PrEP <- rowCumsums(ts_pop_PrEP)
  
  # ---------------------------------------------------------------------------#
  #  4: Healthcare seeking #######
  # ---------------------------------------------------------------------------#
  cum_prep_frac <- as.vector(ts_cum_pop_PrEP / pop)   #: compute once, reuse
  
  ts_exp_PrEP   <- matrix(rbinom(N * horizon, as.vector(rabies_results$ts_exposures), cum_prep_frac), N, horizon)
  ts_exp_noPrEP <- rabies_results$ts_exposures - ts_exp_PrEP
  
  ts_exp_PrEP_seek_care   <- matrix(rbinom(N * horizon, as.vector(ts_exp_PrEP),   pSeek_exposure), N, horizon)
  ts_exp_noPrEP_seek_care <- matrix(rbinom(N * horizon, as.vector(ts_exp_noPrEP), pSeek_exposure), N, horizon)
  ts_exp_seek_care        <- ts_exp_PrEP_seek_care + ts_exp_noPrEP_seek_care
  
  ts_exp_do_not_seek_care <- matrix(rbinom(N * horizon, as.vector(rabies_results$ts_exposures), 1 - pSeek_exposure), N, horizon)
  
  ts_exp_PrEP_do_not_seek_care  <- matrix(rbinom(N * horizon, as.vector(ts_exp_do_not_seek_care), cum_prep_frac), N, horizon)
  ts_exp_noPrEP_do_not_seek_care <- ts_exp_noPrEP - ts_exp_noPrEP_seek_care
  
  ts_healthy_seek_care        <- ts_recorded_bites - ts_exp_seek_care
  ts_healthy_PrEP_seek_care   <- matrix(rbinom(N * horizon, as.vector(ts_healthy_seek_care), cum_prep_frac), N, horizon)
  ts_healthy_noPrEP_seek_care <- ts_healthy_seek_care - ts_healthy_PrEP_seek_care
  
  # ---------------------------------------------------------------------------#
  #  5: Biologicals start & compliance ######
  # ---------------------------------------------------------------------------#
  pNoStart <- 1 - pSeek_exposure * pStart_exposure   # ### OPT: scalar, computed once
  
  ## Exposures — no PrEP
  ts_exp_noPrEP_start       <- matrix(rbinom(N * horizon, as.vector(ts_exp_noPrEP_seek_care), pStart_exposure), N, horizon)
  ts_exp_noPrEP_nostartPEP  <- matrix(rbinom(N * horizon, as.vector(ts_exp_noPrEP),           pNoStart),        N, horizon)
  ts_exp_noPrEP_second_dose <- matrix(rbinom(N * horizon, as.vector(ts_exp_noPrEP_start),      pCompliance_exp), N, horizon)
  ts_exp_noPrEP_completePEP <- matrix(rbinom(N * horizon, as.vector(ts_exp_noPrEP_second_dose),pCompliance_exp), N, horizon)
  ts_exp_noPrEP_incompletePEP <- ts_exp_noPrEP_start - ts_exp_noPrEP_completePEP
  
  ## Exposures — PrEP (only 2 doses needed)
  ts_exp_PrEP_start        <- matrix(rbinom(N * horizon, as.vector(ts_exp_PrEP_seek_care), pStart_exposure), N, horizon)
  ts_exp_PrEP_no_start     <- ts_exp_PrEP - ts_exp_PrEP_start
  ts_exp_PrEP_complete     <- matrix(rbinom(N * horizon, as.vector(ts_exp_PrEP_start),     pCompliance_exp), N, horizon)
  ts_exp_PrEP_incompletePEP <- ts_exp_PrEP_start - ts_exp_PrEP_complete
  
  ts_exp_start    <- ts_exp_noPrEP_start      + ts_exp_PrEP_start
  ts_exp_no_start <- ts_exp_noPrEP_nostartPEP + ts_exp_PrEP_no_start
  ts_exp_incomplete <- ts_exp_noPrEP_incompletePEP + ts_exp_PrEP_incompletePEP
  ts_exp_complete   <- ts_exp_noPrEP_completePEP   + ts_exp_PrEP_complete
  
  ## Healthy bites — no PrEP
  ts_healthy_noPrEP_start       <- matrix(rbinom(N * horizon, as.vector(ts_healthy_noPrEP_seek_care),   pStart_healthy),    N, horizon)
  ts_healthy_noPrEP_second_dose <- matrix(rbinom(N * horizon, as.vector(ts_healthy_noPrEP_start),       pCompliance_healthy),N, horizon)
  ts_healthy_noPrEP_complete    <- matrix(rbinom(N * horizon, as.vector(ts_healthy_noPrEP_second_dose), pCompliance_healthy),N, horizon)
  ts_healthy_noPrEP_incomplete  <- ts_healthy_noPrEP_start - ts_healthy_noPrEP_complete
  
  ## Healthy bites — PrEP (2 doses)
  ts_healthy_PrEP_start    <- matrix(rbinom(N * horizon, as.vector(ts_healthy_PrEP_seek_care), pStart_healthy),    N, horizon)
  ts_healthy_PrEP_complete <- matrix(rbinom(N * horizon, as.vector(ts_healthy_PrEP_start),     pCompliance_healthy),N, horizon)
  ts_healthy_PrEP_incomplete <- ts_healthy_PrEP_start - ts_healthy_PrEP_complete
  
  ts_healthy_start    <- ts_healthy_PrEP_start    + ts_healthy_noPrEP_start
  ts_healthy_complete <- ts_healthy_PrEP_complete + ts_healthy_noPrEP_complete
  ts_healthy_incomplete <- ts_healthy_PrEP_incomplete + ts_healthy_noPrEP_incomplete
  
  ## RIG
  tmp    <- if (ibcm == "no") ts_exp_noPrEP_start + ts_healthy_noPrEP_start else ts_exp_noPrEP_start
  ts_RIG <- matrix(rbinom(N * horizon, as.vector(tmp), RIG_cov), N, horizon)
  
  # ---------------------------------------------------------------------------#
  # 6: IBCM  (Not used in this model) ######
  # ---------------------------------------------------------------------------#
  bites_by_each_rabid_biting_dog <- rabies_results$bites_by_dog
  
  #: biters_sought_care — avoid rep() inside loop by pre-expanding once
  #          per cell; still O(N*horizon) iterations but the body is leaner.
  n_cells <- N * horizon
  biters_sought_care <- vector("list", n_cells)
  exp_seek_vec       <- as.integer(ts_exp_seek_care)
  
  for (x in seq_len(n_cells)) {
    b <- bites_by_each_rabid_biting_dog[[x]]
    if (length(b) == 0L) {
      biters_sought_care[[x]] <- integer(0L)
    } else {
      # rep(dog_id, n_bites) then sample — keep original semantics
      biters_sought_care[[x]] <- sample(
        rep.int(seq_along(b), b),
        exp_seek_vec[x]
      )
    }
  }
  
  ts_rabid_bites_investigated <- matrix(
    rbinom(N * horizon, as.vector(ts_exp_seek_care), pInvestigate), N, horizon)
  
  #: investigated biters — same loop structure, slightly tightened
  biters_sought_care_investigated <- vector("list", n_cells)
  n_investigate_vec <- as.integer(ts_rabid_bites_investigated)
  for (x in seq_len(n_cells)) {
    b <- biters_sought_care[[x]]
    ni <- n_investigate_vec[x]
    if (length(b) == 0L || is.na(ni) || ni <= 0L) {
      biters_sought_care_investigated[[x]] <- integer(0L)
    } else {
      biters_sought_care_investigated[[x]] <- b[sample.int(length(b), min(ni, length(b)))]
    }
  }
  
  ts_rabid_biting_investigated <- matrix(
    vapply(biters_sought_care_investigated, function(x) length(unique(x)), 1L), nrow = N)
  
  ts_rabid_biting_found    <- matrix(rbinom(N * horizon, as.vector(ts_rabid_biting_investigated), pFound),    N, horizon)
  ts_rabid_biting_testable <- matrix(rbinom(N * horizon, as.vector(ts_rabid_biting_found),        pTestable), N, horizon)
  
  ts_healthy_biting_investigated <- matrix(
    rbinom(N * horizon, as.vector(ts_healthy_seek_care), pFS), N, horizon)
  
  # ---------------------------------------------------------------------------#
  # 7: Outcomes 
  ## 7A: deaths ########
  # ---------------------------------------------------------------------------#
  ts_exp_PrEP_protected        <- matrix(rbinom(N * horizon, as.vector(ts_exp_PrEP), PrEP_effectiveness), N, horizon)
  ts_exp_PrEP_not_protected    <- ts_exp_PrEP - ts_exp_PrEP_protected
  
  ts_exp_PrEP_protected_no_start <- matrix(
    rbinom(N * horizon, as.vector(ts_exp_PrEP_protected), pNoStart), N, horizon)
  ts_exp_PrEP_protected_got_PEP  <- ts_exp_PrEP_protected - ts_exp_PrEP_protected_no_start
  
  ts_exp_PrEP_not_protected_no_start <- matrix(
    rbinom(N * horizon, as.vector(ts_exp_PrEP_not_protected), pNoStart), N, horizon)
  
  ts_deaths_no_PEP <- matrix(
    rbinom(N * horizon,
           pmax(0L, as.vector(ts_exp_noPrEP_nostartPEP + ts_exp_PrEP_not_protected_no_start)),
           pDeath), N, horizon)
  
  deaths_incomplete_PEP <- matrix(
    rbinom(N * horizon, as.vector(ts_exp_noPrEP_incompletePEP), 1 - pPrevent_incomplete), N, horizon)
  deaths_complete_PEP   <- matrix(
    rbinom(N * horizon, as.vector(ts_exp_noPrEP_completePEP),   1 - pPrevent_complete),   N, horizon)
  
  ts_deaths <- ts_deaths_no_PEP + deaths_incomplete_PEP + deaths_complete_PEP
  
  ## 7B: Deaths averted  ########
  ts_deaths_averted_PrEP1 <- matrix(
    rbinom(N * horizon, as.vector(ts_exp_PrEP_protected_no_start), pDeath), N, horizon)
  
  tmp2 <- matrix(
    rbinom(N * horizon, as.vector(ts_exp_PrEP_incompletePEP), 1 - pPrevent_incomplete), N, horizon)
  ts_deaths_averted_PrEP2 <- matrix(
    rbinom(N * horizon, as.vector(tmp2), pDeath), N, horizon)
  ts_deaths_averted_PrEP  <- ts_deaths_averted_PrEP1 + ts_deaths_averted_PrEP2
  
  deaths_averted_PEP_complete   <- matrix(
    rbinom(N * horizon, as.vector(ts_exp_complete),   pPrevent_complete   * pDeath), N, horizon)
  deaths_averted_PEP_incomplete <- matrix(
    rbinom(N * horizon, as.vector(ts_exp_incomplete), pPrevent_incomplete * pDeath), N, horizon)
  
  ts_deaths_averted_PEP    <- deaths_averted_PEP_complete + deaths_averted_PEP_incomplete
  ts_deaths_averted_PEPPrEP <- ts_deaths_averted_PrEP + deaths_averted_PEP_complete + deaths_averted_PEP_incomplete
  
  # ---------------------------------------------------------------------------#
  # MDV counterfactual — ### OPT: reuse dog_pop; only recompute what changes
  # ---------------------------------------------------------------------------#

  vax_cov_no_MDV <-
    calc_vax_coverage(0.05, 0.05, NULL, vaccinate_owned_dog_cost, owned_dogs,
                      total_horizon, discount, years_to_target, dog_burnin) * owned_prop +
    calc_vax_coverage(0,    0,    NULL, vaccinate_unowned_dog_cost, unowned_dogs,
                      total_horizon, discount, years_to_target, dog_burnin) * unowned_prop
  
  exposures_no_MDV <- strip_burnin(
    predict_dograbies_split(
      N = N, horizon = total_horizon, vax_cov = vax_cov_no_MDV, dog_pop = dog_pop,
      rabies_inc = rabies_inc, mu = mu, k = k, seed = seed,
      pop_serengeti = 1e5, split_by = "mean")$ts_exposures)
  
  ts_deaths_averted_MDV  <- (exposures_no_MDV - rabies_results$ts_exposures) * pDeath
  ts_deaths_averted_MDV2 <- matrix(
    rbinom(N * horizon, pmax(0L, as.vector(exposures_no_MDV - rabies_results$ts_exposures)), pDeath),
    N, horizon)
  
  expected_deaths_no_intervention <- matrix(
    rbinom(N * horizon, as.vector(exposures_no_MDV), pDeath), N, horizon)
  
  ts_deaths_averted2 <- ts_deaths_averted_PEP + ts_deaths_averted_PrEP + ts_deaths_averted_MDV
  ts_deaths_averted  <- expected_deaths_no_intervention - ts_deaths
  
  # # discount deaths averted/ lives saved -- not discounting deaths averted but will discount QALYs gained outside main function
  # disc <- (1 + discount)^(-(0:(horizon - 1)))
  # ts_deaths_averted_PEP   <- sweep(ts_deaths_averted_PEP, 2, disc, `*`)
  # ts_deaths_averted_PrEP  <- sweep(ts_deaths_averted_PrEP, 2, disc, `*`)
  # ts_deaths_averted_MDV   <- sweep(ts_deaths_averted_MDV, 2, disc, `*`)
  # ts_deaths_averted       <- sweep(ts_deaths_averted, 2, disc, `*`)
  # ts_deaths_averted2 <- ts_deaths_averted_PEP + ts_deaths_averted_PrEP + ts_deaths_averted_MDV
  
  
  # ---------------------------------------------------------------------------#
  #  8: Economics  ########
  # ---------------------------------------------------------------------------#
  ts_complete_PEP   <- ts_exp_complete  + ts_healthy_complete
  ts_incomplete_PEP <- ts_exp_incomplete + ts_healthy_incomplete
  
  PrEP_vials      <- ts_pop_PrEP * PrEP_vials_per_pt
  ts_PEP_vials    <- (ts_exp_start + ts_healthy_start) * PEP_vials_per_pt
  ts_vaccine_vials <- ts_PEP_vials + PrEP_vials
  
  # discount costs
  disc <- (1 + discount)^(-(0:(horizon - 1)))
  ts_cost_PrEP         <- sweep(PrEP_vials    * human_vaccine_cost_per_vial_PrEP, 2, disc, `*`)
  ts_RIG_cost_per_year <- sweep(ts_RIG        * RIG_cost,                    2, disc, `*`)
  ts_cost_PEP_per_year <- sweep(ts_PEP_vials  * human_vaccine_cost_per_vial_PEP, 2, disc, `*`)
  ts_MDV_campaign_cost <- sweep(MDV_campaign_cost,                           2, disc, `*`)
  ts_cost_per_year     <- ts_MDV_campaign_cost + ts_cost_PrEP + ts_cost_PEP_per_year + ts_RIG_cost_per_year
  
  # ---------------------------------------------------------------------------#
  # 9: Collate & return  #######
  # ---------------------------------------------------------------------------#
  
  # All `ts_` objects from the local environment
  my_list      <- ls(pattern = "^ts_")
  out_matrices <- mget(my_list)
  
  # ts_ elements from rabies_results (already stripped)
  rabies_ts <- rabies_results[grep("^ts_", names(rabies_results))]
  dogs_ts   <- vax_results[grep("^ts_", names(vax_results))]
  
  # CHANGE: burn-in already stripped, so NO second stripping here
  out_matrices <- c(out_matrices, rabies_ts, dogs_ts)
  
  return(out_matrices)
}


# ── Single-session example ───────────────────────────────────────────────────

load_rabies_models()   # <-- call once; cached for the whole session

tmp <- decision_tree_PrEP(
  N = 5, pop=35000000, HDR=c(17,18), unowned_prop=0.58, horizon=5, discount=0.03, mu = 0.38, k = 0.72, 
 #pBite_healthy = 0.2, pSeek_healthy = 0.2, 
  pStart_healthy = 0.9, pCompliance_healthy = 0.5, bpi =15.3,
  pSeek_exposure = 0.95, pStart_exposure = 0.9, pCompliance_exp = 0.96, pDeath = 0.16, 
  pPrevent_complete = 0.99, pPrevent_incomplete = 0.98, rabies_inc= c(0.0075,0.0125),
  mdv_unowned_budget = NULL,  mdv_owned_budget = NULL, vaccinate_owned_dog_cost = c(2,4), 
  vaccinate_unowned_dog_cost = c(2,4), base_vax_cov_owned = 0.4, target_vax_cov_owned = 0.4, 
  base_vax_cov_unowned = 0.4, target_vax_cov_unowned = 0.4, years_to_target =3,
  pInvestigate = 0.5, pFound = 0.4, pTestable = 0.2, pFS=0.05, RIG_cov =0.4,
  init_PrEP_cov=0.14, # If zero, means no PrEP. If >0, this is the coverage achieved through school campaigns (relative to entire population)
  PrEP_effectiveness = 1, PrEP_vials_per_pt = 0.6, PEP_vials_per_pt = 0.6,
  human_vaccine_cost_per_vial_PEP=5, human_vaccine_cost_per_vial_PrEP=5, RIG_cost=70, birth_rate = 0.01, routine_vax_cov = 0.7,
  seed = 123, ibcm = "no", dog_burnin = 3
 )

tmp$ts_exposures
tmp$ts_rabid_dogs
tmp$ts_cum_pop_PrEP

