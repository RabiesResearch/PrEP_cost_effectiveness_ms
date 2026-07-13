

# source model
source("./scripts/02_decision_tree_incl_PrEP.R")


## Owned vs unowned & owned vaccination coverage ########
## Owned vs unowned, owned vaccination coverage & care-seeking sensitivity ########

# Known counts & uncertainty (govt livestock census)
owned_dogs_census   <- 836000
unowned_dogs_census <- 289000
unowned_multipliers <- c(1, 2, 3, 4, 5, 6, 7)                   # up to 7× unowned
owned_cov_grid      <- c(0.4,0.5)      # 30–50% # seq(0.30, 0.50, by = 0.10) 
pSeek_grid          <- c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1)        # care-seeking probabilities = c(0.5 -1) 


# Common (SQ-like) arguments
common_args <- list(
  pop = 35e6,
  N = 1000,
  HDR = c(30, 32),          # placeholder; overwritten per grid
  horizon = 10,
  dog_burnin = 1,
  discount = 0,
  mu = 0.38, k = 0.72,
  # pBite_healthy = 0.18,
  # pSeek_healthy = 0.448195,
  bpi = 15.3, # per 1k persons (26.9 vs 15.3)
  pStart_healthy = 0.9684211,
  pCompliance_healthy = 0.9654,
  pSeek_exposure = 0.95,      # overwritten per grid
  pStart_exposure = 0.9684211,
  pCompliance_exp = 0.9654,
  pDeath = 0.17,
  pPrevent_complete = 0.999,
  pPrevent_incomplete = 0.986,
  rabies_inc = c(0.0075, 0.0125),
  mdv_unowned_budget = NULL,
  mdv_owned_budget = NULL,
  vaccinate_owned_dog_cost = c(1, 2),
  vaccinate_unowned_dog_cost = c(2, 3),
  base_vax_cov_unowned = 0.03,
  target_vax_cov_unowned = 0.03,
  base_vax_cov_owned = 0.50,
  target_vax_cov_owned = 0.50,
  pInvestigate = 0.5, pFound = 0.4, pTestable = 0.2, pFS = 0.05,
  RIG_cov = 0.38,
  init_PrEP_cov = 0,
  PrEP_effectiveness = 0,
  PrEP_vials_per_pt = 0.66, PEP_vials_per_pt = 0.6,
  human_vaccine_cost_per_vial_PrEP = 5,
  human_vaccine_cost_per_vial_PEP = 5,
  RIG_cost = 10,
  birth_rate = 0.012,
  routine_vax_cov = 0.78,
  ibcm = "no"
)

# Helper: derive HDR & unowned_prop from counts
derive_dog_params <- function(pop, owned, unowned_est, m) {
  unowned <- m * unowned_est
  total_dogs <- owned + unowned
  hdr <- pop / total_dogs
  unowned_prop <- unowned / total_dogs
  
  list(
    HDR = c(hdr, hdr),
    unowned_prop = unowned_prop,
    dogs_total = total_dogs,
    owned = owned,
    unowned = unowned
  )
}


# Build 3D parameter grid
grid <- tidyr::crossing(
  unowned_mult = unowned_multipliers,
  owned_cov    = owned_cov_grid,
  pSeek_exp    = pSeek_grid
)


# 1) Pre-compute deterministic parts (sequential)
grid2 <- grid %>%
  dplyr::mutate(
    dog_params = map(unowned_mult, 
                     ~ derive_dog_params(
                       pop = common_args$pop,
                       owned = owned_dogs_census,
                       unowned_est = unowned_dogs_census,
                       m = .x
                       )),
    HDR          = map(dog_params, "HDR"),
    unowned_prop = map_dbl(dog_params, "unowned_prop"),
    dogs_total   = map_dbl(dog_params, "dogs_total"),
    
    
    # build args for this cell (a single list of args per row)
    args_base = pmap(
      list(HDR, unowned_prop, owned_cov, pSeek_exp),
      ~ purrr::list_modify(
        common_args,
        HDR = ..1,
        unowned_prop = ..2,
        base_vax_cov_owned   = ..3,
        target_vax_cov_owned = ..3,
        pSeek_exposure       = ..4
      )
    )
  )

# 2) A helper that runs ONE grid cell end-to-end

run_one_cell <- function(args_base) {
  load_rabies_models()
  out <- do.call(decision_tree_PrEP, args_base)
  
  summarise_vec <- function(x) {
    c(LL     = unname(quantile(x, 0.025, na.rm = TRUE)),  # unname() drops "2.5%"
      Median = median(x, na.rm = TRUE),
      UL     = unname(quantile(x, 0.975, na.rm = TRUE)))
  }
  
  deaths_summary    <- summarise_vec(rowSums(out$ts_deaths))
  bites_summary     <- summarise_vec(rowSums(out$ts_recorded_bites))
  exposures_summary <- summarise_vec(rowSums(out$ts_exposures))
  
  list(
    deaths_LL        = deaths_summary[["LL"]],
    deaths_Median    = deaths_summary[["Median"]],
    deaths_UL        = deaths_summary[["UL"]],
    bites_LL         = bites_summary[["LL"]],
    bites_Median     = bites_summary[["Median"]],
    bites_UL         = bites_summary[["UL"]],
    exposures_LL     = exposures_summary[["LL"]],
    exposures_Median = exposures_summary[["Median"]],
    exposures_UL     = exposures_summary[["UL"]]
  )
}



# Run grid

plan(multisession, workers = 8)  # choose workers

# 3) Parallel step: run each grid row
grid3 <- grid2 %>%
  dplyr::mutate(
    summary = future_map(
      args_base,
      run_one_cell,
      .options = furrr_options(seed = TRUE)  # reproducible RNG across workers
    )
  )


# 4) Unpack and finish (sequential)
results <- grid3 %>%
  tidyr::unnest_wider(summary) %>%
  dplyr::mutate(
    exposure_pct_LL     = (exposures_LL     / bites_LL)     * 100,
    exposure_pct_Median = (exposures_Median / bites_Median) * 100,
    exposure_pct_UL     = (exposures_UL     / bites_UL)     * 100
  ) %>%
  dplyr::select(
    unowned_mult, HDR, unowned_prop, owned_cov, pSeek_exp,
    deaths_LL, deaths_Median, deaths_UL,
    bites_LL, bites_Median, bites_UL,
    exposures_LL, exposures_Median, exposures_UL,
    exposure_pct_LL, exposure_pct_Median, exposure_pct_UL
  ) %>%
  dplyr::mutate(
    HDR = map_dbl(HDR, mean),
    unowned_prop = round(unowned_prop, 3),
    owned_vax_cov = owned_cov,
    median_deaths = deaths_Median
  ) %>%
  arrange(pSeek_exp, unowned_mult, owned_vax_cov)



# Preview results
results %>% print(n = 5)

plot2A_data <- results 
# 
# 1) Total seek care/ start PEP (data- 941k); 2) exposures/ all seek care (0.5-3%); 3) deaths (5-26, but may be under reported)

# Save
# write.csv(results, "output/model_calibration2.csv")

# visualization #####
#Supp Figure 1 ####

plot2A_data <- read.csv("./output/model_calibration.csv")


# targets from observed data
# recorded deaths 2012-2025
krl_deaths <- c(13,11, 10,10,5,8,9,8,5,11,27,25,26,33)
mean(krl_deaths)
median(krl_deaths)
# last 4 years
mean(tail(krl_deaths, 4))
median(tail(krl_deaths, 4))


target_deaths <- 15 # all years
target_deaths <- (28/0.7) *10 # last 4 years (Optionally = 28/0.7 as 70% of deaths are reported/ neurological *attacks-- find citation)
# Multiply deaths by 10 years 
target_exposure_pct <- 0.5  # note: this is in percent units (e.g., 0.5-3 where 3 = 3%)

plot2A_scored <- plot2A_data %>%
  dplyr::mutate(
    # squared deviations from targets
    deaths_se   = ((deaths_Median - target_deaths)/target_deaths)^2,
    expos_se    = (exposure_pct_Median - target_exposure_pct)^2,
    
    # optional: rescale each SE to 0–1 so they contribute equally in magnitude
    deaths_se_s = rescale(deaths_se, to = c(0, 1)),
    expos_se_s  = rescale(expos_se,  to = c(0, 1)),
    
    # equal-weight composite loss (lower = better)
    loss = 0.5 * deaths_se_s + 0.5 * expos_se_s
  )

best <- plot2A_scored %>% slice_min(loss, n = 1) %>% 
  dplyr::select(unowned_mult, pSeek_exp, deaths_Median, exposure_pct_Median, loss)

Supp_Fig1 <- ggplot(
  plot2A_scored,
  aes(x = factor(unowned_mult),
      y = factor(pSeek_exp),
      fill = loss)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_point(
    data = best,
    aes(x = factor(unowned_mult), y = factor(pSeek_exp)),
    inherit.aes = FALSE,
    shape = 21, size = 6, fill = "#FFD92F", color = "black", stroke = 1.3
  ) +
  # scale_fill_gradient( 
  #   low = "#e8e4f3", 
  #   high = "#5e3c78", 
  #   limits = c(0, 0.5),
  #   name = "Composite loss\n(lower = better)" 
  #   ) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    limits = c(0, 1),
    oob = scales::squish,
    name = "Composite loss\n(lower = better)"
  ) +
  labs(
    x = "Unowned dog population multiplier",
    y = "Probability of seeking care"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

Supp_Fig1

#Export Supp_Fig1  ####
pdf("./output/manuscript_figures/SupplementaryFigure1.pdf", width = 6, height = 4.5, useDingbats = FALSE)
print(Supp_Fig1)
dev.off()



