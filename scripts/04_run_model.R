
# Run multiple scenarios
#  read parameters   ####
source("./scripts/02_decision_tree_incl_PrEP.R")

##  Common (fixed) arguments
common_args1 <- list(
  N       = 3000,
  pop     = 35000000,
  HDR     = c(15, 16),
  horizon = 10, 
  dog_burnin = 1, # running for a year before the model timeline "begins" to stabilize dog dynamics
  birth_rate = 0.012,
  routine_vax_cov = 0.778,
  vaccinate_unowned_dog_cost = c(4,5), # 380 rupees
  vaccinate_owned_dog_cost = c(0.5, 1), # 80 rupees, 50 no certificate 
  mdv_unowned_budget = NULL,
  mdv_owned_budget = NULL,
  rabies_inc = c(0.0075, 0.0125),
  discount = 0.03,
  unowned_prop = 0.633,
  mu = 0.38,
  k = 0.72,
  pDeath = 0.17, 
  pPrevent_complete = 0.999,
  pPrevent_incomplete = 0.986,
  bpi = 15.3, # dogs only (8.3); including cats 15.3 (per 1000)
  ibcm = "no",
  PrEP_effectiveness = 0
)


##  Scenario-specific parameters
scenarios1 <- read.csv("./data/scenario_parameters_Kerala_general.csv") %>%
  dplyr::mutate(
    across(
      c(
        pStart_healthy, pCompliance_healthy,
        pSeek_exposure, pStart_exposure, pCompliance_exp,
        base_vax_cov_unowned, base_vax_cov_owned,
        target_vax_cov_unowned, target_vax_cov_owned,
        pInvestigate, pFound, pTestable, pFS, RIG_cov,
        init_PrEP_cov, PrEP_vials_per_pt, PEP_vials_per_pt,
        human_vaccine_cost_per_vial_PrEP,
        human_vaccine_cost_per_vial_PEP, RIG_cost
      ),
      as.numeric
    ),
    ibcm = as.character(ibcm)
  )


# =================================================================#
# run one scenario ##########
# =================================================================#

run_scenario <- function(scenario_row, common_args, seed = 100L, discount_rate = 0.03) {
  
  load_rabies_models()
  
  scenario_list <- as.list(scenario_row[setdiff(names(scenario_row), "scenario")])
  scenario_args <- do.call(
    purrr::list_modify,
    c(list(common_args), scenario_list, list(seed = seed))
  )
  
  scen <- do.call(decision_tree_PrEP, scenario_args)
  
  # deaths averted are left undiscounted
  scen$ts_cost_per_life_saved <- scen$ts_cost_per_year / scen$ts_deaths_averted
  
  # undiscounted QALYs gained
  scen$ts_QALYs_gained <- scen$ts_deaths_averted * 62
  
  # discount applied across years (columns)
  n_years <- ncol(scen$ts_QALYs_gained)
  discount_factors <- 1 / ((1 + discount_rate) ^ (0:(n_years - 1)))
  
  scen$ts_QALYs_gained_discounted <- sweep(
    scen$ts_QALYs_gained,
    MARGIN = 2,
    STATS = discount_factors,
    FUN = "*"
  )
  
  # cost per discounted QALY gained
  scen$ts_cost_per_QALY_gained <- scen$ts_cost_per_year / scen$ts_QALYs_gained_discounted
  
  scen
}

# =================================================================#
# parallel execution across scenarios
# =================================================================#

plan(multisession, workers = 8)

# Kerala_general_parallel <- scenarios1 %>%
#   dplyr::mutate(
#     output = future_map2(
#       .x = split(select(., -scenario), seq_len(nrow(.))),
#       .y = 100L + seq_len(nrow(.)),
#       .f = ~ run_scenario(
#         scenario_row = .x,
#         common_args = common_args1,
#         seed = .y
#       ),
#       .options = furrr_options(seed = TRUE)
#     )
#   ) %>%
#   { set_names(.$output, .$scenario) }


# Save
# saveRDS(Kerala_general_parallel, file = "output/Kerala_general.rds")


# =================================================================#
# Post run model######

Kerala_general <- readRDS("output/Kerala_general.rds")

Kerala_general_scenarios <- names(Kerala_general)

vars_to_summarise <- c(
  "ts_exposures",
  "ts_deaths",
  "ts_deaths_averted",
  "ts_deaths_averted_PEP",
  "ts_deaths_averted_PrEP",
  "ts_deaths_averted_MDV",
  "ts_cost_per_year",
  "ts_cost_per_life_saved",
  "ts_QALYs_gained",
  "ts_QALYs_gained_discounted",
  "ts_cost_per_QALY_gained",
  "ts_cost_PEP_per_year",
  "ts_cost_PrEP",
  "ts_RIG_cost_per_year",
  "ts_MDV_campaign_cost"
)

output_table <- summarise_variables_across_scenarios(
  results_list   = Kerala_general,
  variables      = vars_to_summarise,
  scenario_names = names(Kerala_general)
)


output_table$scenario
# order output table to match scenarios as listed in Table 2
scenario_order <- c("SQ" , "SQ_addIndian_PrEP", "WHO_PEP_1ml_IM_PrEP", "WHO_PEP_and_PrEP", 
                    "WHO_PEP_and_PrEP_reducedMDVUnowned", "WHO_PEPandPrEP_MDVscaleupUnowned", "SQ_mod_MDVUnowned",
                    "SQ_MDVscaleupUnowned", "SQ_MDVscaleupOwned", "SQ_MDVscaleupAll", "WHO_PEP_and_SQ_MDV",
                    "WHO_PEP_MDVscaleupAll")


output_table <- output_table %>%
  dplyr::mutate(scenario = factor(scenario, levels = scenario_order),
                # ICERs now not needed as we move to a NHB framework -- but leaving in for the CE plane (supplementary)
                
                # Incremental Cost and Incremental Effect (relative to SQ)
                inc_cost = ts_cost_per_year_Median - (output_table$ts_cost_per_year_Median[output_table$scenario == "SQ"]),
                inc_effect = ts_QALYs_gained_discounted_Median - (output_table$ts_QALYs_gained_discounted_Median[output_table$scenario == "SQ"]),
                
                # ICER = ΔCost / ΔEffect
                #ICER = ifelse(scenario == "SQ", NA_real_, inc_cost / inc_effect),
                
                # RIG is part of PEP!
                ts_cost_PEP_per_year_Median = ts_cost_PEP_per_year_Median + ts_RIG_cost_per_year_Median
                
  ) %>%
  dplyr::select(-ts_RIG_cost_per_year_Median) %>%
  arrange(scenario)



# Helper function to format pretty table

fmt_median_ci <- function(med, ll, ul, digits = 0, scale = 1) {
  med_f <- scales::comma(round(med / scale, digits))
  ll_f  <- scales::comma(round(ll  / scale, digits))
  ul_f  <- scales::comma(round(ul  / scale, digits))
  
  paste0(med_f, " (", ll_f, "–", ul_f, ")")
}




scenario_labels <- c(
  SQ = "A. SQ", 
  SQ_addIndian_PrEP = "B. SQ + PrEP (India regimen)",
  WHO_PEP_1ml_IM_PrEP = "C. WHO PEP + 1-shot PrEP",
  WHO_PEP_and_PrEP = "D. WHO PEP & PrEP",
  WHO_PEP_and_PrEP_reducedMDVUnowned = "E. D + drop in MDV",
  WHO_PEPandPrEP_MDVscaleupUnowned = "F. D + MDV scale-up in unowned",
  SQ_mod_MDVUnowned = "G. SQ + moderate MDV in unowned",
  SQ_MDVscaleupUnowned = "H. SQ + MDV scale-up in unowned",
  SQ_MDVscaleupOwned = "I. SQ + intense MDV in owned",
  SQ_MDVscaleupAll = "J. SQ + MDV scale-up in all", 
  WHO_PEP_and_SQ_MDV = "K. WHO PEP + SQ MDV",
  WHO_PEP_MDVscaleupAll = "L. WHO PEP + MDV scale-up all"
)

# relabel helper
add_scenario_labels <- function(df, label_map = scenario_labels) {
  df %>%
    dplyr::mutate(
      scenario = as.character(scenario),
      scenario_label = recode(scenario, !!!label_map, .default = scenario)
    )
}

# add labels
output_table <- add_scenario_labels(output_table)


# Rewrite pretty_table using LL / Median / UL triplets
pretty_table <- output_table %>%
  transmute(
    Scenario = scenario_label,
    
    Deaths = fmt_median_ci(
      ts_deaths_Median, ts_deaths_LL, ts_deaths_UL
    ),
    
    `Deaths averted` = fmt_median_ci(
      round(ts_deaths_averted_Median, 0),
      round(ts_deaths_averted_LL, 0),
      round(ts_deaths_averted_UL, 0)
    ),
    
    `Total QALYs gained` = fmt_median_ci(
      round(ts_QALYs_gained_Median, 0),
      round(ts_QALYs_gained_LL, 0),
      round(ts_QALYs_gained_UL, 0)
    ),
    
    `Total QALYs gained (discounted)` = fmt_median_ci(
      round(ts_QALYs_gained_discounted_Median, 0),
      round(ts_QALYs_gained_discounted_LL, 0),
      round(ts_QALYs_gained_discounted_UL, 0)
    ),
    
    `Total costs ($M)` = fmt_median_ci(
      ts_cost_per_year_Median, ts_cost_per_year_LL, ts_cost_per_year_UL,
      digits = 2, scale = 1e6
    ),
    
    `Cost per QALY gained ($)` = fmt_median_ci(
      ts_cost_per_QALY_gained_Median,
      ts_cost_per_QALY_gained_LL,
      ts_cost_per_QALY_gained_UL
    )#,
    
    #ICER = round(ICER, 0)
  )


source("./scripts/07_nhb_analysis.R")
merge_table1
pretty_table
pretty_table <- pretty_table %>%
  dplyr::select(-c("Total QALYs gained"#, "ICER"
  )) %>%
  left_join(., merge_table1, by = "Scenario")



# Export
## table 2 ####
write.csv(pretty_table, "./output/Table2.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(output_table, "./output/Kerala_general_out.csv")



# Run scenario where Prep_effectiveness ==1 ##########
## Run model #####
# filter out scenarios without PrEP to save time
scenarios2 <- scenarios1 %>%
  dplyr::filter(stringr::str_detect(scenario, "PrEP"))

# update PrEP_effectiveness to 1
common_args2 <- common_args1
common_args2$PrEP_effectiveness <- 1

## Run parallel
# Kerala_general_prepEff1 <- scenarios2 %>%
#   dplyr::mutate(
#     output = future_map2(
#       .x = split(select(., -scenario), seq_len(nrow(.))),
#       .y = 3000L + seq_len(nrow(.)),
#       .f = ~ run_scenario(
#         scenario_row = .x,
#         common_args = common_args2,
#         seed = .y
#       ),
#       .options = furrr_options(seed = TRUE)
#     )
#   ) %>%
#   { set_names(.$output, .$scenario) }



# Save
# saveRDS(Kerala_general_prepEff1, file = "output/Kerala_general_prepEff1.rds")


## Summarize model #####
Kerala_general_prepEff1 <- readRDS("output/Kerala_general_prepEff1.rds")

scenario_order2 <- c("SQ_addIndian_PrEP", "WHO_PEP_1ml_IM_PrEP", "WHO_PEP_and_PrEP", 
                     "WHO_PEP_and_PrEP_reducedMDVUnowned", "WHO_PEPandPrEP_MDVscaleupUnowned")


prep_eff1_table <- summarise_variables_across_scenarios(
  results_list   = Kerala_general_prepEff1,
  variables      = vars_to_summarise,
  scenario_names = names(Kerala_general_prepEff1)
)

output_table2 <- prep_eff1_table %>%
  dplyr::mutate(scenario = factor(scenario, levels = scenario_order2),
                # Incremental Cost and Incremental Effect (relative to SQ)
                inc_cost = ts_cost_per_year_Median - (output_table$ts_cost_per_year_Median[output_table$scenario == "SQ"]),
                inc_effect = ts_deaths_averted_Median - (output_table$ts_deaths_averted_Median[output_table$scenario == "SQ"]),
                # ICER = ΔCost / ΔEffect
                ICER = ifelse(scenario == "SQ", NA_real_, inc_cost / inc_effect),
                # RIG is part of PEP!
                ts_cost_PEP_per_year_Median = ts_cost_PEP_per_year_Median + ts_RIG_cost_per_year_Median
  ) %>%
  dplyr::select(-ts_RIG_cost_per_year_Median) %>%
  arrange(scenario)


write.csv(output_table2, "./output/Kerala_general_prepEff1_summarized.csv")







