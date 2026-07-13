

library(tidyverse)
library(scales)

Kerala_general <- readRDS("output/Kerala_general.rds")

WTP_list <- c(0, 100, 200, 300, 400, 500, 1000, 2034, 2535, 3092, 4000, 5000)

# Selected scenarios
# Rationale: SQ: current practice, WHO_PEP_and_SQ_MDV: Only cost saving strategy, 
## SQ_addIndian_PrEP: We want to see if prep is CE. This is easist to implement based on national policies
## SQ_MDVscaleupOwned & WHO_PEP_MDVscaleupAll : Other scnearios in the CE efficiency frontier ( CE plane)

selected_scenarios <- c(
  "SQ",
  "WHO_PEP_and_SQ_MDV",
  "SQ_addIndian_PrEP",
  "SQ_MDVscaleupOwned",
  "WHO_PEP_MDVscaleupAll"
)

# Dictionary for fancier labels
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
    mutate(
      scenario_label = recode(scenario, !!!label_map, .default = scenario)
    )
}



# Main function: NHB + CEAC #######
calc_nhb_ceac <- function(Kerala_general, WTP_list,
                          qalys_name = "ts_QALYs_gained_discounted",
                          costs_name = "ts_cost_per_year",
                          scenarios = NULL) {
  
  scenario_names <- if (is.null(scenarios)) names(Kerala_general) else scenarios
  
  nhb_df <- do.call(rbind, lapply(scenario_names, function(scn) {
    qalys <- rowSums(Kerala_general[[scn]][[qalys_name]], na.rm = TRUE)
    costs <- rowSums(Kerala_general[[scn]][[costs_name]], na.rm = TRUE)
    
    do.call(rbind, lapply(WTP_list, function(wtp) {
      nhb <- if (wtp == 0) -costs else qalys - costs / wtp
      
      data.frame(
        scenario = scn,
        iteration = seq_along(nhb),
        WTP = wtp,
        NHB = nhb
      )
    }))
  }))
  
  mean_nhb <- nhb_df %>%
    group_by(scenario, WTP) %>%
    summarise(mean_NHB = mean(NHB, na.rm = TRUE), .groups = "drop")
  
  winners <- nhb_df %>%
    group_by(WTP, iteration) %>%
    mutate(max_NHB = max(NHB, na.rm = TRUE)) %>%
    filter(NHB == max_NHB) %>%
    mutate(weight = 1 / n()) %>%
    ungroup()
  
  ceac <- winners %>%
    group_by(scenario, WTP) %>%
    summarise(prob_optimal = sum(weight) / n_distinct(nhb_df$iteration), .groups = "drop")
  
  ceac_full <- tidyr::complete(
    ceac,
    scenario = scenario_names,
    WTP = WTP_list,
    fill = list(prob_optimal = 0)
  )
  
  summary_long <- mean_nhb %>%
    left_join(ceac_full, by = c("scenario", "WTP")) %>%
    group_by(WTP) %>%
    arrange(desc(mean_NHB), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  
  list(
    nhb_df = nhb_df,
    mean_nhb = mean_nhb,
    ceac = ceac_full,
    summary_long = summary_long
  )
}


# CEACc plot #####
# Define fixed colours

cols5 <- scales::hue_pal()(5)
highlight_cols <- c(
  "A. SQ"                          = cols5[1],
  "B. SQ + PrEP (India regimen)"   = cols5[2],
  "I. SQ + intense MDV in owned"   = cols5[3],
  "K. WHO PEP + SQ MDV"            = cols5[4],
  "L. WHO PEP + MDV scale-up all"  = cols5[5],
  "Other scenarios"                = "grey72"
)



# plot function 
plot_ceac <- function(ceac_df, x_breaks = c(2034, 3092)) {
  
  ggplot(
    ceac_df,
    aes(
      x = WTP,
      y = prob_optimal,
      colour = scenario_group,
      group = scenario_label
    )
  ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1) +
    geom_vline(
      xintercept = x_breaks,
      linetype = "dotted",
      linewidth = 0.5,
      colour = "grey40"
    ) +
    annotate("text", x = 2100, y = 0.5, label = "WTP:\nUS$2,034", size = 3, vjust = 0, hjust=0) +
    annotate("text", x = 3100, y = 0.5, label = "WTP:\nUS$3,092", size = 3, vjust = 0, hjust=0) +
    scale_colour_manual(
      values = highlight_cols,
      breaks = names(highlight_cols),
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "WTP",
      y = "Probability optimal",
      colour = "Scenario"
    ) +
    theme_classic()
}


# Run 
## All scenarios 
res_all <- calc_nhb_ceac( Kerala_general = Kerala_general, 
                          WTP_list = WTP_list ) 
## Selected scenarios 
res_selected <- calc_nhb_ceac( Kerala_general = Kerala_general, 
                               WTP_list = WTP_list, 
                               scenarios = selected_scenarios ) 
# Add labels for plotting 
res_all$ceac <- add_scenario_labels(res_all$ceac) 
res_all$mean_nhb <- add_scenario_labels(res_all$mean_nhb) 
res_all$summary_long <- add_scenario_labels(res_all$summary_long) 

res_selected$ceac <- add_scenario_labels(res_selected$ceac) 
res_selected$mean_nhb <- add_scenario_labels(res_selected$mean_nhb) 
res_selected$summary_long <- add_scenario_labels(res_selected$summary_long)



# Plot
plot_all_df <- res_all$ceac %>%
  dplyr::mutate(
    scenario_group = ifelse(
      scenario_label %in% names(highlight_cols),
      scenario_label,
      "Other scenarios"
    )
  )

plot_ceac_all <- plot_ceac(plot_all_df) +
  theme(axis.title.x = element_blank())

plot_selected_df <- res_selected$ceac %>%
  dplyr::mutate(
    scenario_group = scenario_label
  ) 
plot_ceac_selected <- plot_ceac(plot_selected_df) + theme(
  legend.title = element_blank()
)



plot_ceac_all
plot_ceac_selected

plot_ceac_all/plot_ceac_selected  


pdf("./output/manuscript_figures/Figure3.pdf", width = 7.5, height = 4, useDingbats = FALSE)
print(plot_ceac_selected)
dev.off()



# Tables ######
## Table, all scenarios

# table formating 
format_nhb_table <- function(summary_long, wtps_to_keep = NULL) { 
  out <- summary_long 
  if (!is.null(wtps_to_keep)) { 
    out <- out %>% filter(WTP %in% wtps_to_keep) 
  } 
  out %>% 
    dplyr::mutate( mean_NHB_fmt = comma(round(mean_NHB, 0), 
                                        accuracy = 1), 
                   prob_optimal_pct = percent(prob_optimal, 
                                              accuracy = 1) 
    ) %>% 
    dplyr::select(scenario, scenario_label, WTP, rank, mean_NHB_fmt, prob_optimal_pct) 
}

# wide MS table helper 
make_wide_nhb_table <- function(summary_long, wtps_to_keep = NULL) { 
  out <- summary_long 
  if (!is.null(wtps_to_keep)) { 
    out <- out %>% dplyr::filter(WTP %in% wtps_to_keep) 
  } 
  out %>% 
    dplyr::mutate( scenario_label = recode(scenario, !!!scenario_labels, .default = scenario), 
                   nhb_prob_rank = paste0(comma(round(mean_NHB, 0), accuracy = 1),
                                          " (", percent(prob_optimal, accuracy = 1), #", rank ", rank, 
                                          ")" ) 
    ) %>% 
    dplyr::select(scenario_label, WTP, nhb_prob_rank) %>% 
    pivot_wider(names_from = WTP, values_from = nhb_prob_rank) 
}



nhb_table_all <- format_nhb_table(
  res_all$summary_long,
  wtps_to_keep = c(200, 500, 1000, 2034, 2535, 3092, 4000, 5000)
)
nhb_table_all


# To merge this with pretty_table in run_models.R (Table 3)
merge_table1<- nhb_table_all %>%
  dplyr::filter(WTP == 2535) %>%  # mean weighted WTP for India
  dplyr::select(scenario_label, mean_NHB_fmt, prob_optimal_pct, rank) %>%
  dplyr::rename(
    Scenario = scenario_label, 
    `Mean NHB at US$2,535` = mean_NHB_fmt,
    `Probability optimal` = prob_optimal_pct,
    Ranking = rank
  )

merge_table1


### wide
nhb_wide_all <- make_wide_nhb_table(
  res_all$summary_long,
  wtps_to_keep = c(200, 500, 1000, 2034, 2535, 3092, 4000, 5000)
) %>%
  arrange(scenario_label) %>%
  dplyr::rename_with(
    ~ ifelse(.x == "scenario_label", .x, paste0("US$", .x))
  )

nhb_wide_all

## Table, selected scenarios
nhb_table_selected <- format_nhb_table(
  res_selected$summary_long,
  wtps_to_keep = c(200, 500, 1000, 2034, 2535, 3092, 4000, 5000)
)
nhb_table_selected

### wide
nhb_wide_selected <- make_wide_nhb_table(
  res_selected$summary_long,
  wtps_to_keep = c(200, 500, 1000, 2034, 2535, 3092, 4000, 5000)
) %>%
  arrange(scenario_label) %>%
  dplyr::rename_with(
    ~ ifelse(.x == "scenario_label", .x, paste0("US$", .x))
  )

nhb_wide_selected



# Export Tables

write.csv(nhb_wide_all, "./output/Supp_Table2.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")

write.csv(nhb_wide_selected, "./output/Supp_Table3.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")






