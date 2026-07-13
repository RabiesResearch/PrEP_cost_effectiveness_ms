

run_discounted_prep <- function(prep_price,
                                scenarios1,
                                common_args1,
                                N = 3000,
                                seed = 9001) {
  
  scenario_row <- scenarios1 %>%
    dplyr::filter(scenario == "SQ_addIndian_PrEP") %>%
    dplyr::mutate(
      human_vaccine_cost_per_vial_PrEP = prep_price
    )
  
  common_args_tmp <- common_args1
  common_args_tmp$N <- N
  
  run_scenario(
    scenario_row = scenario_row,
    common_args = common_args_tmp,
    seed = seed,
    discount_rate = common_args1$discount
  )
}

Kerala_general_prep_price <- Kerala_general

# Kerala_general_prep_price[["SQ_addIndian_PrEP_0.5price"]] <-
#   run_discounted_prep(
#     prep_price = 6.25 * 0.50,
#     scenarios1 = scenarios1,
#     common_args1 = common_args1,
#     N = 3000,
#     seed = 9001
#   )
# 
# Kerala_general_prep_price[["SQ_addIndian_PrEP_0.25price"]] <-
#   run_discounted_prep(
#     prep_price = 6.25 * 0.25,
#     scenarios1 = scenarios1,
#     common_args1 = common_args1,
#     N = 3000,
#     seed = 9002
#   )
# 
# saveRDS(Kerala_general_prep_price, "./output/Kerala_general_prep_price.rds")
Kerala_general_prep_price <- readRDS("./output/Kerala_general_prep_price.rds")

names(Kerala_general_prep_price)

selected_scenarios_prep_price <- c(
  "SQ",
  "WHO_PEP_and_SQ_MDV",
  "SQ_addIndian_PrEP",
  "SQ_addIndian_PrEP_0.5price",
  "SQ_addIndian_PrEP_0.25price",
  "SQ_MDVscaleupOwned",
  "WHO_PEP_MDVscaleupAll"
)

scenario_labels_prep_price <- c(
  scenario_labels,
  SQ_addIndian_PrEP_0.5price  = "B2. SQ + PrEP, 50% vial price",
  SQ_addIndian_PrEP_0.25price = "B3. SQ + PrEP, 25% vial price"
)

add_scenario_labels_prep_price <- function(df) {
  df %>%
    dplyr::mutate(
      scenario_label = dplyr::recode(
        scenario,
        !!!scenario_labels_prep_price,
        .default = scenario
      )
    )
}


res_prep_price <- calc_nhb_ceac(
  Kerala_general = Kerala_general_prep_price,
  WTP_list = WTP_list,
  scenarios = selected_scenarios_prep_price
)

res_prep_price$ceac <- add_scenario_labels_prep_price(res_prep_price$ceac)
res_prep_price$summary_long <- add_scenario_labels_prep_price(res_prep_price$summary_long)
res_prep_price$mean_nhb <- add_scenario_labels_prep_price(res_prep_price$mean_nhb)

#. plot

cols7 <- scales::hue_pal()(7)

highlight_cols_prep_price <- c(
  "A. SQ"                            = cols7[1],
  "B. SQ + PrEP (India regimen)"     = cols7[2],
  "B2. SQ + PrEP, 50% vial price"    = cols7[3],
  "B3. SQ + PrEP, 25% vial price"    = cols7[4],
  "I. SQ + intense MDV in owned"     = cols7[5],
  "K. WHO PEP + SQ MDV"              = cols7[6],
  "L. WHO PEP + MDV scale-up all"    = cols7[7]
)

plot_prep_price_df <- res_prep_price$ceac %>%
  dplyr::mutate(scenario_group = scenario_label)

price_plot <- ggplot(
  plot_prep_price_df,
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
    xintercept = c(2034, 3092),
    linetype = "dotted",
    linewidth = 0.5,
    colour = "grey40"
  ) +
  scale_colour_manual(values = highlight_cols_prep_price, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  annotate("text", x = 2100, y = 0.5, label = "WTP:\nUS$2,034", size = 3, vjust = 0, hjust=0) +
  annotate("text", x = 3100, y = 0.5, label = "WTP:\nUS$3,092", size = 3, vjust = 0, hjust=0) +
  
  labs(
    x = "WTP",
    y = "Probability optimal",
    colour = "Scenario"
  ) +
  theme_classic()


pdf("./output/manuscript_figures/SupplementaryFigure4.pdf", width = 7.5, height = 4, useDingbats = FALSE)
print(price_plot)
dev.off()
