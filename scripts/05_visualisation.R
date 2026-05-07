




# Temporal visuals ######
# Prepare data

Kerala_general <- readRDS("output/Kerala_general.rds")

scenarios <- names(Kerala_general)

create_summary_across_horizon(
  variable_name = "ts_deaths",
  matrices = lapply(scenarios, function(s) Kerala_general[[s]]$ts_deaths),
  scenario_names = scenarios
)


# Contribution and cost of interventions ######

output_table # from run_model.R

output_table_stacked_df <- read.csv("./output/Kerala_general_out.csv") %>%
  dplyr::mutate(
    scenario = recode(
      scenario,
      SQ = "A. SQ", 
      SQ_addIndian_PrEP = "B. SQ + PrEP (India regimen)",
      WHO_PEP_1ml_IM_PrEP = "C. WHO PEP + 1-shot PrEP",
      WHO_PEP_and_PrEP = "D. WHO PEP & PrEP",
      WHO_PEP_and_PrEP_reducedMDVUnowned = "E. D + drop in MDV",
      WHO_PEPandPrEP_MDVscaleupUnowned = "F. D + MDV scale-up in unowned",
      SQ_mod_MDVUnowned = "G. SQ  +  moderate MDV in unowned",
      SQ_MDVscaleupUnowned = "H. SQ  +  MDV scale-up in unowned",
      SQ_MDVscaleupOwned = "I. SQ + intense MDV in owned",
      SQ_MDVscaleupAll = "J. SQ + MDV scale-up in all", 
      WHO_PEP_and_SQ_MDV = "K. WHO PEP + SQ MDV",
      WHO_PEP_MDVscaleupAll = "L. WHO PEP  +  MDV scale-up all"
      
    )
  )


# 1. Prepare deaths averted

deaths_long <- output_table_stacked_df %>%
  select(
    scenario,
    ts_deaths_averted_PEP_Median,
    ts_deaths_averted_PrEP_Median,
    ts_deaths_averted_MDV_Median
  ) %>%
  pivot_longer(
    -scenario,
    names_to = "intervention",
    values_to = "value"
  ) %>%
  mutate(
    intervention = recode(
      intervention,
      ts_deaths_averted_PEP_Median  = "PEP",
      ts_deaths_averted_PrEP_Median = "PrEP",
      ts_deaths_averted_MDV_Median  = "MDV"
    ),
    side = "Deaths averted"
  )


# 2. Prepare costs (RIGHT)

costs_long <- output_table_stacked_df %>%
  dplyr::select(
    scenario,
    ts_cost_PEP_per_year_Median,
    ts_cost_PrEP_Median,
    ts_MDV_campaign_cost_Median
  ) %>%
  pivot_longer(
    -scenario,
    names_to = "intervention",
    values_to = "value"
  ) %>%
  mutate(
    intervention = recode(
      intervention,
      ts_cost_PEP_per_year_Median = "PEP",
      ts_cost_PrEP_Median         = "PrEP",
      ts_MDV_campaign_cost_Median = "MDV"
    ),
    value = value / 1e6,   
    side = "Costs"
  )



# 3. Ensure consistent scenario ordering

scenario_order <- rev(output_table_stacked_df %>%
                        mutate(total = ts_deaths_averted_Median) %>%
                        #arrange(total) %>%
                        pull(scenario) )


# Order stacked bars and scenarios
stack_order <- c("MDV", "PEP", "PrEP")

deaths_long <- deaths_long %>%
  mutate(
    intervention = factor(intervention, levels = stack_order),
    scenario     = factor(scenario, levels = scenario_order)
  )

costs_long <- costs_long %>%
  mutate(
    intervention = factor(intervention, levels = stack_order),
    scenario     = factor(scenario, levels = scenario_order)
  )



# 4. Build the two mirrored stacked plots
## 4.1. Define SQ reference values
sq_deaths_averted <- deaths_long %>%
  dplyr::filter(scenario == "A. SQ") %>%
  dplyr::summarise(x = sum(value, na.rm = TRUE)) %>%
  dplyr::pull(x)

sq_costs <- costs_long %>%
  dplyr::filter(scenario == "A. SQ") %>%
  dplyr::summarise(x = sum(value, na.rm = TRUE)) %>%
  dplyr::pull(x)

sq_deaths <- output_table_stacked_df %>%
  dplyr::filter(scenario == "A. SQ") %>%
  dplyr::pull(ts_deaths_Median)

## 4.2. Make plots
## Part A: deaths averted #########
p_left <- ggplot(deaths_long, aes(x = value, y = scenario, fill = intervention)) +
  geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  geom_vline(
    xintercept = sq_deaths_averted,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    breaks = stack_order,
    values = c(
      MDV  = "#C44E52",
      PEP  = "#4C72B0",
      PrEP = "#55A868"
    )
  ) +
  labs(x = "Deaths averted", y = NULL) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 9)
  )

## Part B: costs #########
p_middle <- ggplot(costs_long, aes(x = value, y = scenario, fill = intervention)) +
  geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  geom_vline(
    xintercept = sq_costs,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.5
  ) +
  scale_fill_manual(
    breaks = stack_order,
    values = c(
      MDV  = "#C44E52",
      PEP  = "#4C72B0",
      PrEP = "#55A868"
    )
  ) +
  labs(x = "Costs (million US$)", y = NULL) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


## Part C: deaths #########
p_right <- output_table_stacked_df %>%
  select(scenario, ts_deaths_Median) %>%
  dplyr::mutate(scenario = factor(scenario, levels = scenario_order)) %>%
  ggplot(aes(x = ts_deaths_Median, y = scenario)) +
  geom_col(
    fill = "grey70",
    width = 0.7
  ) +
  geom_vline(
    xintercept = sq_deaths,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.5
  ) +
  scale_x_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(x = "Total deaths", y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )





# Combine

plot2A <- plot_grid(
  p_left,
  p_middle,
  p_right,
  nrow = 1,
  align = "h",
  axis = "tb",
  rel_widths = c(1.7, 1, 1)
)



# Figure 2 ########

## Now exported at the bottom (after supplementary)
Fig2 <- ggdraw(plot2A) 
Fig2

pdf("./output/manuscript_figures/Figure2.pdf", width = 10, height = 5, useDingbats = FALSE)
print(Fig2)
dev.off()


# Supp Fig 3: CE analysis#####

output_table2 <- read.csv("./output/Kerala_general_prepEff1_summarized.csv")
output_table <- read.csv("./output/Kerala_general_out.csv")

output_table # base model (prep_effectiveness = 0)
output_table2 # base model (prep_effectiveness = 1)


# --- helper: scenario labels (reuse in both runs)
scenario_labels <- c(
  SQ = "A. SQ",
  SQ_addIndian_PrEP = "B. SQ + PrEP (India regimen)",
  WHO_PEP_1ml_IM_PrEP = "C. WHO PEP + 1-shot PrEP",
  WHO_PEP_and_PrEP = "D. WHO PEP & PrEP",
  WHO_PEP_and_PrEP_reducedMDVUnowned = "E. D + drop in MDV",
  WHO_PEPandPrEP_MDVscaleupUnowned = "F. D + MDV scale-up in unowned",
  SQ_mod_MDVUnowned = "G. SQ  +  moderate MDV in unowned",
  SQ_MDVscaleupUnowned = "H. SQ  +  MDV scale-up in unowned",
  SQ_MDVscaleupOwned = "I. SQ + intense MDV in owned",
  SQ_MDVscaleupAll = "J. SQ + MDV scale-up in all",
  WHO_PEP_and_SQ_MDV = "K. WHO PEP + SQ MDV",
  WHO_PEP_MDVscaleupAll = "L. WHO PEP  +  MDV scale-up all"
)

# --- ensure output_table2 (PrEP_effectiveness = 1) has scenario as character, and drop the "X" column if present
output_table2_clean <- output_table2 %>%
  as_tibble() %>%
  select(-any_of("X")) %>%
  dplyr::mutate(scenario = as.character(scenario))

output_table_clean <- output_table %>%
  dplyr::mutate(scenario = as.character(scenario))

# --- SQ baseline from base run (PrEP_effectiveness = 0)
SQ0 <- output_table_clean %>%
  dplyr::filter(scenario == "SQ") %>%
  summarise(
    base_effect = ts_deaths_Median,
    base_cost   = ts_cost_per_year_Median,
    .groups = "drop"
  )

# --- base run: x + y and labels
base_df <- output_table_clean %>%
  dplyr::mutate(
    inc_deaths_averted_0 = SQ0$base_effect - ts_deaths_Median,
    inc_cost_musd        = (ts_cost_per_year_Median - SQ0$base_cost) / 1e6,
    PrEP_present         = str_detect(scenario, "PrEP"),
    PrEP_label           = ifelse(PrEP_present, "PrEP", "No PrEP"),
    scenario_label       = recode(scenario, !!!scenario_labels)
  ) %>%
  select(scenario, scenario_label, PrEP_present, PrEP_label,
         inc_deaths_averted_0, inc_cost_musd)

# --- PrEP effectiveness = 1 run: only need the incremental deaths averted vs SQ0 baseline
eff1_df <- output_table2_clean %>%
  dplyr::mutate(
    inc_deaths_averted_1 = SQ0$base_effect - ts_deaths_Median
  ) %>%
  select(scenario, inc_deaths_averted_1)

# --- join and define xmin/xmax (only meaningful for PrEP scenarios)
selected_scenarios <- c("SQ", "WHO_PEP_and_SQ_MDV", "SQ_addIndian_PrEP", "SQ_MDVscaleupOwned", "WHO_PEP_MDVscaleupAll")

ce_df2 <- base_df %>%
  left_join(eff1_df, by = "scenario") %>%
  dplyr::mutate(
    #  compute ranges
    xmin = ifelse(
      PrEP_present,
      pmin(inc_deaths_averted_0, inc_deaths_averted_1, na.rm = TRUE),
      NA_real_
    ),
    xmax = ifelse(
      PrEP_present,
      pmax(inc_deaths_averted_0, inc_deaths_averted_1, na.rm = TRUE),
      NA_real_
    ),
    
    x = inc_deaths_averted_0,
    
    scenario_highlight = ifelse(
      scenario %in% selected_scenarios,
      scenario_label,
      "Other scenarios"
    )
  )



# Plot

highlight_cols <- c(
  "A. SQ"                            = "#F8766D",
  "B. SQ + PrEP (India regimen)"     = "#A3A500",
  "I. SQ + intense MDV in owned"     = "#00BF7D",
  "K. WHO PEP + SQ MDV"              = "#00B0F6",
  "L. WHO PEP  +  MDV scale-up all"  = "#E76BF3",
  "Other scenarios"                  = "grey50"
)


# CE_frontier helper function in Helper.R

# Build ICER frontier data
source("./scripts/01_HelperFun.R")
icer_frontier_df <- get_icer_frontier(ce_df2)

# Plot

ce_plane <- ggplot(ce_df2,aes(
  x = x, y = inc_cost_musd, colour = scenario_highlight, shape  = PrEP_label
)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey70") +
  
  geom_errorbar(
    data = ~ dplyr::filter(.x, PrEP_present),
    aes(xmin = xmin, xmax = xmax),
    orientation = "y",
    width = 0,
    linewidth = 0.6,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  
  geom_point(size = 3.5, alpha = 0.9) +
  
  geom_path(
    data = icer_frontier_df,
    aes(x = x, y = inc_cost_musd),
    inherit.aes = FALSE,
    linewidth = 0.7,
    colour = "black"
  ) +
  
  geom_text_repel(
    aes(label = scenario_label),
    size = 3,
    max.overlaps = 12,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = highlight_cols,
    breaks = c(
      "A. SQ",
      "B. SQ + PrEP (India regimen)",
      "I. SQ + intense MDV in owned",
      "K. WHO PEP + SQ MDV",
      "L. WHO PEP  +  MDV scale-up all",
      "Other scenarios"
    ),
    guide = "none"
  ) +
  
  scale_shape_manual(
    values = c(
      "No PrEP"      = 16,
      "PrEP" = 17
    ),
    guide = guide_legend(
      override.aes = list(
        shape = c(1, 2),   # hollow circle, hollow triangle
        colour = "black",
        fill = NA,
        size = 3,
        alpha = 1
      )
    )
  ) +
  
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  
  labs(
    x = "Incremental deaths averted vs SQ (median)",
    y = "Incremental cost vs SQ (million US$)",
    colour = "Scenario",
    shape  = "Intervention"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10)
  )

ce_plane


## Run Net Health Benefit Analysis
source("./scripts/07_nhb_analysis.R")

plot_ceac_all

ce_plane / plot_ceac_all



combined_plot <- ce_plane / plot_ceac_all + plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    #legend.title = element_text(size = 9),
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )

combined_plot

pdf("./output/manuscript_figures/SupplementaryFigure3.pdf", width = 14, height = 5.5, useDingbats = FALSE)
print(combined_plot)
dev.off()

# Supplementary Figure 2#####

plot_temporal_df<- summarise_multiple_vars(
  variables = c("ts_exposures", "ts_deaths", "ts_cost_per_year"),
  scenario_list = Kerala_general,
  scenario_names = scenarios
) %>%
  dplyr::mutate(
    scenario = recode(
      scenario,
      SQ = "A. SQ", 
      SQ_addIndian_PrEP = "B. SQ + PrEP (India regimen)",
      WHO_PEP_1ml_IM_PrEP = "C. WHO PEP + 1-shot PrEP",
      WHO_PEP_and_PrEP = "D. WHO PEP & PrEP",
      WHO_PEP_and_PrEP_reducedMDVUnowned = "E. D + drop in MDV",
      WHO_PEPandPrEP_MDVscaleupUnowned = "F. D + MDV scale-up in unowned",
      SQ_mod_MDVUnowned = "G. SQ  +  moderate MDV in unowned",
      SQ_MDVscaleupUnowned = "H. SQ  +  MDV scale-up in unowned",
      SQ_MDVscaleupOwned = "I. SQ + intense MDV in owned",
      SQ_MDVscaleupAll = "J. SQ + MDV scale-up in all", 
      WHO_PEP_and_SQ_MDV = "K. WHO PEP + SQ MDV",
      WHO_PEP_MDVscaleupAll = "L. WHO PEP  +  MDV scale-up all"
      
    )
  )




plot_df <- plot_temporal_df %>%
  filter(variable %in% c("ts_exposures", "ts_deaths", "ts_cost_per_year")) %>%
  mutate(
    Median = if_else(variable == "ts_cost_per_year", Median / 1e6, Median),
    LL     = if_else(variable == "ts_cost_per_year", LL / 1e6, LL),
    UL     = if_else(variable == "ts_cost_per_year", UL / 1e6, UL),
    
    metric = recode(
      variable,
      "ts_exposures" = "Exposures",
      "ts_deaths" = "Deaths",
      "ts_cost_per_year" = "Costs"
    ),
    metric = factor(metric, levels = c("Exposures", "Deaths", "Costs")),
    scenario = factor(scenario, levels = unique(scenario))
  )

base_theme <- theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    axis.text.x = element_text(angle = 0)
  )

make_temporal_col <- function(df, metric_name, y_lab) {
  ggplot(
    df %>% filter(metric == metric_name),
    aes(x = year, y = Median, colour = scenario, fill = scenario)
  ) +
    geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 0.8) +
    facet_grid(rows = vars(scenario), scales = "fixed") +
    scale_y_continuous(labels = scales::comma) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    labs(
      title = metric_name,
      x = "Year",
      y = y_lab
    ) +
    base_theme
}

p_exp <- make_temporal_col(plot_df, "Exposures", "Exposures") +
  theme(
    strip.text.y = element_blank(),
    strip.background = element_blank()
  )

p_deaths <- make_temporal_col(plot_df, "Deaths", "Deaths") +
  theme(
    strip.text.y = element_blank(),
    axis.title.y = element_text()
  )

p_costs <- make_temporal_col(plot_df, "Costs", "Costs (US$ million)") +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.background = element_blank()
  )

myplot <- p_exp + p_deaths + p_costs +
  plot_layout(ncol = 3, widths = c(1, 1, 1)) +
  patchwork::plot_annotation(tag_levels = "A")

myplot

# Export to 2 (pdf) pages

scenarios <- levels(plot_df$scenario)
chunks <- split(scenarios, ceiling(seq_along(scenarios) / 6))

make_page <- function(df, scen_subset) {
  df_sub <- df %>% filter(scenario %in% scen_subset)
  
  p_exp <- make_temporal_col(df_sub, "Exposures", "Exposures") +
    theme(strip.text.y = element_blank(), strip.background = element_blank())
  
  p_deaths <- make_temporal_col(df_sub, "Deaths", "Deaths") +
    theme(strip.text.y = element_blank())
  
  p_costs <- make_temporal_col(df_sub, "Costs", "Costs (US$ million)") +
    theme(strip.text.y = element_text(angle = 0, hjust = 0),
          strip.background = element_blank())
  
  p_exp + p_deaths + p_costs +
    plot_layout(ncol = 3) +
    plot_annotation(tag_levels = "A")
}

pdf("./output/manuscript_figures/SupplementaryFigure2.pdf", width = 10, height = 12) 
for (i in seq_along(chunks)) {
  print(make_page(plot_df, chunks[[i]]))
}
dev.off()

# pdf("./output/manuscript_figures/SupplementaryFigure2.pdf", width = 8, height = 8, useDingbats = FALSE)
# print(myplot)
# dev.off()




# 
# ## Extra ####
# # Exposures temporally
# 
# plot_exp_temporal <-plot_temporal_df %>%
#   dplyr::filter(variable == "ts_exposures") %>%
#   ggplot(., aes(x = year, y = Median)) +
#   geom_line() +
#   geom_ribbon(
#     aes(ymin = LL, ymax = UL),
#     fill = "purple",
#     alpha = 0.4,
#     color = NA
#   ) +
#   labs(x = "Year", y = "Exposures") +
#   theme_bw() +
#   scale_y_continuous(labels = scales::comma) +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
#   #scale_color_manual(values = line_palette) +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
#   facet_wrap(~scenario) +
#   theme_classic() +
#   theme(strip.background = element_blank())
# 
# plot_exp_temporal
# 
# 
# # Deaths temporally
# 
# plot_deaths_temporal<- plot_temporal_df %>%
#   dplyr::filter(variable == "ts_deaths") %>%
#   ggplot(., aes(x = year, y = Median)) +
#   geom_line() +
#   geom_ribbon(
#     aes(ymin = LL, ymax = UL),
#     fill = "purple",
#     alpha = 0.4,
#     color = NA
#   ) +
#   labs(x = "Year", y = "Deaths") +
#   theme_bw() +
#   scale_y_continuous(labels = scales::comma) +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
#   #scale_color_manual(values = line_palette) +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
#   facet_wrap(~scenario) +
#   theme_classic() +
#   theme(strip.background = element_blank())
# 
# plot_deaths_temporal
# 
# 
# # Costs temporally
# 
# plot_costs_temporal<- plot_temporal_df %>%
#   dplyr::filter(variable == "ts_cost_per_year") %>%
#   ggplot(., aes(x = year, y = Median)) +
#   geom_line() +
#   geom_ribbon(
#     aes(ymin = LL, ymax = UL),
#     fill = "purple",
#     alpha = 0.4,
#     color = NA
#   ) +
#   labs(x = "Year", y = "Annual costs") +
#   theme_bw() +
#   scale_y_continuous(labels = scales::comma) +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 5))+
#   #scale_color_manual(values = line_palette) +
#   theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
#   facet_wrap(~scenario) +
#   theme_classic() +
#   theme(strip.background = element_blank())
# 
# 
# plot_costs_temporal
# 
# plot_exp_temporal/ plot_deaths_temporal +
#   plot_annotation(tag_levels = "A")
# 



