
## Run analyses ##### 
# explore a range of parameters

# Some values/ R objects from run_model.R are used here eg output_table

# 1. PrEP effectiveness
# 2. pSeek
# 3. Number of unwoned dogs
# 4. MDV -- in unowned dogs
## Everything varied alongside PrEP coverage (first year, to be topped up by birth cohort)

# Dog numbers ???? Unowned

# we perform a 2 way sensitivity analysis ie everything relative to PrEP coverage 

# 1. Establish base case (needs to have PrEP active)
# scenarios from run_model.R
SQ_base <- scenarios1 %>%
  dplyr::filter(scenario == "WHO_PEP_and_PrEP") %>%
  dplyr::select(-scenario) %>%
  as.list()


# 2. Define sensitivity ranges

unowned_lookup <- tibble::tribble(
  ~unowned_mult, ~HDR,         ~unowned_prop,
  1,             31.11111111,  0.257,
  2,             24.75247525,  0.409,
  3,             20.55196712,  0.509,
  4,             17.57028112,  0.580,
  5,             15.34414730,  0.633#,
  # 6,             13.61867704,  0.675,
  # 7,             12.24204267,  0.708
)


ranges <- list(
  init_PrEP_cov      = seq(0, 0.9, by = 0.1),
  PrEP_effectiveness = seq(0, 1, by = 0.1),
  pSeek_exposure     = seq(0.5, 0.95, by = 0.05),
  mdv_cov            = seq(0, 0.7, by = 0.1),
  unowned_mult       = unowned_lookup$unowned_mult
)

# ranges <- list(
#   init_PrEP_cov      = seq(0, 0.9, by = 0.1),
#   PrEP_effectiveness = seq(0, 1, by = 0.1),
#   pSeek_exposure     = seq(0.5, 0.95, by = 0.05),
#   mdv_cov            = seq(0, 0.7, by = 0.1)
# )


# 3. Build two-way sensitivity design table
sensitivity_grid <- bind_rows(
  
  # PrEP effectiveness × PrEP coverage
  tidyr::expand_grid(
    analysis = "PrEP_eff_x_PrEP_cov",
    init_PrEP_cov = ranges$init_PrEP_cov,
    PrEP_effectiveness = ranges$PrEP_effectiveness
  ),
  
  # pSeek × PrEP coverage
  tidyr::expand_grid(
    analysis = "pSeek_x_PrEP_cov",
    init_PrEP_cov = ranges$init_PrEP_cov,
    pSeek_exposure = ranges$pSeek_exposure
  ),
  
  # MDV × PrEP coverage
  tidyr::expand_grid(
    analysis = "MDV_x_PrEP_cov",
    init_PrEP_cov = ranges$init_PrEP_cov,
    mdv_cov = ranges$mdv_cov
  ),
  
  # unowned dog multiplier × PrEP coverage
  tidyr::expand_grid(
    analysis = "unowned_mult_x_PrEP_cov",
    init_PrEP_cov = ranges$init_PrEP_cov,
    unowned_mult = ranges$unowned_mult
  )
  
)


# 4. Translate grid → model arguments (with enforced rules)

build_2wsa_args <- function(
    analysis,
    init_PrEP_cov,
    PrEP_effectiveness = NA,
    pSeek_exposure     = NA,
    mdv_cov            = NA,
    unowned_mult       = NA, # label only — not written to args
    SQ_base,
    common_args,
    unowned_lookup
) {
  
  args <- do.call(
    purrr::list_modify,
    c(list(common_args), SQ_base)
  )
  
  args$N <- 1000
  args$init_PrEP_cov <- init_PrEP_cov
  
  if (!is.na(PrEP_effectiveness)) {
    args$PrEP_effectiveness <- PrEP_effectiveness
  }
  
  if (!is.na(pSeek_exposure)) {
    args$pSeek_exposure <- pSeek_exposure
  }
  
  if (!is.na(mdv_cov)) {
    args$target_vax_cov_unowned <- mdv_cov
  }
  
  if (!is.na(unowned_mult)) {
    this_row <- unowned_lookup %>%
      dplyr::filter(unowned_mult == !!unowned_mult)
    
    if (nrow(this_row) != 1) {
      stop("unowned_mult did not match exactly one row in unowned_lookup")
    }
    
    args$HDR <- c(this_row$HDR, this_row$HDR)
    args$unowned_prop <- this_row$unowned_prop
  }
  
  args
}





# 5. Run the Two-Way SA in parallel


library(furrr)
plan(multisession, workers = 8)


# Seeds and burn-in
burnin <- common_args1$dog_burnin
seeds  <- 300L + seq_len(nrow(sensitivity_grid))

## argument construction
args_list <- sensitivity_grid %>%
  pmap(
    build_2wsa_args,
    SQ_base = SQ_base,
    common_args = common_args1,
    unowned_lookup = unowned_lookup
  )


output_list <- future_map2(
  .x = args_list,
  .y = seeds,
  .f = ~ {
    load_rabies_models()                          # EACH WORKER/ CORE NEEDS A COPY OF THE MODELS IN MEMORY
    out <- do.call(decision_tree_PrEP, purrr::list_modify(.x, seed = .y))
    out   # burn-in already stripped inside decision_tree_PrEP
  },
  .options = furrr::furrr_options(seed = TRUE)
)



TWSA_results <- sensitivity_grid %>%
  dplyr::mutate(
    seed = seeds,
    args = args_list,
    output = output_list
  )


## Summarise and plot #########

# save 
# saveRDS(TWSA_results, file = "output/two_way_sensitivity.rds")

TWSA_results <- readRDS("output/two_way_sensitivity.rds")

# Visualize 2 way sensitivity results
## Calculate ICERs
# we need the SQ costs and effects from output table (run_model.R) to be comparator.
## Used in compute_TWSA_ICERs function
names(TWSA_results)
output_table

# 1. Helpers
summarise_across_horizon<- function(
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


# 2. Summarise TWSA results
summarise_TWSA_results <- function(
    TWSA_results,
    vars_to_summarise = c(
      "ts_deaths",
      "ts_deaths_averted",
      "ts_cost_per_year"
    ),
    probs = c(0.025, 0.5, 0.95)
) {
  
  out <- TWSA_results %>%
    dplyr::select(-args) %>%   # args not needed downstream
    dplyr::mutate(
      
      summaries = purrr::map(
        output,
        function(res) {
          
          purrr::imap(
            vars_to_summarise,
            function(v, nm) {
              
              if (!v %in% names(res)) {
                stop("Variable ", v, " not found in TWSA output.")
              }
              
              summarise_across_horizon(
                res[[v]],
                probs = probs
              ) %>%
                dplyr::rename_with(
                  ~ paste0(v, "_", .x)
                )
            }
          ) %>%
            dplyr::bind_cols()
        }
      )
    ) %>%
    tidyr::unnest(summaries) %>%
    dplyr::select(-output)
  
  out
}


twsa_summary <- summarise_TWSA_results(
  TWSA_results,
  vars_to_summarise = c(
    "ts_deaths",
    "ts_deaths_averted",
    "ts_cost_per_year"
  )
)

# 3. ICERs relative to SQ (main analysis)

compute_TWSA_ICERs <- function(
    twsa_summary,
    output_table,
    baseline_scenario = "SQ"
) {
  
  ## Extract SQ baseline
  base_cost <- output_table %>%
    dplyr::filter(scenario == baseline_scenario) %>%
    dplyr::pull(ts_cost_per_year_Median)
  
  base_effect <- output_table %>%
    dplyr::filter(scenario == baseline_scenario) %>%
    dplyr::pull(ts_deaths_averted_Median)
  
  if (length(base_cost) != 1 || length(base_effect) != 1) {
    stop("Baseline scenario must resolve to exactly one row.")
  }
  
  twsa_summary %>%
    dplyr::mutate(
      inc_cost   = ts_cost_per_year_Median - base_cost,
      inc_effect = ts_deaths_averted_Median - base_effect,
      
      ICER = dplyr::if_else(
        inc_effect > 0,
        inc_cost / inc_effect,
        NA_real_
      )
    )
}


output_table <- read.csv("./output/Kerala_general_out.csv")
twsa_icers <- compute_TWSA_ICERs(
  twsa_summary,
  output_table
)


prepare_deaths_surface <- function(
    twsa_summary,
    x_var,
    y_var,
    deaths_col = "ts_deaths_Median"  # 95th percentile
) {
  
  twsa_summary %>%
    dplyr::filter(
      !is.na(.data[[x_var]]),
      !is.na(.data[[y_var]])
    ) %>%
    dplyr::select(
      x = !!rlang::sym(x_var),
      y = !!rlang::sym(y_var),
      deaths = !!rlang::sym(deaths_col)
    )
}

plot_deaths_surface <- function( 
    surface_df, xlab, ylab ) { 
  ggplot(surface_df, 
         aes(x = x, y = y, z = deaths)) + 
    geom_contour_filled() + 
    labs( x = xlab, y = ylab, fill = "Deaths (95th PI)" ) + 
    theme_classic() }


make_surface_for_analysis <- function(
    twsa_summary,
    analysis_name,
    x_var,
    y_var,
    deaths_col = "ts_deaths_Median"
) {
  
  twsa_summary %>%
    dplyr::filter(analysis == analysis_name) %>%
    prepare_deaths_surface(
      x_var = x_var,
      y_var = y_var,
      deaths_col = deaths_col
    )
}


surface_PrEP_eff <- make_surface_for_analysis(
  twsa_summary,
  analysis_name = "PrEP_eff_x_PrEP_cov",
  x_var = "init_PrEP_cov",
  y_var = "PrEP_effectiveness",
  deaths_col = "ts_deaths_Median"
)

surface_pSeek <- make_surface_for_analysis(
  twsa_summary,
  analysis_name = "pSeek_x_PrEP_cov",
  x_var = "init_PrEP_cov",
  y_var = "pSeek_exposure",
  deaths_col = "ts_deaths_Median"
)

surface_MDV <- make_surface_for_analysis(
  twsa_summary,
  analysis_name = "MDV_x_PrEP_cov",
  x_var = "init_PrEP_cov",
  y_var = "mdv_cov",
  deaths_col = "ts_deaths_Median"
)


surface_unowned_dog <- make_surface_for_analysis(
  twsa_summary,
  analysis_name = "unowned_mult_x_PrEP_cov",
  x_var = "init_PrEP_cov",
  y_var = "unowned_mult",
  deaths_col = "ts_deaths_Median"
)


p1 <- plot_deaths_surface(
  surface_PrEP_eff,
  xlab = "PrEP coverage",
  ylab = "PrEP effectiveness"
)

p2 <- plot_deaths_surface(
  surface_pSeek,
  xlab = "PrEP coverage",
  ylab = "Healthcare seeking probability"
)

p3 <- plot_deaths_surface(
  surface_MDV,
  xlab = "PrEP coverage",
  ylab = "MDV coverage in unowned dogs"
)


p4 <- plot_deaths_surface(
  surface_unowned_dog,
  xlab = "PrEP coverage",
  ylab = "Unowned dog multiplier"
)

( p1 + p2 ) / (p3 + p4)

myplot<- ( p1 + p2 ) / (p3 + p4) +
  plot_annotation(tag_levels = "A")

myplot


pdf("./output/manuscript_figures/Figure4.pdf", width = 9, height = 7, useDingbats = FALSE)
print(myplot)
dev.off()



