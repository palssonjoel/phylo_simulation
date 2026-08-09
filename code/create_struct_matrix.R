

# This script works in creating a matrix which ONLY concerns Denmark.
# But I need to incorporate the trade data of the other countries between each other too.
# Create a function for this, and call it on the relevant dfs.

library(tidyverse)

trade_data <- read.csv("output/pig_stats_filtered.csv")
herd_data <- read.csv("data/eu_herd_data/herd_sizes_eu.csv")

# Convert trade data from metric tonnes to no. of heads
# For now, use some values found online:
#   Breeding sows = 250 kg/pig
#   Non-breeding pigs > 50kg = 120 kg
#   Non-breeding pigs < 50kg = 30 kg

weight_breeding <- 250
weight_large    <- 120
weight_small    <- 30

trade_data <- trade_data |>
  mutate(
    quantity_kg = quantity * 1000,
    quantity_heads = case_when(
      description ==
        "Swine: live, other than pure-bred breeding animals, weighing 50kg or more" ~
        quantity_kg / weight_large,
      
      description ==
        "Swine: live, other than pure-bred breeding animals, weighing less than 50kg" ~
        quantity_kg / weight_small,
      
      description ==
        "Swine: live, pure-bred breeding animals" ~
        quantity_kg / weight_breeding,
      
      TRUE ~ NA_real_
    )
  )

# Extract top 5 countries Denmark exports to. All three years added.
dk_top_exports <- trade_data |> 
  filter(exporter == "Denmark") |> 
  group_by(importer) |> 
  summarise(imp_sum = sum(quantity_heads)) |> 
  arrange(desc(imp_sum)) |> 
  slice(1:5)

countries <- c("Denmark", all_of(dk_top_exports$importer)) 

# Keep relevant countries and clean dataframe
trade_data_filtered <- trade_data |> 
  filter(exporter %in% countries,
         importer %in% countries) |> 
  select(-product, -value, -year, -description, -quantity_kg, -quantity, -X) |> 
  group_by(exporter, importer) |> 
  summarize(quantity_heads = sum(quantity_heads)) |> 
  ungroup()

# Prepare relevant herd data
herd_data_filtered <- herd_data |> 
  filter(geo %in% countries) |> 
  filter(TIME_PERIOD == 2024) |>                # Only 2024 for now
  filter(month == "November - December") |>     # End-of-year data only
  select(geo, OBS_VALUE, unit, animals) |> 
  group_by(geo) |> 
  summarize(herd_size = sum(OBS_VALUE) * 1000)  # Numbers are in 1000s of heads

# Merge trade and herd data
complete_data <- trade_data_filtered |> 
  left_join(herd_data_filtered,
            by = join_by(exporter == geo)) |> 
  rename(exporter_herd_size = herd_size,
         no_exported_animals = quantity_heads) |> 
  mutate(
    no_exported_animals = no_exported_animals / (365 * 3), # avg. daily trade
    
    # Per-pig daily movement rate
    rate_exp_to_imp = no_exported_animals / exporter_herd_size,
    
    # Probability of each individual moving to another country
    p_exp_to_imp = 1 - exp(-rate_exp_to_imp)
  )

# Create empty matrix
trade_matrix <- matrix(0,
                       nrow = length(countries),
                       ncol = length(countries),
                       dimnames = list(countries, countries))

# Populate the matrix
trade_matrix <- matrix(
  0,
  nrow = length(countries),
  ncol = length(countries),
  dimnames = list(countries, countries)
)

trade_matrix[
  cbind(
    complete_data$exporter,
    complete_data$importer
  )
] <- complete_data$p_exp_to_imp

# Save
write.csv(
  trade_matrix,
  "output/trade_matrix_daily_probabilities.csv",
  row.names = TRUE)

write.csv(
  complete_data,
  "output/trade_movement_rates.csv",
  row.names = FALSE)
