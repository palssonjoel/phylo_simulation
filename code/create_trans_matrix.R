# ==============================================================================
# Pig movement matrix construction
# ==============================================================================
#
# PURPOSE
# -------
# Construct a country-to-country pig movement structure matrix for use in
# a swIAV simulation in Nosoi, focusing on Denmark and its five largest
# export destinations.
#
# APPROACH
# --------
# 1. Convert reported live-swine trade from metric tonnes to estimated numbers
#    of pigs using assumed mean live weights for each trade category.
# 2. Aggregate trade over the available three-year period.
# 3. Convert three-year trade totals to average daily trade by dividing by
#    3 * 365.
# 4. Divide daily trade from exporter i to importer j by the total pig
#    population in exporter i to obtain a per-pig daily movement rate.
# 5. Within each exporter, normalize movement rates so that they sum to one.
#    These relative probabilities represent the destination of a pig
#    conditional on an international movement occurring.
# 6. Construct the resulting destination-probability matrix for use as the
#    Nosoi structure.matrix.
#
# IMPORTANT DISTINCTION
# ---------------------
# The structure matrix does NOT represent the absolute probability that a
# pig moves between countries. Rows represent the probability of destination
# conditional on a movement occurring and therefore sum to one.
#
# The absolute movement rate is retained separately in `complete_data` as
# `rate_exp_to_imp`.
#
# KEY ASSUMPTIONS
# ---------------
# - Denmark is used as the reference country because the study concerns
#   asymmetric trade and potential sampling bias associated with Denmark's
#   substantially greater exports than imports.
# - The simulation includes Denmark and its five largest export destinations.
# - Trade outside this selected network is excluded from the simulation.
# - Trade reported in metric tonnes is converted to numbers of pigs using
#   assumed mean weights:
#       breeding pigs:       250 kg
#       non-breeding >=50kg: 120 kg
#       non-breeding <50kg:   30 kg
# - The three-year trade total is assumed to represent a constant average
#   daily movement rate.
# - The November–December 2024 "Live swine, domestic species" population is
#   used as the denominator for each country's movement rate.
# - "Live swine, domestic species" is assumed to represent the total domestic
#   pig population, avoiding double-counting of its component categories.
#
# OUTPUTS
# -------
# trade_matrix_destination_probabilities.csv
#   Destination probabilities for the Nosoi structure.matrix.
#
# trade_movement_rates.csv
#   Underlying daily trade and per-pig movement rates used to construct
#   the matrix.
#
# ============================================================================== 

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

trade_data |>
  filter(is.na(quantity_heads)) |>
  count(description, sort = TRUE)

# Extract top 5 countries Denmark exports to. All three years added.
#dk_top_exports <- trade_data |> 
#  filter(exporter == "Denmark") |> 
#  group_by(importer) |> 
#  summarise(imp_sum = sum(quantity_heads)) |> 
#  arrange(desc(imp_sum)) |> 
#  slice(1:5)

#countries <- c(
#  "Denmark",
#  dk_top_exports$importer
#)

# Keep only relevant countries and clean dataframe
#trade_data_filtered <- trade_data |> 
#  filter(exporter %in% countries,
#         importer %in% countries) |> 
#  select(-product, -value, -year, -description, -quantity_kg, -quantity, -X) |> 
#  group_by(exporter, importer) |> 
#  summarise(
#    quantity_heads = sum(quantity_heads, na.rm = TRUE),
#    .groups = "drop"
#  )

# Instead of using top 5 importers of Denmark, use only imports/exports in EU
# anything else is considered "other"

# EU countries
countries <- c(
  "Albania",
  "Andorra",
  "Armenia",
  "Austria",
  "Azerbaijan",
  "Belarus",
  "Belgium",
  "Bosnia and Herzegovina",
  "Bulgaria",
  "Croatia",
  "Cyprus",
  "Czechia",
  "Denmark",
  "Estonia",
  "Finland",
  "France",
  "Georgia",
  "Germany",
  "Greece",
  "Hungary",
  "Iceland",
  "Ireland",
  "Italy",
  "Kazakhstan",
  "Kosovo",
  "Latvia",
  "Liechtenstein",
  "Lithuania",
  "Luxembourg",
  "Malta",
  "Moldova",
  "Monaco",
  "Montenegro",
  "Netherlands",
  "North Macedonia",
  "Norway",
  "Poland",
  "Portugal",
  "Romania",
  "Russia",
  "San Marino",
  "Serbia",
  "Slovakia",
  "Slovenia",
  "Spain",
  "Sweden",
  "Switzerland",
  "Turkey",
  "Ukraine",
  "United Kingdom",
  "Vatican City",
  "Other"
)

# Filter EU countries. Anything else is "other"
trade_data_filtered <- trade_data |> 
  mutate(
    importer = case_when(
      importer %in% countries ~ importer,
      TRUE ~ "Other"
    ),
    exporter = case_when(
      exporter %in% countries ~ exporter,
      TRUE ~ "Other"
    )
  ) |> 
  select(-product, -value, -year, -description, -quantity_kg, -quantity, -X) |> 
  group_by(exporter, importer) |> 
  summarise(
    quantity_heads = sum(quantity_heads, na.rm = TRUE),
    .groups = "drop"
  )

# Prepare relevant herd data. I'm doing some assumptions in the filtering here:
#   - Only 2024 data is used, and only the end-of-year values. Assumption is that
#     individual animals exist in both time periods, so using both would inflate
#     the numbers.
#   - Only "Live swine, domestic species" is used, assuming that this contains
#     a complete sum of all animals in the population. 

herd_data_filtered <- herd_data |> 
  filter(
    geo %in% countries,
    TIME_PERIOD == 2024,
    month == "November - December",
    animals == "Live swine, domestic species") |>   
  select(geo, OBS_VALUE, unit, animals) |> 
  group_by(geo) |> 
  summarize(herd_size = sum(OBS_VALUE) * 1000) |> # Numbers are in 1000s of heads
  ungroup()

# Merge trade and herd data
complete_data <- trade_data_filtered |> 
  left_join(herd_data_filtered,
            by = join_by(exporter == geo)) |> 
  rename(exporter_herd_size = herd_size,
         no_exported_animals = quantity_heads) |> 
  mutate(
    no_exported_animals = no_exported_animals / (365 * 3), # avg. daily trade
    
    # Per-pig daily movement rate
    rate_exp_to_imp = no_exported_animals / exporter_herd_size
  ) |> 
  group_by(exporter) |> 
  mutate(
    
    # Relative movement probability, sums to 1
    rel_p_exp_to_imp =
      rate_exp_to_imp / sum(rate_exp_to_imp)
  ) |> 
  ungroup() |> 
  drop_na(exporter_herd_size)
  

# Create empty matrix
trade_matrix <- matrix(0,
                       nrow = length(countries),
                       ncol = length(countries),
                       dimnames = list(countries, countries))

# Populate the matrix
trade_matrix[
  cbind(
    complete_data$exporter,
    complete_data$importer
  )
] <- complete_data$rel_p_exp_to_imp

# Keep countries with at least some outgoing trade
active_countries <- countries[rowSums(trade_matrix) > 0]

# Remove countries with no outgoing trade
trade_matrix <- trade_matrix[
  active_countries,
  active_countries,
  drop = FALSE
]

# Normalize
trade_matrix <- trade_matrix / rowSums(trade_matrix)

# Check matrix requirements. SHould all be TRUE, and rowSums be 1.
is.matrix(trade_matrix) 
nrow(trade_matrix) == ncol(trade_matrix) 
identical(
  rownames(trade_matrix),
  colnames(trade_matrix) 
)
rowSums(trade_matrix) 

# These should be FALSE
anyNA(trade_matrix)
any(trade_matrix < 0)

# Save
write.csv(
  trade_matrix,
  "output/trade_matrix_daily_probabilities_eu.csv",
  row.names = TRUE)

write.csv(
  complete_data,
  "output/trade_movement_rates_eu.csv",
  row.names = FALSE)
