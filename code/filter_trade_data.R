# BACI
#
# Version: 202601

# Release Date: 2026 01 22

# Weblink: http://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37

# Content:
#  Trade flows at the year - exporter - importer - product level.
#  Products in Harmonized System 6-digit nomenclature.
#  Values in thousand USD and quantities in metric tons.

# List of Variables:
#  t: year
#  i: exporter
#  j: importer
#  k: product
#  v: value
#  q: quantity

# Reference:
#  Gaulier, G. and Zignago, S. (2010)
#  BACI: International Trade Database at the Product-Level. The 1994-2007 Version.
#  CEPII Working Paper, N°2010-23

library(tidyr)
library(stringr)
library(dplyr)
library(readr)
library(recode)

# Load data and codebooks
data_2022 <- read_csv("data/trade_data/BACI_HS22_Y2022_V202601.csv")
data_2023 <- read_csv("data/trade_data/BACI_HS22_Y2023_V202601.csv")
data_2024 <- read_csv("data/trade_data/BACI_HS22_Y2024_V202601.csv") 

product_codes <- read_csv("data/trade_data/product_codes_HS22_V202601.csv")
country_codes <- read_csv("data/trade_data/country_codes_V202601.csv")

# Helpers
relabel_data <- function(data, country_codes, product_codes) {
  # Relabel data for better readability. Returns labelled dataframe
  
  data <- data |> 
    rename(
      year = t,
      exporter = i,
      importer = j,
      product = k,
      value = v,
      quantity = q
    ) |> 
    mutate(
      exporter = recode_values(
        exporter,
        from = country_codes$country_code,
        to = country_codes$country_name
      ),
      importer = recode_values(
        importer,
        from = country_codes$country_code,
        to = country_codes$country_name
      ))
  
  return(data)
}

# Live swine product codes are: 010310, 010391, 010392
data_2022 <- relabel_data(data_2022, country_codes, product_codes) |> 
  filter(
    product %in% c("010310", "010391", "010392")
  )

data_2023 <- relabel_data(data_2023, country_codes, product_codes) |> 
  filter(
    product %in% c("010310", "010391", "010392")
  )

data_2024 <- relabel_data(data_2024, country_codes, product_codes) |> 
  filter(
    product %in% c("010310", "010391", "010392")
  )

# Merge
all_data <- merge(data_2022, data_2023, all = TRUE) |> 
  merge(data_2024, all = TRUE)

# Quick check
if(sum(nrow(data_2022), nrow(data_2023),nrow(data_2024)) == nrow(all_data)) 
  print("Merged row sums equivalent.")

all_data <- all_data |> 
  left_join(product_codes,
            by = join_by(product == code)) |> 
  rename()






