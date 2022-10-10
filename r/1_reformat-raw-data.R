library(tidyverse)
library(fs)
library(readxl)
library(data.table)
library(janitor)

# Read in and re-format published datasets for analysis

## Naming conventions for datasets:
### Raw data from each supplementary: raw_*dataframe-id*
### Uniformly formatted dataframes: format_*dataframe-id*
# Columns for format_df:
column_names_reformat <- c("protein_id", ## Uniprot accession number
                           "AA", ## Amino acid
                           "position", ## Residue number (from structuremap output)
                           "HalfShell", ## pPSE value
                           "quality", ## AlphaFold residue prediction quality
                           "structure_group") ## Predicted secondary structure


### Weerapana/Abo optimised caged-probes (REF)
raw_caged.cik4 <- dir_ls("data/caged_weerapana", regexp = "cbic") %>%
  readxl::read_excel(sheet = "CIK4_AllReps") 

reformat_caged.cik4 <- raw_caged.cik4 %>%
  row_to_names(row_number = 1) %>%
  rename(protein_id = "id",
         position = "cysteine position") %>%
  select(protein_id, position) %>%
  distinct() %>%
  mutate(position = as.integer(position)) %>%
  mutate(Dataset = "Caged-iodoketone",
         Probe = "200uM CIK4", 
         Reference = "Ref.6",
         Labelling = "Live cell")

## Write function to check the AF and uniprot residue numbers -- resolve/remove conflicts

