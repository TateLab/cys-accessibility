### Reformat data for analysis - Figure 2 (all cysteines found across 12 datasets + uniprot Cys annotations)

library(tidyverse)
library(fs)
library(data.table)
library(magrittr)
library(janitor)
library(readxl)
library(here)


### Set absolute paths to data directories
path.to.data.folder <- here("data") ## Set to wherever the directory of the unzipped folder of raw/resource/processed data 

output.file.date <- gsub("-", "", Sys.Date())


## Read in human fasta
res.fasta <- dir_ls(here("resources"), regexp = "FASTA") %>%
  fread() %>%
  mutate(protein_id = sub(Gene, 
                          pattern = "..\\|(.*)\\|.*",
                          replacement = "\\1"))

#### Human AlphaFold predictions
res.af <- dir_ls(here("resources"), regexp = "AlphaFoldPredicted.*hsapiens") %>% 
  fread() %>%
  rename(pPSE = "nAA_12_70_pae")



####################################################################################
######################## UniProt Cys PTM assignments ###############################
####################################################################################



#Uniprot PTM annotations
res.uniprot.annotations <- dir_ls(here("resources"), regexp = "allann") %>%
  fread()

### Make data frame for each cys modification

### Nucleophiles
uniprot.nucleophile <- res.uniprot.annotations %>%
  filter(str_detect(`Active site`, "Nucleophile")) %>%
  mutate(position = as.integer(sub(`Active site`, 
                                   pattern = "ACT_SITE (\\d*)[^\\d].*",
                                   replacement = "\\1"))) %>%
  rename(protein_id = "Entry") %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  select(protein_id, quality, AA, position, pPSE, structure_group) %>%
  filter(AA == "C",
         !is.na(pPSE),
         !is.na(position)) %>%
  mutate(Class = "Active Site Nucleophile") %>%
  distinct()


#### All active site residues (not just cysteines)
uniprot.activesite <- res.uniprot.annotations %>%
  filter(!is.na(`Active site`)) %>%
  mutate(position = as.integer(sub(`Active site`, 
                                   pattern = "ACT_SITE (\\d*)[^\\d].*",
                                   replacement = "\\1"))) %>%
  rename(protein_id = "Entry") %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  select(protein_id, quality, AA, position, pPSE, structure_group) %>%
  filter(AA == "C",
         !is.na(position),
         !is.na(pPSE)) %>%
  mutate(Class = "All Active Site") %>%
  distinct()


### SwissPalm
swiss.palmitoyl <- here("resources") %>%
  dir_ls(regexp="wiss") %>%
  fread() %>% 
  clean_names() %>%
  filter(str_detect(organism, "sapiens"),
         acyl_type == "palmitate",
         uncertain_pos == "No") %>%
  select(uniprot_ac, pos) %>%
  rename(protein_id = "uniprot_ac",
         position = "pos") %>%
  left_join(res.af, by=c("protein_id", "position")) %>%
  select(protein_id, quality, AA, position, pPSE, structure_group) %>%
  filter(!is.na(pPSE),
         AA == "C",
         !is.na(position)) %>%
  mutate(Class = "Palmitoyl - SwissPalm") %>%
  distinct()


### Palmitoylated cysteines
uniprot.palmitoyl <- res.uniprot.annotations %>%
  mutate(All.Residues = str_extract_all(Lipidation, 'LIPID (\\d*)[^\\d]* /note="S-palmitoyl cysteine'),
         All.Residues = str_extract_all(as.character(All.Residues), 'LIPID (\\d*)'),
         All.Residues = sub(All.Residues, 
                            pattern = 'c\\(',
                            replacement = ""),
         All.Residues = gsub(All.Residues, 
                             pattern = '[\\(\\;\\"\\)]',
                             replacement = ""),
         All.Residues = gsub(All.Residues, 
                             pattern = 'LIPID ',
                             replacement = "")) %>%
  separate(All.Residues, into = paste("Residue", 1:6), sep = ", ") %>%
  select(c("Entry", "Gene Names", paste("Residue", 1:6))) %>%
  filter(`Residue 1` != "character0") %>%
  pivot_longer(cols = 3:8,
               names_to = "Residue.X",
               values_to = "position") %>%
  select(-Residue.X, -`Gene Names`) %>%
  rename(protein_id = "Entry") %>%
  mutate(position = as.integer(position)) %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  select(protein_id, quality, AA, position, pPSE, structure_group) %>%
  filter(!is.na(pPSE),
         AA == "C",
         !is.na(position)) %>%
  mutate(Class = "Palmitoyl - UniProt") %>%
  distinct()


### Farnesylated proteins
uniprot.prenyl <- res.uniprot.annotations %>%
  filter(str_detect(Lipidation, 
                    pattern = "farnesyl") | 
           str_detect(Lipidation, pattern = "geranyl")) %>%
  mutate(All.Residues.F = str_extract_all(Lipidation, 'LIPID (\\d*)[^\\d]* /note="S-farnesyl cysteine'),
         All.Residues.F = str_extract_all(as.character(All.Residues.F), 'LIPID (\\d*)'),
         All.Residues.GG = str_extract_all(Lipidation, 'LIPID (\\d*)[^\\d]* /note="S-geranylgeranyl cysteine'),
         All.Residues.GG = str_extract_all(as.character(All.Residues.GG), 'LIPID (\\d*)'),
         All.Residues.F = sub(All.Residues.F, 
                              pattern = 'c\\(',
                              replacement = ""),
         All.Residues.F = gsub(All.Residues.F, 
                               pattern = '[\\(\\;\\"\\)]',
                               replacement = ""),
         All.Residues.F = gsub(All.Residues.F, 
                               pattern = 'LIPID ',
                               replacement = ""),
         All.Residues.GG = sub(All.Residues.GG, 
                               pattern = 'c\\(',
                               replacement = ""),
         All.Residues.GG = gsub(All.Residues.GG, 
                                pattern = '[\\(\\;\\"\\)]',
                                replacement = ""),
         All.Residues.GG = gsub(All.Residues.GG, 
                                pattern = 'LIPID ',
                                replacement = "")) %>%
  separate(All.Residues.GG, into = paste0("All.Residues.GG.", 1:2), sep = ", ") %>%
  select(Entry, `Gene Names`, 24:26) %>%
  pivot_longer(cols = 3:5,
               names_to = "Residue.X",
               values_to = "position") %>%
  filter(position != "character0",
         !is.na(position)) %>%
  mutate(position = as.numeric(position)) %>%
  rename(protein_id = "Entry") %>%
  select(protein_id, position) %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  select(protein_id,quality,  AA, position, pPSE, structure_group) %>%
  filter(!is.na(pPSE),
         !is.na(position),
         AA == "C") %>%
  mutate(Class = "Prenyl") %>%
  distinct()

## Metal binding cysteines
uniprot.metalcys <- res.uniprot.annotations %>%
  filter(!is.na(`Metal binding`)) %>%
  mutate(All.Residues = str_extract_all(`Metal binding`, 'METAL (\\d*)[^\\d]* /note=".*?"'),
         #All.Metals = str_extract_all(`Metal binding`, 'note=\\".*\\"\\; \\/')) %>%
         All.Residues = gsub(All.Residues, 
                             pattern = '[\\(\\;\\"\\)]',
                             replacement = "")) %>%
  separate(All.Residues, into = paste("Residue", 1:34), sep = ", ") %>%
  select(c("Entry", "Gene Names", paste("Residue", 1:34))) %>%
  filter(`Residue 1` != "character0") %>%
  pivot_longer(cols = 3:36,
               names_to = "Residue.X",
               values_to = "Residue") %>%
  filter(!is.na(Residue)) %>%
  select(-Residue.X) %>%
  mutate(position = as.integer(sub(Residue, 
                                   pattern = ".*METAL (\\d.*) /note=.*",
                                   replacement = "\\1"))) %>%
  rename(protein_id = "Entry") %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  filter(!is.na(pPSE),
         AA == "C",
         !is.na(position)) %>%
  select(protein_id,quality, AA, position, pPSE, structure_group) %>%
  mutate(Class = "Metal Binding") 

  
## Disulfides
uniprot.disulfide <- res.uniprot.annotations %>%
  filter(str_detect(`Disulfide bond`, "DISULF")) %>%
  select(Entry, `Gene Names`, `Disulfide bond`) %>%
  separate(`Disulfide bond`, into = paste0("DS",1:160), sep = "DISULFID") %>%
  select(-DS1) %>%
  pivot_longer(cols=3:161,
               names_to = "Type",
               values_to = "Info") %>%
  filter(!is.na(Info)) %>%
  separate(Info, into = c("Residues", "Annotation"), sep = '; /note="') %>%
  mutate(Residues.Sub = sub(Residues, pattern = "; /evidence=.*", replacement = ""),
         Class = case_when(str_detect(Annotation, "edox") ~ "Redox-active",
                           str_detect(Annotation, "nterchain") ~ "Interchain"),
         Class = ifelse(is.na(Class), 
                        yes = "Intrachain",
                        no = Class)) %>%
  separate(Residues.Sub, into = paste0("Residue",1:2), sep = "\\.\\.") %>%
  pivot_longer(cols = c("Residue1", "Residue2"),
               names_to = "Res",
               values_to = "position") %>%
  mutate(position = sub(position, pattern = " ", replacement = ""),
         position = as.numeric(position)) %>%
  filter(!is.na(position)) %>%
  rename(protein_id = "Entry") %>%
  select(protein_id, position, Class) %>% 
  distinct() %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  filter(!is.na(pPSE),
         !is.na(AA),
         AA == "C") %>%
  select(protein_id, quality, AA, position, pPSE, Class, structure_group) %>%
  mutate(Class = paste0(Class, " disulfide")) %>%
  distinct()



## Electrophilic metabolite binding Cys
hne.match <- dir_ls("resources", regexp = "HNE") %>%
  fread() %>%
  rename(uniprot_name = "From",
         protein_id = "Entry") %>%
  dplyr::select(uniprot_name, protein_id) %>%
  distinct()

hne.cys <- dir_ls(here("resources"), regexp = "/41592_2014_BFnmeth2759_MOESM638_ESM") %>%
  readxl::read_excel() %>%
  left_join(hne.match) %>%
  left_join(res.fasta) %>%
  na.omit() %>%
  mutate(Cys.Peptide.Position = str_locate(sequence, "\\*")[,1] - 1,
         Peptide.Sequence = gsub(sequence, pattern="\\*", replacement=""),
         Peptide.Prot.Position = str_locate(Sequence, Peptide.Sequence)[,1],
         position = Cys.Peptide.Position + Peptide.Prot.Position - 1,
         ID = paste0(protein_id, "_", position)) %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  filter(!is.na(pPSE),
         !is.na(AA),
         AA == "C") %>%
  select(protein_id, quality, AA, position, pPSE, structure_group) %>%
  mutate(Class = "Electrophilic metabolite-modified") %>%
  distinct()


## Combine to one dataframe for plotting
uniprot.combined <- rbind(uniprot.nucleophile,
                          uniprot.palmitoyl,
                          swiss.palmitoyl,
                          uniprot.metalcys,
                          uniprot.activesite,
                          uniprot.prenyl,
                          uniprot.disulfide) %>%
  add_count(Class, name = "Total") %>%
  mutate(Facet = paste0(Class, "\nn = ", Total))



write_csv(uniprot.combined, paste0(here("formatted_data"), "/", output.file.date, "_figure2_uniprot.csv"))
