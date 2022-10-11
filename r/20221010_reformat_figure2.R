### Reformat data for analysis - Figure 2 (fragment screening datasets)

## Read in human fasta
res.fasta <- dir_ls("resources", regexp = "FASTA") %>%
  fread() %>%
  mutate(protein_id = sub(Gene, 
                          pattern = "..\\|(.*)\\|.*",
                          replacement = "\\1")) %>%
  select(protein_id, Sequence)

#### Human AlphaFold predictions
res.af <- dir_ls("resources", regexp = "AlphaFoldPredicted.*hsapiens") %>% 
  fread() %>%
  rename(pPSE = "nAA_12_70_pae")


######################################################################
############# In vitro vs. in situ fragment treatment ################
######################################################################


## Fragments treated both in vitro and in situ
treatment.overlapping.conditions <- c(2,4,8,9,10,11,12,13,14,27,21,28,29,31,33,38,41,45,51,56)


### Read in vitro fragment treatment
treatment.invitro.raw <- dir_ls("data/backus_cys", regexp = "/41586.*xlsx") %>%
  readxl::read_excel(sheet = 1) %>% 
  select(Identifier, 3:88) %>% 
  pivot_longer(2:87, 
               names_to = "Col.Name",
               values_to = "R") %>%
  separate(Col.Name, into = c("Fragment", "Concentration", "Treatment", "Cell.Line"), sep = "_") %>%
  separate(Identifier, into = c("protein_id", "position"), sep = "_") %>%
  mutate(Treatment = "Liganded in vitro",
         position = as.numeric(sub(position,
                                   pattern = "C", replacement = "")),
         R = as.numeric(R)) %>%
  filter(Fragment %in% treatment.overlapping.conditions, ## Keep only fragments shared
         !str_detect(position, ","), ## Check no multiply annotated cys
         !is.na(R)) %>% 
  group_by(protein_id, position, Treatment) %>%
  summarise(R = max(R), .groups = "drop") %>% ## Take maximum R value per residue across fragments/cell lines
  mutate(Liganded = ifelse(R > 4, 
                           yes = 1,
                           no = 0)) %>%
  select(protein_id, position, Liganded, Treatment) %>% 
  distinct() %>%
  filter(Liganded == 1)


### Read in situ fragment treatment
treatment.insitu.raw <- dir_ls("data/backus_cys", regexp = "/41586.*xlsx") %>%
  readxl::read_excel(sheet = 7) %>% 
  select(Identifier, 3:27) %>% 
  pivot_longer(2:26, 
               names_to = "Col.Name",
               values_to = "R") %>%
  separate(Col.Name, into = c("Fragment", "Concentration", "Treatment", "Cell.Line"), sep = "_") %>%
  separate(Identifier, into = c("protein_id", "position"), sep = "_") %>%
  mutate(Treatment = "Liganded in situ",
         position = as.numeric(sub(position,
                                   pattern = "C", replacement = "")),
         R = as.numeric(R)) %>%
  filter(Fragment %in% treatment.overlapping.conditions, ## Keep only fragments shared
         !str_detect(position, ","), ## Check no multiply annotated cys
         !is.na(R)) %>% 
  group_by(protein_id, position, Treatment) %>%
  summarise(R = max(R), .groups = "drop") %>% ## Take maximum R value per residue across fragments/cell lines
  mutate(Liganded = ifelse(R > 4, 
                           yes = 1,
                           no = 0)) %>%
  select(protein_id, position, Liganded, Treatment) %>% 
  distinct() %>%
  filter(Liganded == 1)



### Read in all fragment screening data to find non-liganded cysteines
treatment.nonliganded.raw <- dir_ls("data/backus_cys", regexp = "/41586.*xlsx") %>%
  readxl::read_excel(sheet = 2) %>%
  select(Identifier, 3:88) %>% 
  pivot_longer(2:87, 
               names_to = "Col.Name",
               values_to = "R") %>%
  separate(Col.Name, into = c("Fragment", "Concentration", "Treatment", "Cell.Line"), sep = "_") %>%
  separate(Identifier, into = c("protein_id", "position"), sep = "_") %>%
  mutate(Treatment = "Not liganded in vitro",
         position = as.numeric(sub(position,
                                   pattern = "C", replacement = "")),
         R = as.numeric(R)) %>%
  filter(Fragment %in% treatment.overlapping.conditions, ## Keep only fragments shared
         !str_detect(position, ","), ## Check no multiply annotated cys
         !is.na(R)) %>% 
  group_by(protein_id, position, Treatment) %>%
  summarise(R = max(R), .groups = "drop") %>% ## Take maximum R value per residue across fragments/cell lines
  mutate(Liganded = ifelse(R > 4, 
                           yes = 1,
                           no = 0)) %>%
  select(protein_id, position, Liganded, Treatment) %>% 
  distinct() %>%
  filter(Liganded == 0)

### Combine in vitro and in situ fragments
treatment.combined <- rbind(treatment.invitro.raw,
                            treatment.insitu.raw,
                            treatment.nonliganded.raw) %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  filter(AA == "C",
         quality > 70,
         structure_group != "unstructured",
         !is.na(Liganded),
         !is.na(pPSE))

write_csv(treatment.combined, "formatted_data/20221010_figure2_treatment.csv")


### Import each fragment dataset, parse protein_id (UniProt accession), Cys position, fragment ID and R (competition) value

######################################################################
############### Scout fragments (KB02, KB05) #########################
######################################################################

### Reformat into tidy table with: protein_id, position, max R value (across two scout fragments)

## Vinogradova et al., 2020
scout.vinogradova.cell <- dir_ls("data/vinogradova_cell", regexp = "/1-.*S0092867420308230.*xlsx") %>%
  readxl::read_excel(sheet = 7) %>%
  select(2,4,24) %>%
  filter(!is.na(`...2`)) %>%
  row_to_names(row_number = 1) %>%
  rename(protein_id = "Uniprot",
         position = "Residues",
         R = "Max") %>%
  mutate(position = as.integer(sub(position, ## Also removes multiply annotated Cys, as cannot be coerced to integer
                                   pattern = "C", 
                                   replacement = "")), 
         R = as.numeric(R)) %>%
  filter(!is.na(R),
         !is.na(position)) %>%
  mutate(Dataset="Vinogradova et al., 2020\nisoTOP/TMT-ABPP") %>%
  distinct()

### Kuljanin et al., Nature Biotechnology 2021
scout.kuljanin.natbiotech <- dir_ls("data/gygi_cys", regexp = "/41587.*MOESM6.*xlsx") %>%
  readxl::read_excel(sheet = 2) %>%
  mutate(protein_id = sub(`Uniprot ID`, 
                          pattern = "..\\|(.*)\\|.*",
                          replacement = "\\1"),
         KB02 = as.numeric(`Comp ratio (KBO2)`),
         KB05 = as.numeric(`Comp ratio (KBO5)`)) %>%
  rename(position = "Site Position") %>% 
  mutate(R = pmax(KB02, KB05)) %>%
  select(protein_id, position, R) %>%
  filter(!is.na(R),
         !is.na(position)) %>%
  mutate(Dataset="Kuljanin et al., 2021\nSLC-ABPP") %>%
  distinct()

### Yang et al., 2021
scout.yang.jacs <- dir_ls("data/dia_abpp", regexp = "pp/ja") %>%
  readxl::read_excel(sheet = 2, na = "--") %>%
  filter(!str_detect(Proteins, ",")) %>% ## Multiple protein annotations removed
  rename(protein_id = "Proteins") %>%
  left_join(res.fasta %>% select(Sequence, protein_id)) %>%
  mutate(Stripped.Sequence = sub(Peptides, 
                                 pattern = "([A-Z]*)\\*([A-Z]*)", 
                                 replacement = "\\1\\2"),
         Peptide.Cys = str_locate(Peptides, ## Location of Cys in peptide
                                  "\\*")[,1] - 1, # account for +1 of matching
         Peptide.Start = str_locate(string = Sequence, ## find peptide location in protein
                                    pattern = Stripped.Sequence)[,1], ## [1] takes start of seq
         position = Peptide.Start + Peptide.Cys - 1) %>%
  mutate(R = pmax(F2, F5)) %>%
  select(protein_id, position, R) %>%
  filter(!is.na(R),
         !is.na(position)) %>%
  mutate(Dataset="Yang et al., 2022\nDIA-ABPP") %>%
  distinct()

  
### Backus et al., 2016
scout.backus.nature <- dir_ls("data/backus_cys", regexp = "/41586.*xlsx") %>%
  readxl::read_excel(sheet = 2) %>%
  select(1,3,4,9,10) %>% ## Fragments 2 and 5
  separate(1, into = c("protein_id", "position"), sep = "_") %>%
  filter(!str_detect(position, ",")) %>% ## Check no multiply annotated cys
  mutate(position = as.integer(sub(position, pattern = "C", replacement = "")),
         R = pmax(`2_500uM_invitro_231`, `2_500uM_invitro_ramos`,
                  `5_500uM_invitro_231`, `5_500uM_invitro_ramos`)) %>%
  select(protein_id, position, R) %>%
  filter(!is.na(R),
         !is.na(position)) %>%
  mutate(Dataset="Backus et al., 2016\nisoTOP-ABPP") %>%
  distinct()

scout.combined <- rbind(scout.backus.nature,
                        scout.kuljanin.natbiotech,
                        scout.vinogradova.cell,
                        scout.yang.jacs) %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  filter(!is.na(pPSE),
         !is.na(position),
         AA == "C",
         quality > 70,
         structure_group != "unstructured") %>%
  group_by(protein_id, position, Dataset) %>%
  filter(R == max(R)) %>%
  ungroup() %>%
  select(protein_id,position, R, Dataset, quality, structure_group, AA, pPSE) 
  

write_csv(scout.combined, "formatted_data/20221010_figure2_scout.csv")

######################################################################
######################## Fragment screenings #########################
######################################################################

### Manual annotation of chloroacetamide/acrylamide warheads from screened fragments
### Chloroacetamide/acrylamide annotation for each fragment
fragment.warhead.backus <- data.frame(
  Fragment = c(2:15,17,20:56),
  Warhead = c(1,1,1,2,0,1, ## Taken from Extended Data Figure 1 - Backus et al., Nature, 2016
              1,1,1,1,1,1,
              2,2,0,1,1,1,0,
              0,1,2,1,1,1,
              1,2,1,1,1,1,1,
              2,2,1,2,0,0,
              1,1,1,2,2,2,1,
              1,1,1,2,1,1,2)) %>% ## 1 = Cl, 2 = Ac, 0 = other
  mutate(Frag.Class = ifelse(Warhead == 2, 
                             yes = "Acrylamide",
                             no = ifelse(Warhead == 1, 
                                         yes= "Chloroacetamide",
                                         no = "Other"))) %>%
  select(Fragment, Frag.Class) 

fragment.warhead.yang <- data.frame(Fragment = c("F2","F3","F4","F5","F7","F8","F9","F10","F11",
                                                 "F12","F13","F14","F20","F21","F23","F27","F28",
                                                 "F30","F31","F32", "F33","F38","F52","F56")) %>%
  mutate(Frag.Class = ifelse(Fragment %in% paste0("F", c(5,14,23,31,38,56)),
                             yes = "Acrylamide",
                             no = "Chloroacetamide")) %>% ## Taken from Supplementary Information - Yang et al., JACS, 2021
  select(Fragment, Frag.Class)


### Import each fragment dataset, parse protein_id, Cys position, fragment ID and R (competition) value

## Vinogradova et al., 2020
fragment.vinogradova.cell <- dir_ls("data/vinogradova_cell", regexp = "/1-.*S0092867420308230.*xlsx") %>%
  readxl::read_excel(sheet = 7) %>%
  select(2,4,11:23) %>%
  filter(!is.na(`...2`)) %>%
  row_to_names(row_number = 1) %>%
  pivot_longer(cols = 3:15,
               names_to = "Fragment",
               values_to = "R") %>%
  rename(protein_id = "Uniprot",
         position = "Residues") %>%
  mutate(position = as.integer(sub(position, ## Also removes multiply annotated Cys, as cannot be coerced to integer
                                   pattern = "C", 
                                   replacement = "")), 
         R = as.numeric(R)) %>%
  filter(!is.na(R),
         !is.na(position)) %>%
  mutate(Dataset = "Vinogradova et al., 2020\nisoTOP/TMT-ABPP",
         Frag.Class = "Elaborated electrophile") %>%
  group_by(protein_id, position, Fragment) %>%
  filter(R == max(R)) %>% ## Keep only single R value for protein_id, position, Fragment combinations
  ungroup() 


## Kuljanin et al., 2021
fragment.kuljanin.natbiotech <- rbind(
  ## Converted to CSVs to remove slow loading by readxl
  dir_ls("data/gygi_cys", regexp = "41587_2020_778_MOESM9_ESM")[1] %>%
    fread() %>%
    mutate(Cell.Line = "HEK293T"),
  
  dir_ls("data/gygi_cys", regexp = "41587_2020_778_MOESM10_ESM")[1] %>%
    fread() %>%
    mutate(Cell.Line = "PaTu-8988T"),
  
  dir_ls("data/gygi_cys", regexp = "41587_2020_778_MOESM8_ESM")[1] %>%
    fread() %>%
    mutate(Cell.Line = "HCT116")) %>%
  mutate(protein_id = sub(`Uniprot ID`, 
                          pattern = "..\\|(.*)\\|.*",
                          replacement = "\\1")) %>%
  rename(position = "Site Position") %>%
  select(protein_id, position, 6:290) %>%
  pivot_longer(cols = 3:287,
               names_to = "Fragment",
               values_to = "R") %>%
  mutate(position = as.integer(sub(position, ## Also removes multiply annotated Cys, as cannot be coerced to integer
                                   pattern = "C", 
                                   replacement = "")), 
         R = as.numeric(R),
         Frag.Class = ifelse(str_detect(Fragment, "CL"),
                             yes = "Chloroacetamide",
                             no = ifelse(str_detect(Fragment, "AC"),
                                         yes = "Acrylamide",
                                         no = "Other"))) %>%
  filter(!is.na(R),
         !is.na(position),
         !is.na(Frag.Class)) %>%
  mutate(Dataset = "Kuljanin et al., 2021\nSLC-ABPP") %>%
  distinct() %>%
  group_by(protein_id, position, Fragment, Frag.Class) %>%
  filter(R == max(R)) %>% ## Keep highest ligandability value for each fragment
  ungroup() ## across 3 cell lines profiles


### Yang et al., 2022
fragment.yang.jacs <- dir_ls("data/dia_abpp", regexp = "pp/ja") %>%
  readxl::read_excel(sheet = 2, na = "--") %>%
  filter(!str_detect(Proteins, ",")) %>%
  rename(protein_id = "Proteins") %>%
  left_join(res.fasta %>% select(Sequence, protein_id)) %>%
  mutate(Stripped.Sequence = sub(Peptides, 
                                 pattern = "([A-Z]*)\\*([A-Z]*)", 
                                 replacement = "\\1\\2"),
         Peptide.Cys = str_locate(Peptides, ## Location of Cys in peptide
                                  "\\*")[,1] - 1, # account for +1 of matching
         Peptide.Start = str_locate(string = Sequence, ## find peptide location in protein
                                    pattern = Stripped.Sequence)[,1], ## [1] takes start of seq
         position = Peptide.Start + Peptide.Cys - 1) %>%
  select(protein_id, position, 3:26) %>%
  pivot_longer(cols = 3:26,
               names_to = "Fragment",
               values_to = "R") %>%
  mutate(position = as.integer(sub(position, ## Also removes multiply annotated Cys, as cannot be coerced to integer
                                   pattern = "C", 
                                   replacement = "")), 
         R = as.numeric(R)) %>%
  filter(!is.na(R),
         !is.na(position)) %>%
  left_join(fragment.warhead.yang, by = "Fragment") %>%
  mutate(Dataset="Yang et al., 2022\nDIA-ABPP") %>%
  distinct() %>%
  group_by(protein_id, position, Fragment) %>%
  filter(R == max(R)) %>% ## Keep only single R value for protein_id, position, Fragment combinations
  ungroup() 


### Backus et al., 2016
fragment.backus.nature <- dir_ls("data/backus_cys", regexp = "/41586.*xlsx") %>%
  readxl::read_excel(sheet = 2) %>%
  select(1,3:88) %>%
  separate(1, into = c("protein_id", "position"), sep = "_") %>%
  filter(!str_detect(position, ",")) %>% ## Check no multiply annotated cys
  pivot_longer(cols = 3:88,
               values_to = "R",
               names_to = "Fragment") %>%
  mutate(position = as.integer(sub(position, pattern = "C", replacement = "")),
         R = as.numeric(R),
         Fragment = sub(Fragment, pattern = "(\\d*)_.*", replacement = "\\1"),
         Fragment = as.numeric(Fragment)) %>%
  distinct() %>%
  filter(!is.na(R),
         !is.na(position)) %>%
  left_join(fragment.warhead.backus, by = "Fragment") %>%
  mutate(Dataset="Backus et al., 2016\nisoTOP-ABPP") %>%
  distinct() %>%
  group_by(protein_id, position, Fragment, Frag.Class) %>%
  filter(R == max(R)) %>% ## Keep only single R value for protein_id, position, Fragment combinations (across 2 cell lines)
  ungroup() 



fragment.combined <- rbind(fragment.backus.nature, 
                           fragment.kuljanin.natbiotech, 
                           fragment.vinogradova.cell, 
                           fragment.yang.jacs) %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  filter(!is.na(pPSE),
         !is.na(position),
         AA == "C",
         quality > 70,
         structure_group != "unstructured") %>%
  group_by(protein_id, position, Fragment, Frag.Class, Dataset) %>%
  filter(R == max(R)) %>%
  ungroup() %>%
  select(protein_id, position, Fragment, Frag.Class, R, Dataset, AA, quality, structure_group, pPSE)

write_csv(fragment.combined, "formatted_data/20221010_figure2_fragment-all-f-c-interactions.csv")


######################################################################
######## Max ligandability per Cys from fragment screenings ##########
######################################################################

fragment.singlesite <- fragment.combined %>%
  group_by(protein_id, position, Dataset) %>%
  filter(R == max(R)) %>%
  ungroup()

write_csv(fragment.singlesite, "formatted_data/20221010_figure2_fragment_max-R-per-cys.csv")

  
  
  
  