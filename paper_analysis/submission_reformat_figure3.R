### Reformat data for analysis - Figure 2 (all cysteines found across 16 datasets + uniprot Cys annotations)

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



### Import each coverage dataset, parse protein_id (UniProt accession) and Cys position

### Yang et al., 2022, Curr Res Chem Bio (TOP-ABPP)
### IA-alkyne - 100uM - Lysate
reformat.yang.crcb <- dir_ls(here("data", "yang_crcb_2022"), regexp = "xls") %>%
  readxl::read_excel(sheet = 2) %>%
  filter(!str_detect(Protein, ","),## Remove proteins with multiple mapped accessions
         `UK-alkyne` == "+") %>% ## Unsaturated ketone
  rename(protein_id = "Protein") %>%
  left_join(res.fasta %>% select(Sequence, protein_id)) %>%
  mutate(Stripped.Sequence = sub(Peptide, 
                                 pattern = "([A-Z]*)\\*([A-Z]*)", 
                                 replacement = "\\1\\2"),
         Peptide.Cys = str_locate(Peptide, ## Location of Cys in peptide
                                  "C")[,1] - 1, # account for +1 of matching
         Peptide.Start = str_locate(string = Sequence, ## find peptide location in protein
                                    pattern = Stripped.Sequence)[,1], ## [1] takes start of seq
         position = Peptide.Start + Peptide.Cys,
         Dataset = "a,b-unsaturated ketone",
         Probe = "UK-alkyne",
         Year = "2022", 
         Concentration = "100uM",
         Reference = "Yang et al.", 
         Labelling = "Lysate",
         N.Matches.In.Protein = str_count(Sequence, pattern = Stripped.Sequence)) %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct()

### Lai et al., 2022, ChemRxiv
### NAIA - 10um - Lysate
reformat.lai.chemrxiv <- dir_ls(here("data", "lai_chemrxiv_2022"), regexp = "/supp") %>%
  readxl::read_excel(sheet =3) %>%
  select(1) %>%
  na.omit() %>%
  separate(1, into = c("Genes", "position"), sep = "_") %>%
  left_join(res.fasta) %>%
  mutate(position = as.numeric(sub(position, pattern="C",replacement="")),
         Dataset = "N-acryloylindole",
         Probe = "NAIA",
         Reference = "Lai et al.",
         Concentration = "10uM",
         Year = "2022",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct()

### Kemper et al., 2022, Nature Methods
### IA-DTB - 100uM - Lysate
reformat.kemper.natmethods <- dir_ls(here("data", "kemper_natmethods_2022"), regexp = "1398") %>%
  readxl::read_excel(sheet = 2) %>%
  filter(!is.na(...3)) %>%
  row_to_names(1) %>%
  separate(gene_res, into = c('gene', 'position'), sep = "_") %>%
  select(accession, position) %>%
  rename(protein_id = "accession") %>%
  mutate(position = as.numeric(position),
         Dataset = "TMT-ABPP",
         Probe = "IA-DTB",
         Year = "2022", 
         Concentration = "100uM",
         Reference = "Kemper et al.", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct()
  
  
### Yang et al., 2022, JACS (DIA-ABPP)
### IA-alkyne - 100uM - Lysate
reformat.yang.jacs <- dir_ls(here("data", "yang_jacs_2022"), regexp = "2022/ja") %>%
  readxl::read_excel(sheet = 2, na = "--") %>%
  filter(!str_detect(Proteins, ",")) %>% ## Remove proteins with multiple mapped accessions
  rename(protein_id = "Proteins") %>%
  left_join(res.fasta %>% select(Sequence, protein_id)) %>%
  mutate(Stripped.Sequence = gsub(Peptides, 
                                 pattern = "([A-Z]*)\\*([A-Z]*)", 
                                 replacement = "\\1\\2"),
         Peptide.Cys = str_locate(Peptides, ## Location of Cys in peptide
                                  "\\*")[,1] - 1, # account for +1 of matching
         Peptide.Start = str_locate(string = Sequence, ## find peptide location in protein
                                    pattern = Stripped.Sequence)[,1], ## [1] takes start of seq
         position = Peptide.Start + Peptide.Cys - 1,
         Dataset = "DIA-ABPP",
         Probe = "IA-alkyne",
         Year = "2022", 
         Concentration = "100uM",
         Reference = "Yang et al.", 
         Labelling = "Lysate",
         N.AA = nchar(Stripped.Sequence),
         N.Matches.In.Protein = str_count(Sequence, pattern = Stripped.Sequence)) %>%
  select(protein_id, position, 
         #Stripped.Sequence, 
         Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  #na.omit() %>%
  distinct()

  
  
### Kuljanin et al., 2021, Nature Biotechnology (SLC-ABPP)
### DBIA - 500uM - Lysate
reformat.kuljanin.natbiotech <- rbind(
  ## Converted to CSVs to remove slow loading by readxl
  dir_ls(here("data", "kuljanin_natbiotech_2021"), regexp = "41587_2020_778_MOESM9_ESM")[1] %>%
    fread() %>%
    mutate(Cell.Line = "HEK293T"),
  
  dir_ls(here("data", "kuljanin_natbiotech_2021"), regexp = "41587_2020_778_MOESM10_ESM")[1] %>%
    fread() %>%
    mutate(Cell.Line = "PaTu-8988T"),
  
  dir_ls(here("data", "kuljanin_natbiotech_2021"), regexp = "41587_2020_778_MOESM8_ESM")[1] %>%
    fread() %>%
    mutate(Cell.Line = "HCT116")) %>%
  mutate(protein_id = sub(`Uniprot ID`, 
                          pattern = "..\\|(.*)\\|.*",
                          replacement = "\\1")) %>%
  rename(position = "Site Position") %>%
  mutate(Dataset = "SLC-ABPP",
         Probe = "DBIA",
         Concentration = "500uM",
         Reference = "Kuljanin et al.",
         Year = "2021", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 

### Yan et al., ChemMedChem, 2021
### IA-alkyne - Unknown conc. - Lysate (2M Urea)
reformat.yan.cbc <- dir_ls(here("data", "yan_cbc_2021"), regexp = "cbic202000870") %>%
  readxl::read_excel(sheet=2) %>% 
  separate(identifier, into = c("protein_id", "position")) %>%
  filter(Backus == 1,
         !str_detect(protein_id, "ntaminant")) %>% ## Remove contaminants
  mutate(Dataset = "SP3-FAIMS",
         Probe = "IA-alkyne",
         Year = "2021", 
         Concentration = "NA",
         Reference = "Yan et al.",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 

### Cao et al., Anal. Chem., 2021
### Iodophenyl-IA - 200uM - Lysate
reformat.cao.analchem <- dir_ls(here("data", "cao_ac_2021"), regexp = "xls") %>%
  readxl::read_excel(sheet = "Suzuki_CuAAC_w_sp3") %>%
  separate(identifier, into = c("protein_id", "position"), sep = "_") %>%
  mutate(position = as.numeric(sub(position, pattern = "C", replacement = ""))) %>%
  filter(`identified_in_suzuku_dataset (0: no | 1: yes)` == 1) %>% ## Only Suzuki-identified Cys
  mutate(Dataset = "mCSCP",
         Probe = "Iodophenyl-IA",
         Year = "2021", 
         Concentration = "200uM",
         Reference = "Cao et al.", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 

### Zanon et al., ChemRxiv, 2021
### Ethynylbenziodoxolone (EBX2) - 100uM - Lysate
reformat.zanon.chemrxiv <- dir_ls(here("data", "zanon_chemrxiv_2021/"), regexp = "table-4") %>%
  readxl::read_excel(sheet = "HS EBX2-alkyne") %>%
  rename(protein_id = "UniProt code") %>%
  left_join(res.fasta, by = "protein_id") %>%
  filter(!is.na(protein_id)) %>%
  mutate(Peptide.Cys = str_locate(string = `Modified peptide`,
                                  pattern = "\\*")[,1] - 2,
         Peptide.Start = str_locate(string = Sequence, ## find peptide location in protein
                                    pattern = `Peptide sequence`)[,1], ## [1] takes start of seq
         position = Peptide.Start + Peptide.Cys,
         N.Matches.In.Protein = str_count(Sequence, pattern = `Peptide sequence`)) %>%
  mutate(Dataset = "Ethynylbenziodoxolone",
         Probe = "EBX2",
         Year = "2021", 
         Concentration = "100uM",
         Reference = "Zanon et al.",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 


  
### Vinogradova et al., Cell, 2020
### IA-DTB (100uM) - Lysate
reformat.vinogradova.cell <- dir_ls(here("data", "vinogradova_cell_2020"), regexp = "/1.*xlsx") %>%
  readxl::read_excel(sheet = 4) %>%
  filter(!is.na(...3)) %>%
  row_to_names(row_number = 1) %>%
  select(uniprot, residue_all) %>%
  separate(residue_all, into = paste0("residue_", 1:54)) %>% #### Separate reisdues and pivot longer for 1 row each
  pivot_longer(cols = 2:55,
               names_to = "residue",
               values_to = "position") %>%
  rename(protein_id = "uniprot") %>%
  mutate(position = as.numeric(sub(position, pattern = "C(\\d*)", replacement = "\\1")),
         Dataset = "TMT-ABPP",
         Probe = "IA-DTB",
         Year = "2020", 
         Concentration = "100uM",
         Reference = "Vinogradova et al.", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 


### Abo et al., ChemMedChem, 2016
### Caged-iodoketone (CIK4) - 200uM - Live cell
reformat.abo.cbc.cik4 <- dir_ls(here("data", "abo_cbc_2016"), regexp = "cbic") %>%
  readxl::read_excel(sheet = "CIK4_AllReps") %>%
  row_to_names(row_number = 1) %>%
  rename(protein_id = "id",
         position = "cysteine position") %>%
  filter(!str_detect(protein_id, "Reverse"),
         !str_detect(position, "Bad ID")) %>%
  mutate(position = as.numeric(position),
         Dataset = "Caged-iodoketone",
         Probe = "CIK4", 
         Year = "2016", 
         Concentration = "200uM",
         Reference = "Abo et al.",
         Labelling = "Live cell") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 


### Abo et al., ChemMedChem, 2016
### Caged-bromoketone (CBK1) - 200uM - Live cell
reformat.abo.cbc.cbk1 <- dir_ls(here("data", "abo_cbc_2016"), regexp = "cbic") %>%
  readxl::read_excel(sheet = "CBK1_AllReps") %>%
  row_to_names(row_number = 1) %>%
  rename(protein_id = "id",
         position = "cysteine position") %>%
  filter(!str_detect(protein_id, "Reverse"),
         !str_detect(position, "Bad ID")) %>%
  mutate(position = as.numeric(position),
         Dataset = "Caged-bromoketone",
         Probe = "CBK1", 
         Year = "2016", 
         Concentration = "200uM",
         Reference = "Abo et al.",
         Labelling = "Live cell") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 


### Abo et al., JACS, 2015
### Bromoketone (BK1) - 100uM - Lysate
reformat.abo.cbc.bk1 <- dir_ls(here("data", "abo_jacs_2015"), regexp = "/ja5b04350_si_001") %>%
  readxl::read_excel(sheet = 1) %>%
  row_to_names(row_number = 1) %>%
  rename(protein_id = "ID", 
         position = "Cysteine position") %>%
  filter(!str_detect(protein_id, "Reverse"),
         !str_detect(position, "Bad ID")) %>%
  mutate(position = as.numeric(position),
         Dataset = "Bromoketone",
         Probe = "BK1", 
         Year = "2015",
         Concentration = "100uM",
         Reference = "Abo et al.",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 


### Backus et al., Nature, 2016
### IA-alkyne - 100uM - Lysate
reformat.backus.nature <- dir_ls(here("data", "backus_nature_2016"), regexp = "/41586.*xlsx") %>%
  readxl::read_excel(sheet = 2) %>%
  select(1) %>%
  separate(1, into = c("protein_id", "position")) %>%
  mutate(position = as.integer(sub(position, pattern = "C", replacement = "")),
         Dataset = "isoTOP-ABPP",
         Probe = "IA-alkyne",
         Reference = "Backus et al.", 
         Year = "2016",
         Concentration = "100uM",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 

### Weerapana et al., Nature, 2010
### IA-alkyne - 100uM - Lysate
reformat.weerapana.nature <- dir_ls(here("data", "weerapana_nature_2010"), regexp = "S4") %>%
  readxl::read_excel() %>%
  rename(protein_id = "UniProt",
         position = "Residue") %>%
  filter(!str_detect(position, ",")) %>%
  mutate(position = as.integer(sub(position, pattern = "C", replacement = "")),
         Dataset = "isoTOP-ABPP",
         Probe = "IA-alkyne", 
         Concentration = "100uM",
         Reference = "Weerapana et al.",
         Year = "2010",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 


### Darabedian et al., Nature Chem. Bio., 2023
### IA-Phos (CPT) - 10mM - Lysate
reformat.darabedian.nchembio <- dir_ls(here("data", "darabedian_nchembio_2023"), regexp="/41589_2023_1273_MOESM3_ESM") %>%
  readxl::read_excel(sheet="Table S1") %>%
  row_to_names(row_number=1) %>%
  clean_names() %>%
  rename(protein_id = "uniprot_id", 
         position = "site") %>%
  filter(!str_detect(position, ";")) %>%
  mutate(Dataset = "CPT (CKi)",
         Probe = "IA-Phos", 
         Concentration = "10mM",
         Reference = "Darabedian et al.",
         Year = "2023",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct() 

### Motiwala et al., 2020, JACS
reformat.motiwala.jacs <- dir_ls(here("data", "motiwala_jacs_2020"), regexp = "ja9b08831_si_003") %>%
  readxl::read_excel(sheet = "BT-D2_200uM") %>%
  filter(str_detect(Modifications, "BT-D2"),
         !str_detect(`Master Protein Accessions`, ";")) %>%
  mutate(N.Cys = sub(Modifications, 
                     pattern = ".*(\\d)xBT-D2.*",
                     replacement = "\\1"),
         Peptide.Cys = sub(Modifications, 
                           pattern = ".*xBT-D2 \\[(.*?)\\].*",
                           replacement = "\\1"),
         Peptide.Cys = sub(Peptide.Cys,pattern = "C", replacement = "")) %>%
  filter(N.Cys == 1,
         !str_detect(Peptide.Cys, "K"),
         str_detect(Peptide.Cys, "\\d")) %>%
  rename(protein_id = "Master Protein Accessions",
         Peptide.Sequence="Sequence") %>%
  left_join(res.fasta, by="protein_id") %>%
  mutate(Peptide.Start = str_locate(string = Sequence, ## find peptide location in protein
                                    pattern = Peptide.Sequence)[,1], ## [1] takes start of seq
         Peptide.Start = as.numeric(Peptide.Start),
         Peptide.Cys = as.numeric(Peptide.Cys),
         position = Peptide.Start + Peptide.Cys - 1) %>%
  mutate(position = as.numeric(sub(position, pattern="C",replacement="")),
         Dataset = "Methylsulfonylbenzothiazole",
         Probe = "BT-D2",
         Reference = "Motiwala et al.",
         Concentration = "200uM",
         Year = "2020",
         Labelling = "Live cell") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct()


### Abegg et al., 2021, JACS
reformat.abegg.jacs <- dir_ls(here("data", "abegg_jacs_2021"), regexp = "21/ja1c") %>%
  readxl::read_excel(sheet = "Table S8",
                     col_names = T) %>%
  row_to_names(row_number=1) %>%
  clean_names() %>%
  mutate(Peptide.Cys = sub(modifications, pattern = ".*Desthiobiotin \\[(.*)\\].*", replacement = "\\1"),
         Peptide.Cys = gsub(Peptide.Cys, pattern = "C", replacement = ""),
         Peptide.Cys = gsub(Peptide.Cys, pattern = " ", replacement = "")) %>%
  separate(Peptide.Cys, into = paste0("Peptide.Cys_", 1:2), sep=";") %>%
  pivot_longer(cols=c("Peptide.Cys_1", "Peptide.Cys_2"), 
               names_to = "temp",
               values_to = "Peptide.Cys") %>%
  filter(!is.na(Peptide.Cys),
         Peptide.Cys != "") %>%
  select(1,2,7) %>%
  rename(protein_id = 1,
         peptide_seq = 2) %>%
  left_join(res.fasta) %>%
  mutate(peptide_seq_reformat = sub(peptide_seq, pattern = "\\[[RK]\\]\\.([A-Z]*)\\..*", replacement = "\\1"),
         Peptide.Start = str_locate(string = Sequence, ## find peptide location in protein
                                    pattern = peptide_seq_reformat)[,1], ## [1] takes start of seq
         Peptide.Start = as.numeric(Peptide.Start),
         Peptide.Cys = as.numeric(Peptide.Cys),
         position = Peptide.Start + Peptide.Cys - 1) %>%
  select(protein_id, position) %>%
  mutate(position = as.numeric(sub(position, pattern="C",replacement="")),
         Dataset = "Tetrafluoroalkylbenziodoxole",
         Probe = "TFBX",
         Reference = "Abegg et al.",
         Concentration = "100uM",
         Year = "2021",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct()


### Liu et al., 2023, ACS Chem. Bio.
reformat.motiwala.jacs <- dir_ls(here("data", "motiwala_jacs_2020"), regexp = "ja9b08831_si_003") %>%
  readxl::read_excel(sheet = "BT-D2_200uM") %>%
  filter(str_detect(Modifications, "BT-D2"),
         !str_detect(`Master Protein Accessions`, ";")) %>%
  mutate(N.Cys = sub(Modifications, 
                     pattern = ".*(\\d)xBT-D2.*",
                     replacement = "\\1"),
         Peptide.Cys = sub(Modifications, 
                           pattern = ".*xBT-D2 \\[(.*?)\\].*",
                           replacement = "\\1"),
         Peptide.Cys = sub(Peptide.Cys,pattern = "C", replacement = "")) %>%
  filter(N.Cys == 1,
         !str_detect(Peptide.Cys, "K"),
         str_detect(Peptide.Cys, "\\d")) %>%
  rename(protein_id = "Master Protein Accessions",
         Peptide.Sequence="Sequence") %>%
  left_join(res.fasta, by="protein_id") %>%
  mutate(Peptide.Start = str_locate(string = Sequence, ## find peptide location in protein
                                    pattern = Peptide.Sequence)[,1], ## [1] takes start of seq
         Peptide.Start = as.numeric(Peptide.Start),
         Peptide.Cys = as.numeric(Peptide.Cys),
         position = Peptide.Start + Peptide.Cys - 1) %>%
  mutate(position = as.numeric(sub(position, pattern="C",replacement="")),
         Dataset = "Methylsulfonylbenzothiazole",
         Probe = "BT-D2",
         Reference = "Motiwala et al.",
         Concentration = "200uM",
         Year = "2020",
         Labelling = "Live cell") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct()


reformat.liu.acscb <- here("data", "liu_acscb_2023") %>%
  dir_ls() %>%
  readxl::read_excel(sheet = "Table S4-2") %>%
  filter(!is.na(`...2`)) %>%
  row_to_names(row_number=1) %>%
  clean_names() %>%
  mutate(protein_id = sub(protein_i_ds, 
                          pattern = "^..\\|(.*?)\\|.*",
                          replacement = "\\1")) %>%
  filter(!str_detect(protein_id, ";"),
         !str_detect(positions_within_proteins, ";")) %>%
  rename(position = "positions_within_proteins") %>%
  mutate(position = as.numeric(position),
         Dataset = "Nitrile oxide",
         Probe = "W1",
         Reference = "Liu et al.",
         Concentration = "500uM",
         Year = "2023",
         Labelling = "Live cell") %>%
  select(protein_id, position, Dataset, Probe, Reference, Concentration, Labelling, Year) %>%
  distinct()

reformat.combined <- rbind(
  reformat.abo.cbc.bk1, 
  reformat.abo.cbc.cbk1, 
  reformat.abo.cbc.cik4,
  reformat.backus.nature, 
  reformat.cao.analchem,
  reformat.kemper.natmethods,
  reformat.kuljanin.natbiotech,
  reformat.lai.chemrxiv,
  reformat.vinogradova.cell,
  reformat.yan.cbc,
  reformat.yang.jacs,
  reformat.weerapana.nature,
  reformat.liu.acscb,
  reformat.abegg.jacs,
  reformat.motiwala.jacs,
  reformat.darabedian.nchembio,
  reformat.yang.crcb,
  reformat.zanon.chemrxiv) %>%
  mutate(position = as.numeric(position)) %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  filter(!is.na(protein_id),
         !is.na(position),
         !is.na(pPSE),
         quality > 70,
         AA == "C") %>%
  select(protein_id, position, Dataset, Probe, Concentration, Reference, Year, Labelling, AA, quality, pPSE, structure_group) %>%
  distinct()

write_csv(reformat.combined, paste0(here("formatted_data"), "/", output.file.date, "_figure3_sites.csv"))



#### Without fltering for AlphaFold prediction or quality
reformat.unfiltered <- rbind(reformat.abo.cbc.bk1, 
                             reformat.abo.cbc.cbk1, 
                             reformat.abo.cbc.cik4,
                             reformat.backus.nature, 
                             reformat.cao.analchem,
                             reformat.kemper.natmethods,
                             reformat.abegg.jacs,
                             reformat.kuljanin.natbiotech,
                             reformat.lai.chemrxiv,
                             reformat.vinogradova.cell,
                             reformat.yan.cbc,
                             reformat.liu.acscb,
                             reformat.yang.jacs,
                             reformat.weerapana.nature,
                             reformat.motiwala.jacs,
                             reformat.darabedian.nchembio,
                             reformat.yang.crcb,
                             reformat.zanon.chemrxiv) %>%
  mutate(position = as.numeric(position)) %>%
  distinct()

write_csv(reformat.unfiltered, paste0(here("formatted_data"), "/", output.file.date, "_figure3_unfiltered-sites.csv"))
