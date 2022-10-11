### Reformat data for analysis - Figure 1 (all cysteines found across 12 datasets + uniprot Cys annotations)


## Read in human fasta
res.fasta <- dir_ls("resources", regexp = "FASTA") %>%
  fread() %>%
  mutate(protein_id = sub(Gene, 
                          pattern = "..\\|(.*)\\|.*",
                          replacement = "\\1"))

#### Human AlphaFold predictions
res.af <- dir_ls("resources", regexp = "AlphaFoldPredicted.*hsapiens") %>% 
  fread() %>%
  rename(pPSE = "nAA_12_70_pae")



### Import each coverage dataset, parse protein_id (UniProt accession) and Cys position

### Lai et al., 2022, ChemRxiv
### NAIA - 10um - Lysate
reformat.lai.chemrxiv <- dir_ls("data/nomura_naia", regexp = "/supp") %>%
  readxl::read_excel(sheet =3) %>%
  select(1) %>%
  na.omit() %>%
  separate(1, into = c("Genes", "position"), sep = "_") %>%
  left_join(res.fasta) %>%
  mutate(position = as.numeric(sub(position, pattern="C",replacement="")),
         Dataset = "N-acryloylindole",
         Probe = "10uM NAIA",
         Reference = "ChemRxiv, 2022",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct()

### Kemper et al., 2022, Nature Methods
### IA-DTB - 100uM - Lysate
reformat.kemper.natmethods <- dir_ls("data/kemper_phos/", regexp = "1398") %>%
  readxl::read_excel(sheet = 2) %>%
  filter(!is.na(...3)) %>%
  row_to_names(1) %>%
  separate(gene_res, into = c('gene', 'position'), sep = "_") %>%
  select(accession, position) %>%
  rename(protein_id = "accession") %>%
  mutate(position = as.numeric(position),
         Dataset = "TMT-ABPP",
         Probe = "100uM IA-DTB",
         Reference = "Nature Methods, 2022", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct()
  
  
### Yang et al., 2022, JACS (DIA-ABPP)
### IA-alkyne - 100uM - Lysate
reformat.yang.jacs <- dir_ls("data/dia_abpp", regexp = "pp/ja") %>%
  readxl::read_excel(sheet = 2, na = "--") %>%
  filter(!str_detect(Proteins, ",")) %>% ## Remove proteins with multiple mapped accessions
  rename(protein_id = "Proteins") %>%
  left_join(res.fasta %>% select(Sequence, protein_id)) %>%
  mutate(Stripped.Sequence = sub(Peptides, 
                                 pattern = "([A-Z]*)\\*([A-Z]*)", 
                                 replacement = "\\1\\2"),
         Peptide.Cys = str_locate(Peptides, ## Location of Cys in peptide
                                  "\\*")[,1] - 1, # account for +1 of matching
         Peptide.Start = str_locate(string = Sequence, ## find peptide location in protein
                                    pattern = Stripped.Sequence)[,1], ## [1] takes start of seq
         position = Peptide.Start + Peptide.Cys - 1,
         Dataset = "DIA-ABPP",
         Probe = "100uM IA-alkyne",
         Reference = "JACS, 2022", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct()


### Kuljanin et al., 2021, Nature Biotechnology (SLC-ABPP)
### DBIA - 500uM - Lysate
reformat.kuljanin.natbiotech <- rbind(
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
  mutate(Dataset = "SLC-ABPP",
         Probe = "500uM DBIA",
         Reference = "Nature Biotech., 2021", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct() 

### Yan et al., ChemMedChem, 2021
### IA-alkyne - Unknown conc. - Lysate (2M Urea)
reformat.yan.cmc <- dir_ls("data/yan_faims") %>%
  readxl::read_excel(sheet=2) %>% 
  separate(identifier, into = c("protein_id", "position")) %>%
  filter(Backus == 1,
         !str_detect(protein_id, "ntaminant")) %>% ## Remove contaminants
  mutate(Dataset = "SP3-FAIMS",
         Probe = "-- IA-alkyne",
         Reference = "ChemMedChem, 2021", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct() 

### Cao et al., Anal. Chem., 2021
### Iodophenyl-IA - 200uM - Lysate
reformat.cao.analchem <- dir_ls("data/backus_suzuki", regexp = "xls") %>%
  readxl::read_excel(sheet = "Suzuki_CuAAC_w_sp3") %>%
  separate(identifier, into = c("protein_id", "position"), sep = "_") %>%
  mutate(position = as.numeric(sub(position, pattern = "C", replacement = ""))) %>%
  filter(`identified_in_suzuku_dataset (0: no | 1: yes)` == 1) %>% ## Only Suzuki-identified Cys
  mutate(Dataset = "mCSCP",
         Probe = "200uM Iodophenyl-IA",
         Reference = "Anal. Chem., 2021", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct() 

### Zanon et al., ChemRxiv, 2021
### Benziodoxole (EBX1) - 100uM - Lysate
reformat.zanon.chemrxiv <- dir_ls("data/zanon_raw/", regexp = "table-4") %>%
  readxl::read_excel(sheet = "HS EBX2-alkyne") %>%
  rename(protein_id = "UniProt code") %>%
  left_join(res.fasta, by = "protein_id") %>%
  filter(!is.na(protein_id)) %>%
  mutate(Peptide.Cys = str_locate(string = `Modified peptide`,
                                  pattern = "\\*")[,1] - 2,
         Peptide.Start = str_locate(string = Sequence, ## find peptide location in protein
                                    pattern = `Peptide sequence`)[,1], ## [1] takes start of seq
         position = Peptide.Start + Peptide.Cys) %>%
  mutate(Dataset = "Benziodoxole",
         Probe = "100uM EBX2",
         Reference = "ChemRxiv, 2021",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct() 
  
### Vinogradova et al., Cell, 2020
### IA-DTB (100uM) - Lysate
reformat.vinogradova.cell <- dir_ls("data/vinogradova_cell", regexp = "/1.*xlsx") %>%
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
         Probe = "100uM IA-DTB",
         Reference = "Cell, 2020", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct() 


### Abo et al., ChemMedChem, 2016
### Caged-iodoketone (CIK4) - 200uM - Live cell
reformat.abo.cmc.cik4 <- dir_ls("data/caged_weerapana", regexp = "cbic") %>%
  readxl::read_excel(sheet = "CIK4_AllReps") %>%
  row_to_names(row_number = 1) %>%
  rename(protein_id = "id",
         position = "cysteine position") %>%
  filter(!str_detect(protein_id, "Reverse"),
         !str_detect(position, "Bad ID")) %>%
  mutate(position = as.numeric(position),
         Dataset = "Caged-iodoketone",
         Probe = "200uM CIK4", 
         Reference = "ChemMedChem, 2016",
         Labelling = "Live cell") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct() 


### Abo et al., ChemMedChem, 2016
### Caged-bromoketone (CBK1) - 200uM - Live cell
reformat.abo.cmc.cbk1 <- dir_ls("data/caged_weerapana", regexp = "cbic") %>%
  readxl::read_excel(sheet = "CBK1_AllReps") %>%
  row_to_names(row_number = 1) %>%
  rename(protein_id = "id",
         position = "cysteine position") %>%
  filter(!str_detect(protein_id, "Reverse"),
         !str_detect(position, "Bad ID")) %>%
  mutate(position = as.numeric(position),
         Dataset = "Caged-bromoketone",
         Probe = "200uM CBK1", 
         Reference = "ChemMedChem, 2016",
         Labelling = "Live cell") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct() 


### Abo et al., JACS, 2015
### Bromoketone (BK1) - 100uM - Lysate
reformat.abo.cmc.bk1 <- dir_ls("data/caged_weerapana", regexp = "/ja5b04350_si_001") %>%
  readxl::read_excel(sheet = 1) %>%
  row_to_names(row_number = 1) %>%
  rename(protein_id = "ID", 
         position = "Cysteine position") %>%
  filter(!str_detect(protein_id, "Reverse"),
         !str_detect(position, "Bad ID")) %>%
  mutate(position = as.numeric(position),
         Dataset = "Bromoketone",
         Probe = "100uM BK1", 
         Reference = "JACS, 2015",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct() 


### Backus et al., Nature, 2016
### IA-alkyne - 100uM - Lysate
reformat.backus.nature <- dir_ls("data/backus_cys", regexp = "/41586.*xlsx") %>%
  readxl::read_excel(sheet = 2) %>%
  select(1) %>%
  separate(1, into = c("protein_id", "position")) %>%
  mutate(position = as.integer(sub(position, pattern = "C", replacement = "")),
         Dataset = "isoTOP-ABPP",
         Probe = "100uM IA-alkyne",
         Reference = "Nature, 2016", 
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct() 

### Weerapana et al., Nature, 2010
### IA-alkyne - 100uM - Lysate
reformat.weerapana.nature <- dir_ls("data/weerapana_2010_raw", regexp = "S4") %>%
  readxl::read_excel() %>%
  rename(protein_id = "UniProt",
         position = "Residue") %>%
  filter(!str_detect(position, ",")) %>%
  mutate(position = as.integer(sub(position, pattern = "C", replacement = "")),
         Dataset = "isoTOP-ABPP",
         Probe = "100uM IA-alkyne", 
         Reference = "Nature, 2010",
         Labelling = "Lysate") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling) %>%
  distinct() 


reformat.combined <- rbind(
  reformat.abo.cmc.bk1, 
  reformat.abo.cmc.cbk1, 
  reformat.abo.cmc.cik4,
  reformat.backus.nature, 
  reformat.cao.analchem,
  reformat.kemper.natmethods,
  reformat.kuljanin.natbiotech,
  reformat.lai.chemrxiv,
  reformat.vinogradova.cell,
  reformat.weerapana.nature,
  reformat.yan.cmc,
  reformat.yang.jacs,
  reformat.zanon.chemrxiv) %>%
  mutate(position = as.numeric(position)) %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  filter(!is.na(protein_id),
         !is.na(position),
         !is.na(pPSE),
         quality > 70,
         AA == "C") %>%
  select(protein_id, position, Dataset, Probe, Reference, Labelling, AA, quality, pPSE, structure_group) %>%
  distinct()

write_csv(reformat.combined, "formatted_data/20221010_figure1_sites.csv")

####################################################################################
######################## UniProt Cys PTM assignments ###############################
####################################################################################


#Uniprot PTM annotations
res.uniprot.annotations <- dir_ls("resources", regexp = "allann") %>%
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
  mutate(Class = "Palmitoyl") %>%
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
  select(c("Entry", "Gene Names", 24:26)) %>%
  pivot_longer(cols = 3:5,
               names_to = "Residue.X",
               values_to = "position") %>%
  mutate(position = as.numeric(position)) %>%
  filter(!is.na(position)) %>%
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
                                   replacement = "\\1")),
         Metal = sub(Residue, 
                     pattern = ".*METAL \\d.* /note=([A-z]*)[^A-z].*",
                     replacement = "\\1"),
         Metal = sub(Metal, 
                     pattern = ".*METAL \\d.* /note=([A-z]*)$",
                     replacement = "\\1"),
         Metal = ifelse(Metal == "Cu",
                        yes = "Copper",
                        no = ifelse(Metal == "Fe",
                                    yes = "Iron",
                                    no = ifelse(Metal == "Zn",
                                                yes = "Zinc",
                                                no = ifelse(Metal == "Li",
                                                            yes = "Lithium",
                                                            no = ifelse(Metal == "Ni",
                                                                        yes = "Nickel",
                                                                        no = Metal)))))) %>%
  rename(protein_id = "Entry") %>%
  select(protein_id, position, Metal) %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  filter(!is.na(pPSE),
         AA == "C",
         !is.na(position)) %>%
  select(protein_id,quality,  AA, position, pPSE, Metal, structure_group) %>%
  mutate(Class = "Metal Binding") 



## Disulfides
uniprot.disulfide <- res.uniprot.annotations %>%
  filter(str_detect(`Disulfide bond`, "DISULF")) %>%
  select(Entry, `Gene Names`, `Disulfide bond`) %>%
  mutate(Start = sub(`Disulfide bond`, 
                     pattern = "DISULFID (\\d*)[^\\d].*",
                     replacement = "\\1"),
         End = sub(`Disulfide bond`, 
                   pattern = "DISULFID \\d*\\.\\.(\\d*).*",
                   replacement = "\\1")) %>%
  pivot_longer(cols = 4:5, 
               names_to = "Type",
               values_to = "Residue") %>%
  rename(protein_id = "Entry") %>%
  filter(str_detect(Residue, "\\d"),
         !is.na(Residue)) %>%
  mutate(position = as.integer(Residue)) %>%
  left_join(res.af, by = c("protein_id", "position")) %>%
  filter(!is.na(pPSE),
         !is.na(AA),
         AA == "C") %>%
  select(protein_id, quality, AA, position, pPSE, structure_group) %>%
  mutate(Class = "Disulfide") %>%
  distinct()



## Combine to one dataframe for plotting
uniprot.combined <- rbind(
  uniprot.nucleophile,
  uniprot.palmitoyl,
  uniprot.prenyl,
  uniprot.disulfide) %>%
  add_count(Class, name = "Total") %>%
  mutate(Facet = paste0(Class, "\nn = ", Total))



write_csv(uniprot.combined, "formatted_data/20221010_figure1_uniprot.csv")
