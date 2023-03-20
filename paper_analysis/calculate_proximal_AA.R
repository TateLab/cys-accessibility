## Load requisite packages
library(tidyverse)
library(magrittr)
library(data.table)
library(janitor)
library(scales)
library(broom)
library(fs)
library(here)


## Read in human fasta
res.fasta <- here("resources") %>%
  dir_ls(regexp = "FASTA") %>%
  fread() %>%
  mutate(protein_id = sub(Gene, 
                          pattern = "..\\|(.*)\\|.*",
                          replacement = "\\1")) %>%
  select(protein_id, Sequence)


## Human AlphaFold predictions
res.af <- here("resources") %>%
  dir_ls(regexp = "20230315_Supplementary_AlphaFold_StructureMapOutput.csv$") %>% 
  fread(drop="V1") %>%
  select(protein_id, AA, position, quality, colnames(.)[str_detect(colnames(.), "coord_c[ab]")])


## Function for getting the 3D distance of a residue from AlphaFold structure
calculate.proximal.aa.hq <- function(uniprot_id,
                                  residue_number,
                                  n_or_dist = c("n", "dist"),
                                  n_residue = 10,
                                  max_dist = 10){
  
  temp.af <- res.af %>% ### Looped over individual Cys on individual proteins
    filter(protein_id == uniprot_id) %>%
    arrange(position)
  
  temp.quality <- temp.af %>%
    select(protein_id, position, quality)
  
  if(nrow(temp.af) == 0){
    print("Protein not in existing database")
  }
  
  temp.aa.residue <- paste0(temp.af$AA, 
                            temp.af$position)
  
  ## Use position of beta-carbon from residue of interest (must be Cys, or at least not G)
  
  ## For G, use alpha-carbon
  temp.af$x_coord_cb[temp.af$AA == "G"] = temp.af$x_coord_ca[temp.af$AA == "G"]
  temp.af$y_coord_cb[temp.af$AA == "G"] = temp.af$y_coord_ca[temp.af$AA == "G"]
  temp.af$z_coord_cb[temp.af$AA == "G"] = temp.af$z_coord_ca[temp.af$AA == "G"]
  
  temp.dist.matrix <- temp.af %>%
    mutate(residue = paste0(AA, position)) %>%
    select(residue, x_coord_cb, y_coord_cb, z_coord_cb) %>%
    column_to_rownames(var = "residue") %>%
    dist() %>%
    as.matrix() %>%
    as.data.frame()
  
  colnames(temp.dist.matrix) <- temp.aa.residue
  
  rownames(temp.dist.matrix) <- temp.aa.residue
  
  ## Gives distance matrix for all residues in protein from AlphaFold coordinates

  temp.cys.with.aa.info <- temp.dist.matrix %>%
    select(residue_number) %>% ## Pick only Cys of interest
    rownames_to_column(var = "Residue") %>%
    rename(Distance = 2) %>%
    mutate(position = sub(Residue, pattern = "[A-Z](\\d*)", replacement = "\\1"),
           position = as.numeric(position),
           protein_id = uniprot_id) %>%
    left_join(temp.quality, by = c("protein_id", "position")) %>%
    filter(quality > 70) %>%
    slice_min(order_by = Distance, n = n_residue + 1)
  

  if(nrow(temp.cys.with.aa.info) == n_residue + 1){
    
    return(temp.cys.with.aa.info)
  
  }else{
    return(data.frame())
  }

}



## For all AF-predicted Cys with high-quality prediction (pLDDT > 70)
calculate.all.af.protein_id <- res.af %>%
  mutate(ID = paste0(protein_id, "_", position)) %>%
  filter(AA=="C",
         quality > 70) %>%
  select(protein_id, position, ID) %>%
  distinct() %>%
  pull(protein_id)

calculate.all.af.position <- res.af %>%
  mutate(ID = paste0(protein_id, "_", position)) %>%
  filter(AA=="C",
         quality > 70) %>%
  select(protein_id, position, ID) %>%
  distinct() %>%
  pull(position)

out.file <- paste0(here("resources"), "/20230313_HighQual_ProximalCys_178k.csv")

time.now <- Sys.time()

### For loop across all Cys sites
for(i in 1:length(calculate.all.af.position)){
  
  write_csv(calculate.proximal.aa.hq(uniprot_id = calculate.all.af.protein_id[i],
                                     residue_number = calculate.all.af.position[i],
                                     n_or_dist = "n",
                                     n_residue = 10),
            file = out.file,
            append = T)

  counter <- i / 500
  
  
  if(round(counter,0) == counter){
    
    message(i, " Cys in total, ", (Sys.time() - time.now), " for last 500")
    
    time.now <- Sys.time()
  }

}

all.cys.predicted <- out.file %>%
  fread() 


all.cys.export <- all.cys.predicted %>%
  rename(residue = 1, 
         distance = 2, 
         position = 3, 
         protein_id = 4,
         quality = 5) %>%
  distinct() %>%
  mutate(Central.Cys = ifelse(distance==0,
                              yes=1,
                              no=0),
         unique_cys_number = cumsum(Central.Cys)) %>%
  group_by(unique_cys_number) %>%
  mutate(cys_position = position[distance == 0]) %>% ## First row of each group is the central Cys
  ungroup() %>%
  add_count(protein_id, cys_position) %>%
  select(-Central.Cys, -n)

## Check correct number of rows/values


### Write table out to resources
write_tsv(all.cys.export, 
          "resources/20230315_Proximity_HighQualCys_NearestResidues_Long.tsv")


## All calculated with 10 nearby AA
proximal.master <- all.cys.export %>%
  group_by(cys_position, protein_id) %>%
  mutate(residue_rank = rank(distance, 
                             ties.method = "first",
                             na.last = NA),
         residue_rank = residue_rank - 1,
         AA = sub(residue, pattern = "([A-Z])\\d*", replacement = "\\1")) %>% 
  filter(max(distance) < 10, ## Remove Cys with few proximal Cys
         residue_rank %in% 1:5, ## Keep 5 closest for downstream analysis
         distance != 0) %>% ## Remove central Cys
  ungroup() %>%
  select(protein_id, cys_position, residue_rank, AA) %>%
  pivot_wider(names_from = "residue_rank",
              values_from = "AA") %>%
  unite(col = "Seq3D", colnames(.)[str_detect(colnames(.), "^\\d$")], sep = "", na.rm = T)


### Put a column with y/n for each AA presence within 5AA
AAs <- unique(res.af$AA)

for(i in 1:length(AAs)){
  
  proximal.master <- proximal.master %>%
    mutate(temp = ifelse(str_detect(Seq3D, AAs[i]),
                         yes=1,
                         no=0)) 
  
  colnames(proximal.master)[i+3] = AAs[i]
  
}

write_tsv(proximal.master,"resources/20230315_Proximity_HighQualCys_Sequences_Wide.tsv")


