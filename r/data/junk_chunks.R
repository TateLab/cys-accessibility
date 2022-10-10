
``` {R Junk functions}

find.residue <- function(protein, ## Uniprot ID
                         residue.nr){ ## Residue number from Uniprot annotation
  
  protein <- as.character(protein)
  residue.nr <- as.numeric(residue.nr)
  
  human.fasta[human.fasta$protein_id == protein,] %>%
    pull(Sequence) %>%
    substr(residue.nr, residue.nr)
  
}

#### JUNK CHUNKS

af.dubs <- dir_ls("/Users/mew21/Documents/GitHub/structuremap/nbs", regexp = "DUBs.csv") %>%
  fread() %>%
  rename(HalfShell = "nAA_12_70_pae",
         WholeShell = "nAA_24_180_pae") %>%
  left_join(gene.protein.id.match)

af.dubs %>%
  group_by(AA) %>%
  mutate(Total = n()) %>%
  ungroup() %>%
  group_by(AA, Total, HalfShell) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Part = n / Total) %>%
  ggplot(aes(x = HalfShell,
             y = Part)) + 
  geom_col() +
  facet_wrap(~AA)

```  ###########################################################


```{R 20220714 - Make list of ~10k downloaded AF structures}

cif.dir <- dir_ls("~/Documents/GitHub/structuremap/nbs", regexp = "cif")
pae.dir <- dir_ls("~/Documents/GitHub/structuremap/nbs", regexp = "pae")

cif.files <- dir_ls(cif.dir) %>%
  as.character() 

cif.ids <- cif.files %>%
  sub(pattern = "/Users/mew21/Documents/GitHub/structuremap/nbs/tutorial_cif/(.*).cif", 
      replacement = "\\1")

pae.files <- dir_ls(pae.dir) %>%
  as.character() 

pae.ids <- pae.files %>%
  sub(pattern = "/Users/mew21/Documents/GitHub/structuremap/nbs/tutorial_pae/pae_(.*).hdf", 
      replacement = "\\1")

## Overlap
overlap.ids <- intersect(pae.ids, 
                         cif.ids)

## Copy only overlapping files to new dir
overlap.cif.files <- paste0("/Users/mew21/Documents/GitHub/structuremap/nbs/tutorial_cif/", overlap.ids, ".cif")

length(overlap.cif.files)

overlap.pae.files <- paste0("/Users/mew21/Documents/GitHub/structuremap/nbs/tutorial_pae/pae_", overlap.ids, ".hdf")

length(overlap.pae.files)

head(overlap.cif.files)

## Make list of files *NOT* in both dirs
cif.delete <- cif.files[!cif.files %in% overlap.cif.files]

pae.delete <- pae.files[!pae.files %in% overlap.pae.files]

## delete extra cif files
file.remove(cif.delete)


## Make human fasta list with overlap
human.fasta.list <- fread("/Users/mew21/Desktop/HumanFASTA_Headers_AlphaFold.txt") 

human.fasta.overlap <- human.fasta.list %>%
  t() %>%
  as.data.frame() %>%
  filter(V1 %in% overlap.ids)

## Write to .txt
fileConn<-file("HumanFasta_AF_Headers_20220714-Overlap.txt")
writeLines(human.fasta.overlap$V1, fileConn)
close(fileConn)

## Read in human fasta
human.fasta <- dir_ls("/Users/mew21/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DUBs/Resources", regexp = "FASTA") %>%
  fread() %>%
  mutate(protein_id = sub(Gene, 
                          pattern = "..\\|(.*)\\|.*",
                          replacement = "\\1"))


########### Re-format human fasta headers to try downloading all in one go
hsap.headers <- human.fasta$protein_id %>%
  paste(collapse = "\t")

## Write to .txt
fileConn<-file("/Users/mew21/Downloads/hsapiens_headers.txt")
writeLines(hsap.headers, fileConn)
close(fileConn)

```  ###########################################################


``` {R Amino acid enrichments in each structural featuer}

rbind(
  af.pred %>%
    filter(secondary_structure == "TURN_TY1_P") %>%
    mutate(Total = nrow(.)) %>%
    group_by(AA, Total) %>%
    summarise(Part=n(), .groups = "drop") %>%
    mutate(Part = Part / Total,
           Sample = "Turns")
  ,
  af.pred %>%
    #filter(secondary_structure == "TURN_TY1_P") %>%
    mutate(Total = nrow(.)) %>%
    group_by(AA, Total) %>%
    summarise(Part=n(), .groups = "drop") %>%
    mutate(Part = Part / Total,
           Sample = "All")) %>%
  select(-Total) %>%
  pivot_wider(names_from = Sample,
              values_from = Part) %>%
  mutate(`Turns.Enrichment` = Turns - All) %>% View

```
###################################################################

## Histograms of adjacent AAs (1 / accessibility)
af.pred %>%
  ggplot(aes(x = WholeShell)) + 
  geom_histogram(binwidth = 1) +
  facet_wrap(~AA) + 
  xlim(0,50)

af.pred %>%
  filter(quality > 80) %>%
  ggplot(aes(x = WholeShell,
             y = HalfShell)) + 
  geom_point(alpha = 0.05) + 
  facet_wrap(~AA)

###################################################################


``` {R Compare cysteine and whole proteome distrbutions}
##


## Cys
af.cys <- af.pred %>%
  filter(AA == "C",
         #quality > 80
  ) %>%
  select(protein_id, AA, position, HalfShell) %>%
  mutate(Total = length(AA),
         Class = "Cysteine")

## Ser
af.ser <- af.pred %>%
  filter(AA == "S",
         #quality > 80
  ) %>%
  select(protein_id, AA, position, HalfShell) %>%
  mutate(Total = length(AA),
         Class = "Serine")

## Phe
af.phe <- af.pred %>%
  filter(AA == "F",
         #quality > 80
  ) %>%
  select(protein_id, AA, position, HalfShell) %>%
  mutate(Total = length(AA),
         Class = "Phenylalanine")


#### All amino acids vs cysteine

rbind(af.ser,
      af.phe,
      af.cys) %>%
  group_by(Class, Total, HalfShell) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Part = n / Total) %>%
  left_join(af.all.aa.ref) %>%
  mutate(Part.Norm = Part - Ref.Part) %>%
  ggplot(aes(x = HalfShell,
             y = 100*Part.Norm,
             fill = Class)) + 
  geom_col(position = position_identity()) +
  facet_wrap(~factor(Class,
                     levels = c(#"All amino acids", 
                       "Serine", "Phenylalanine", "Cysteine")), 
             ncol = 2) + 
  xlab("Amino Acid Depth") + 
  xlim(-0.5,16) +
  theme_linedraw() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        legend.position = "none") + 
  ylab("Proportion relative to whole proteome") + 
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)))

#ggsave("AF/Plots/CysvsAllAA_Distributions.png", height = 4, width = 5)


```


###################################################################


``` {R Distributions of metal-binding cysteines by metal}

uniprot.metalcys %>% 
  #filter(quality > 80) %>%
  group_by(Metal) %>%
  mutate(Total = n()) %>%
  ungroup() %>%
  group_by(Class, Total, HalfShell, Metal) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Part = n / Total) %>%
  left_join(human.cys.wp, by = "HalfShell") %>% 
  mutate(Part.Norm = Part - Ref.Part) %>%
  ggplot(aes(x = HalfShell,
             y = 100*Part.Norm,
             fill = Class)) + 
  geom_col(position = position_identity()) +
  facet_wrap(~paste0(Metal, " (n=", Total, ")")) + 
  xlab("Amino Acid Depth") + 
  #xlim(-0.5,16) +
  theme_linedraw() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        legend.position = "none") + 
  ylab("% enrichment vs. whole proteome") + 
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)))


```



###################################################################


vino.scout.cys %>%
  filter(quality > 80) %>%
  group_by(R.bin, structure_group) %>%
  mutate(Total = n()) %>%
  ungroup() %>%
  group_by(R.bin, Total, HalfShell, structure_group) %>%
  summarise(n = n(),
            .groups = "drop") %>% 
  mutate(Part = 100*n/Total) %>% 
  select(structure_group, HalfShell, R.bin, Part) %>%
  pivot_wider(names_from = R.bin,
              values_from = Part) %>% 
  mutate(Liganded = ifelse(is.na(Liganded),
                           yes=0,
                           no = Liganded),
         Difference = Liganded-`Not liganded`) %>%
  ggplot(aes(x = HalfShell,
             y = Difference)) +
  geom_col() +
  facet_wrap(~structure_group) + 
  theme_linedraw() + 
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(colour = "black")) + 
  ylab("% of residues") + 
  xlab("Amino Acid Depth")

###################################################################


##### as heatmap
vino.scout.cys %>%
  filter(quality > 80) %>%
  group_by(R.bin, structure_group) %>%
  mutate(Total = n()) %>%
  ungroup() %>%
  group_by(R.bin, Total, HalfShell, structure_group) %>%
  summarise(n = n(),
            .groups = "drop") %>% 
  mutate(Part = 100*n/Total) %>% 
  select(structure_group, HalfShell, R.bin, Part) %>%
  pivot_wider(names_from = R.bin,
              values_from = Part) %>% 
  mutate(Liganded = ifelse(is.na(Liganded),
                           yes=0,
                           no = Liganded),
         Difference = Liganded-`Not liganded`) %>%
  ggplot(aes(x = HalfShell,
             y = structure_group,
             fill = Difference)) + 
  geom_tile() + 
  scale_fill_gradient2() + 
  xlab("Amino acid depth") + 
  ylab("Structural domain") + 
  theme_linedraw() +
  theme(panel.background = element_rect(fill = "light grey"),
        panel.grid = element_blank()) + 
  scale_x_continuous(expand = expansion()) +
  scale_y_discrete(expand = expansion())

###################################################################


###################################################################


###################################################################


###################################################################


###################################################################


###################################################################


###################################################################


###################################################################


###################################################################


###################################################################


###################################################################


###################################################################



