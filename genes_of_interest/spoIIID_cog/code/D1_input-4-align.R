library(here)
library(tidyverse)
cDir <- "genes_of_interest/spoIIID_cog"

# import sequences collected previously
data_dir <- here(cDir ,"data/")
d.bact <- read_csv(file = paste0(data_dir,"faa_COG_clustered.csv"))
d.vir <- read_csv(file = paste0(data_dir,"faa_virome_clustered.csv"))

d <- bind_rows(
  d.vir %>% 
    filter(protClstr_rep) %>% 
    select(SeqName, seqAA = sequence) %>% 
    mutate(SeqName = str_c("virome_",SeqName)),
  d.bact %>% 
    filter(protClstr_rep) %>% 
    select(SeqName  = protein_id, seqAA = seq) %>% 
    mutate(SeqName = str_c("cog_",SeqName))
)  

# > make multi fasta ------------------------

faa_file = paste0(data_dir, "seq2align.faa")

d %>%
  # prepare as fasta
  mutate(fasta = paste0(">",SeqName,"\n", seqAA)) %>%
  pull(fasta) %>%
  # write to file
  write_lines(., file = faa_file)
