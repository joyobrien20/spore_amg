
# ------------------------------------------------------------------------
# Get AA sequences for a KEGG Ortholog 
# ------------------------------------------------------------------------


# load packaes ------------------------------------------------------------
library(here)
library(tidyverse)
# BiocManager::install("KEGGREST")
library(KEGGREST)
library(foreach)


# get list of genes ----------------------------------------------------------------
ko <- "ko:K07699" # spo0A

# get list of genes in KO
l.ko <- keggGet(ko)
# extract genes
genes <- l.ko[[1]]$GENES

# clean up gene names
g <-
  tibble(g = genes) %>% 
  # separate organism code from genes
  separate(g, ":| ", into = letters[1:10]) %>% 
  #some organisms have multiple genes
    pivot_longer(-1) %>% 
    filter(!is.na(value)) %>% 
    filter(value != "") %>% 
  #organim code to lower case
    mutate(a = tolower(a)) %>% 
  # some genes have unnecessary name in parenthesis. remove
    mutate(value = str_remove(value, "\\(.*?\\)")) %>%
  # pair back organim code with gene name
    mutate(g = paste(a,value,sep = ":")) %>% 
  # get clean result as vector
  pull(g) 



# Get data for each gene --------------------------------------------------

# initialize empty tibble
d <- tibble()

for(i in g){
  # get data from KEGG
  seqDB <- keggGet(i)
  # extract relavent data into table
  d <- 
    tibble(
      kegg_gene = i,  
      entry = seqDB[[1]]$ENTRY,
      organism = seqDB[[1]]$ORGANISM,
      DBlinks = paste(seqDB[[1]]$DBLINKS, collapse = ";"),
      seqAA = toString(seqDB[[1]]$AASEQ)) %>% 
    bind_rows(d,.)
  
  print(paste(which(g %in% i),"/",length(g)))
}


# export data -------------------------------------------------------------
write_csv(d, here("genes_of_interest","spo0A","data","spo0A_KEGG_sequences.csv"))
  