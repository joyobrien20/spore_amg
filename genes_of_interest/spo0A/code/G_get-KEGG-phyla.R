
# ------------------------------------------------------------------------
# Get taxa data for  KEGG Orthologs 
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
  #organim code to lower case
    mutate(a = tolower(a)) %>% 
  # get clean result as vector
  pull(a) 
g <- unique(g)
keggGet("genome:mpc")

# Get data for each gene --------------------------------------------------

# initialize empty tibble
d <- tibble()

for(i in g){
  # get data from KEGG
  seqDB <- keggGet(paste0("genome:",i))
  # extract relavent data into table
  d <- 
    tibble(
      organism = i,  
      keggID = seqDB[[1]]$ENTRY,
      taxID = seqDB[[1]]$TAXONOMY$TAXONOMY,
      lineage = seqDB[[1]]$TAXONOMY$LINEAGE,
      seqID = seqDB[[1]]$CHROMOSOME$SEQUENCE
      ) %>% 
    bind_rows(d,.)
  
  print(paste(which(g %in% i),"/",length(g)))
}



# parse lineage -----------------------------------------------------------
# str_count(d$lineage,";") %>% range() # 3 6
# # 4-7 taxa levels

d1 <- d %>% 
  select(lineage) %>% 
  separate(lineage, letters[1:7], sep = ";")

d1$a %>% unique()
# [1] "Bacteria" >>> domain
d1$b %>% unique()
# [1] " Proteobacteria" " Firmicutes"     " Actinobacteria" >>> phylum
d1$c %>% unique()
# >>> class/order
# Micrococcales is an order in class Actinomycetia iin phylum Actinobacteria
# all others are class
d1$d %>% unique()
# >>> class/other
# Micrococcaceae is a family, in order Micrococcales
# Sedimentibacter is a genus: Bacteria; Firmicutes; Tissierellia; Sedimentibacter
# 
d1$e %>% unique()
# >>> family/others
d1$f %>% unique()
# >>> genus


# export data -------------------------------------------------------------
write_csv(d, here("genes_of_interest","spo0A","data","spo0A_KEGG_taxa.csv"))
# d <- read_csv(here("genes_of_interest","spo0A","data","spo0A_KEGG_taxa.csv"))
  