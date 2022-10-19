
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
  separate(lineage, letters[1:7], sep = ";", remove = F)

d1$a %>% unique()
# [1] "Bacteria" 
  # >>> domain
d1$b %>% unique()
# [1] " Proteobacteria" " Firmicutes"     " Actinobacteria" 
  # >>> phylum
d1$c %>% unique()
# " Gammaproteobacteria" " Bacilli"             " Clostridia"         
# " Tissierellia"        " Negativicutes"       " Erysipelotrichia"   
# " Limnochordia"        " Micrococcales"   
# >>> All are class, except one.
# Micrococcales is an order in class Actinomycetia iin phylum Actinobacteria
# all others are class

# add class above Micrococcales
  # Remove row with missing class
  d1 <- d1 %>% 
    filter(!str_detect(c, "Micrococcales"))
  # correct and add back
    d1 <- d %>% 
      select(lineage) %>% 
      filter(str_detect(lineage, "Micrococcales")) %>% 
      mutate(lineage = str_replace(lineage, "Micrococcales", "Actinomycetia; Micrococcales")) %>% 
      separate(lineage, letters[1:7], sep = ";", remove = F) %>%
      bind_rows(d1,.)
  # check
    d1$c %>% unique()
  # looks good!

d1$d %>% unique()
# All end with -ales" except one 
  # >>> Order
# Sedimentibacter is a genus: Bacteria; Firmicutes; Tissierellia; Sedimentibacter
# This entry is missing both Order and Family
# according to Lawson
# "Sedimentibacter remains unaffiliated at the family level"
# Lawson, P.A. (2022). Tissierellaceae. In Bergey's Manual of Systematics of Archaea and Bacteria (eds M.E. Trujillo, S. Dedysh, P. DeVos, B. Hedlund, P. Kämpfer, F.A. Rainey and W.B. Whitman). https://doi.org/10.1002/9781118960608.fbm00275
# according to LPSN (https://lpsn.dsmz.de/family/eubacteriales-no-family)
# Order is Eubacteriales, Family is not assigned.
# class is not Tissierellia but Clostridia
# using this assignment

# add Order and Family, and change Class above Sedimentibacter
  # Remove row with missing class
    d1 <- d1 %>% 
      filter(!str_detect(d, "Sedimentibacter"))
  # correct and add back
    d1 <- d %>% 
      select(lineage) %>% 
      filter(str_detect(lineage, "Sedimentibacter")) %>% 
      mutate(lineage = str_replace(lineage, "Tissierellia.*", "Clostridia; Eubacteriales; Fam_NA; Sedimentibacter")) %>% 
      separate(lineage, letters[1:7], sep = ";", remove = F) %>%
      bind_rows(d1,.)
  # check
    d1$d %>% unique()
  # looks good!

d1$e %>% unique()
# >>> family
d1$f %>% unique()
# >>> genus

d <- d1 %>% 
  rename(domain = a,
         phylum = b,
         class = c,
         order = d,
         family = e,
         genus = f) %>% 
  select(-g) %>% 
  left_join(d,., by = "lineage")


# export data -------------------------------------------------------------
write_csv(d, here("genes_of_interest","spo0A","data","spo0A_KEGG_taxa.csv"))
# d <- read_csv(here("genes_of_interest","spo0A","data","spo0A_KEGG_taxa.csv"))
  