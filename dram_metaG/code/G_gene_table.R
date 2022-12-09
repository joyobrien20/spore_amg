library(here)
library(tidyverse)
library(cowplot)
library(KEGGREST)


# get data from reviewed AMGs ----------------------------------------------

data_dir <- ("dram_metaG/data/")
d <- read_csv(here(data_dir,"viral_sporGenes.csv"))


# get sporulation gene list -----------------------------------------------

d_spor <- read_csv(here("spor_gene_list","data", "dram_spore_genes.csv"))

# only enriched-validated
d_spor <- d_spor %>%
  filter(gene_id.ko %in% d$gene_id)

cd_spor <- read_csv(here("spor_gene_list","data", "cdif_spor_KOs.csv"))
#B. subtilis
bs_spor <- read_csv(here("spor_gene_list","data", "dram_spore_genes_RS.csv"))
  # to add locus tags
  bs_tags <- read_csv(here("spor_gene_list","data", "batch.csv"))
  bs_tags <- select(bs_tags, locus_tag = locus, symbol = name)

d.spor_genes <- bind_rows(

  cd_spor %>%
  filter(ko %in% d$gene_id) %>% 
  select(ko, locus_tag, symbol, description) %>% 
  mutate (strain = "Cd"),

  bs_spor %>% 
  filter(gene_id.ko %in% d$gene_id) %>% 
  separate(gene_description, into = c("symbol", "description"), sep = ";") %>% 
  select(ko = gene_id.ko, symbol, description) %>% 
  left_join(.,bs_tags) %>% 
  mutate (strain = "Bs")
  
)

# gett al genes that correspond to each KO
# Cd KEGG ----------------------------------------------------------------

# all kegg genes
raw.kegg.cdf <- keggFind("genes", "cdf:CD630") 
d.kegg.cdf <- raw.kegg.cdf %>% 
  enframe(name = "kegg", value = "kegg.txt") %>% 
  separate(kegg, into = c("strain", "locus_tag"), sep = ":") %>% 
  separate(kegg.txt, into = c("symbol", "description"), sep = ";", fill = "left", extra = "merge") 


# KOs
ko.cd <- keggLink("cdf", "ko")
ko.cd <- enframe(ko.cd, name = "ko", value = "cdf")
ko.cd$cdf <-  str_replace(ko.cd$cdf,pattern = "cdf:",replacement = "")
ko.cd$ko <-  str_replace(ko.cd$ko,pattern = "ko:",replacement = "")

d.kegg.cdf <- left_join(d.kegg.cdf, ko.cd, by = c("locus_tag"="cdf"))


# B. subtilis KEGG genes --------------------------------------------------


raw.kegg.bsu <- keggFind("genes", "bsu:BSU") 
d.kegg.bsu <- raw.kegg.bsu %>% 
  enframe(name = "kegg", value = "kegg.txt") %>% 
  separate(kegg, into = c("strain", "locus_tag"), sep = ":") %>% 
  separate(kegg.txt, into = c("symbol", "description"), sep = ";", fill = "left", extra = "merge") 

# KOs
ko.bs <- keggLink("bsu", "ko")
ko.bs <- enframe(ko.bs, name = "ko", value = "bsu")
ko.bs$bsu <-  str_replace(ko.bs$bsu,pattern = "bsu:",replacement = "")
ko.bs$ko <-  str_replace(ko.bs$ko,pattern = "ko:",replacement = "")

d.kegg.bsu <- left_join(d.kegg.bsu, ko.bs, by = c("locus_tag"="bsu"))


# all genes of viral spore KOs --------------------------------------------


d.genes <- bind_rows(
  d.kegg.bsu %>%  filter(ko %in% d$gene_id),
  d.kegg.cdf %>%  filter(ko %in% d$gene_id)) %>% 
  relocate(ko, .before = 1) %>% 
  arrange(ko)

# mark sporulation genes
d.genes <- d.genes %>% 
  #adjust to match
  mutate(locus_tag = str_replace(locus_tag, "BSU", "BSU_")) %>% 
  mutate(spor.gene = if_else(
    locus_tag %in% d.spor_genes$locus_tag,TRUE,FALSE))

# write tables

  #sporulation genes
d.genes %>% 
  filter(spor.gene) %>% 
  group_by(ko, spor.gene, strain ) %>%
  summarise(loci = str_c(locus_tag, collapse = "\n"),
            symbols = str_c(symbol, collapse = "\n"),
            descriptions = str_c(description, collapse = "\n"),
            spor.gene = str_c(spor.gene, collapse = "\n")) %>%
  # make nice
  mutate(strain = str_replace(strain, "bsu", "BSU")) %>% 
  mutate(strain = str_replace(strain, "cdf", "CD630")) %>% 
  mutate(loci = str_remove_all(loci, "BSU_|CD630_")) %>% 
  mutate(spor.gene = spor.gene %>%  as.logical()) %>% 
  arrange(ko, strain, loci) %>% 
  ungroup() %>% 
  select(-spor.gene) %>%
  write_csv(here(data_dir,"genes_spor_out.csv"))

  #NONsporulation genes
d.genes %>% 
  filter(!spor.gene) %>% 
  group_by(ko, spor.gene, strain ) %>%
  summarise(loci = str_c(locus_tag, collapse = "\n"),
            symbols = str_c(symbol, collapse = "\n"),
            descriptions = str_c(description, collapse = "\n"),
            spor.gene = str_c(spor.gene, collapse = "\n")) %>%
  # make nice
  mutate(strain = str_replace(strain, "bsu", "BSU")) %>% 
  mutate(strain = str_replace(strain, "cdf", "CD630")) %>% 
  mutate(loci = str_remove_all(loci, "BSU_|CD630_")) %>% 
  mutate(spor.gene = spor.gene %>%  as.logical()) %>% 
  arrange(ko, strain, loci) %>% 
  ungroup() %>% 
  select(-spor.gene) %>%
  write_csv(here(data_dir,"genes_nonSPOR_out.csv"))
