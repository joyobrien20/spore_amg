library(here)
library(tidyverse)
library(cowplot)

data_dir <- ("dram_metaG/data/")

# Mapping the list of genes from the enrichment analysis back to locus tags

d.enriched <- read_csv(here(data_dir, "enrichment","gvd_gpd_enrich.csv"))

p.signif <- 1e-6

k_treshold <- 32

d.enriched %>% 
  mutate(enriched = if_else(significant, "enriched", "not enriched")) %>% 
  group_by(spor_gene, enriched) %>% 
  summarise(n=n()) %>% 
  pivot_wider(names_from = enriched, values_from = n)

# spor_gene        enriched `not enriched`
# <chr>               <int>          <int>
# 1 other                   5            795
# 2 sporulation_gene       54            332


# map back to genes -------------------------------------------------------

# Get gene lists

bs.genes <- read_csv(here("spor_gene_list/data/dram_spore_genes_RS.csv")) %>% 
  separate(gene_description, into = c("symbol", "description"), sep = ";") %>% 
  rename(ko = gene_id.ko)

cd.genes <- read_csv(here("spor_gene_list/data/cdif_spor_KOs.csv"))

# enriched spor KOs
spor.enriched.ko <- 
  d.enriched %>% 
  filter(spor_gene == "sporulation_gene") %>%
  filter(significant) %>% pull(gene_id) %>% unique()

genes_enriched <- 
  bs.genes %>% 
  filter(ko %in% spor.enriched.ko) %>% 
  select(ko , symbol, description ) %>% 
  mutate(sp = "bs")

genes_enriched <-
  cd.genes %>% 
  filter(ko %in% spor.enriched.ko) %>% 
  select(ko , symbol, description ) %>% 
  mutate(sp = "cd") %>% 
  bind_rows(genes_enriched, .)

# shared.enriched <- 
#   genes_enriched %>% 
#   group_by(ko,sp) %>% 
#   summarise(n=n()) %>% 
#   select(ko, sp) %>% 
#   mutate(presence = 1) %>% 
#   pivot_wider(names_from = sp, values_from = presence, values_fill = 0) %>% 
#   mutate(shared = (bs+cd)>1) %>% 
#   filter(shared) %>% 
#   pull(ko) %>% 
#   unique()
# 
# genes_enriched %>% 
#   filter(ko %in% shared.enriched) 


# add to amg data frame
d <- genes_enriched %>% 
  # select(-description) %>% 
  mutate(description = str_c(symbol," [",description,"]")) %>% 
  group_by(ko, sp) %>% 
  summarise(genes = str_c(description, collapse = ";")) %>% 
  pivot_wider(names_from = sp, values_from = genes ) %>% 
  left_join(d.enriched %>% filter(gene_id %in% spor.enriched.ko), .,
            by = c("gene_id" = "ko")) %>% 
  filter(spor_gene == "sporulation_gene") 

write_csv(d, here(data_dir, "enrichment","spor_enriched.csv"))
