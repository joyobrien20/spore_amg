# Gene presence across samples
library(here)
library(tidyverse)
 
# Previously, we identified homologs of sporulation genes that were enriched in
# phages predicted to infect sporulator hosts.

# Here we will describe the distribution of those genes across samples.


# metagenomic data sets ---------------------------------------------------

data_dir <-"dram_metaG/data/dram_output"
sets <- list.dirs(here(data_dir), full.names = F, recursive = F)

# remove extra data
sets <- sets[sets!="extra_data"]


# enriched genes ----------------------------------------------------------

d.enriched <- read_csv(here("metaG/data/gvd_enrich.csv")) %>% 
  filter(K > 8) %>% 
  filter(adj.p < 0.001) %>% 
  filter(spor_gene == "sporulation_gene")

#enriched KOs by P-value
enriched_ko <- d.enriched %>%
  arrange(adj.p) %>% pull(gene_id)

# count amgs --------------------------------------------------------------
d.sum_amg <- tibble(gene_id = NA)

for (set in sets){
  #load amgs
  d.amg <- read_tsv(here(data_dir, set, "amg_summary.tsv"))
  
  d.sum_amg <-
    d.amg %>% 
    group_by(gene_id) %>% 
    summarise(n=n()) %>% 
    rename({{set}} := n) %>% #assign column name from loop var
        # # transpose
        # pivot_wider(names_from = gene_id, values_from = n) %>% 
        # mutate(set = set) %>% 
        # relocate(set)%>% 
        # bind_rows(d.sum_amg,.) 
    full_join(d.sum_amg, .)
  
}

# replace NAs with 0
d.sum_amg <- d.sum_amg %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)))

# total amg sum
  total_amg <- 
    d.sum_amg %>% 
    summarise(across(where(is.numeric), sum)) %>% 
    pivot_longer(everything(),names_to = "set", values_to = "n_total")

  d.sum_amg %>% 
    filter(gene_id %in% enriched_ko) %>% 
  pivot_longer(cols = !gene_id, names_to = "set", values_to = "n") %>% 

  #normalize by totals
  left_join(., total_amg) %>% 
  mutate(perc = 100*n/n_total) %>% 
    
  # order of Kos
    mutate(gene_id = fct_relevel(gene_id, enriched_ko)) %>% 
  ggplot(aes(gene_id, set))+
  geom_tile(aes(fill  = perc), color = "white")+
  scale_fill_gradient(low = "white", high = "black")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
