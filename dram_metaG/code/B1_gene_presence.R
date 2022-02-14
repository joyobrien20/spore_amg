# Gene presence across samples
library(here)
library(tidyverse)
library(egg)
 
# Previously, we identified homologs of sporulation genes that were enriched in
# phages predicted to infect sporulator hosts.

# Here we will describe the distribution of those genes across samples.


# metagenomic data sets ---------------------------------------------------

data_dir <-"dram_metaG/data/dram_output"
sets <- list.dirs(here(data_dir), full.names = F, recursive = F)

# remove extra data
sets <- sets[sets!="extra_data"]


# Sporulation genes -------------------------------------------------------

spor_genes <- read_csv(here("spor_gene_list/data/dram_spore_genes.csv")) %>% 
  # genes with KO
  filter(!is.na(gene_id.ko)) %>% 
  separate(gene_description, into = c("symbol", "description"), sep = ";") %>% 
  rename(ko = gene_id.ko)

spor_ko <- 
  unique(spor_genes$ko)

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

# total amg sum per set
  total_amg <- 
    d.sum_amg %>% 
    summarise(across(where(is.numeric), sum)) %>% 
    pivot_longer(everything(),names_to = "set", values_to = "n_total")
  
  # total observations per gene
  obs_amg <- 
    d.sum_amg %>% 
    filter(!is.na(gene_id)) %>% 
    filter(gene_id %in% spor_ko) %>%
    rowwise() %>% 
    mutate(obs = sum(across(where(is.numeric)))) %>% 
    arrange(desc(obs)) %>% 
    select(gene_id, obs)

  
  set_order <- total_amg %>% 
    arrange(n_total) %>% 
    pull(set)

p <-   d.sum_amg %>% 
  filter(!is.na(gene_id)) %>% 
    filter(gene_id %in% spor_ko) %>%
  pivot_longer(cols = !gene_id, names_to = "set", values_to = "n") %>% 

  #normalize by totals
  left_join(., total_amg) %>% 
  mutate(perc = 100*n/n_total) %>% 
    
  # order of Kos
  mutate(gene_id = fct_relevel(gene_id, obs_amg$gene_id)) %>%
  mutate(set = fct_relevel(set, set_order)) %>%
  
  #mark enriched
  mutate(enriched = if_else(gene_id %in% enriched_ko, TRUE, FALSE)) %>% 
  
  #plot
  ggplot(aes(x = gene_id, y =set))+
  geom_tile(aes(fill  = log10(perc)), color = "white")+
  # scale_fill_gradient(low = "pink", high = "red", na.value = "grey")+
  #mark enriched
  geom_point(aes(color = enriched), y=Inf, shape = 25, size=1, show.legend = F)+
  scale_fill_viridis_b(na.value = "white")+
  scale_color_manual(values = c("transparent", "red"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        axis.ticks = element_blank())

p1 <- obs_amg %>% 
  mutate(gene_id = fct_inorder(gene_id)) %>%
  ggplot(aes(gene_id, obs))+
  geom_col()+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank())

p2 <- total_amg %>% 
  mutate(set = fct_relevel(set, set_order)) %>%
  ggplot(aes(set, n_total))+
  geom_col()+
  coord_flip()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

ggarrange(p1, ggplot()+ theme_void(),
          p+theme(axis.text.x = element_blank())+ xlab("sporulation KOs"),
          p2+ ylab("N amg"), 
          nrow = 2, widths = c(5,1), heights = c(1,5)) %>% 
ggsave(here("dram_metaG/plots/","PA.png"), ., width = 8, height = 4)

