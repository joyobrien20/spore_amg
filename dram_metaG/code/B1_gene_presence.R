# Gene presence across samples
library(here)
library(tidyverse)
library(egg)
library(cowplot)
 
# Previously, we identified homologs of sporulation genes that were enriched in
# phages predicted to infect sporulator hosts.

# Here we will describe the distribution of those genes across samples.


# metagenomic data sets ---------------------------------------------------

data_dir <-"dram_metaG/data/dram_output"
sets <- list.dirs(here(data_dir), full.names = F, recursive = F)


# ecosystem classification
d.eco <- read_csv(here(data_dir, "../ecosystem_details.csv"))

# remove gvd and gpd
sets <- sets[ ! sets %in% c("Gregory_gvd","Camarillo-Guerrero")]
d.eco <- d.eco %>% 
  filter(! set %in% c("Gregory_gvd","Camarillo-Guerrero"))

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


# genomes and genes per dataset ------------------------------------------------
#In annotation file each row is an annotated gene, and scaffold is a genome 

d.genomes <- tibble()

for (cur_set in sets){

  annot_file <- 
    list.files(here(data_dir,cur_set), pattern = "annotation", full.names = T)
  
  # annot <- read_tsv(annot_file, col_select = "scaffold") %>% 
    annot <- read_tsv(annot_file) %>% 
      select("scaffold") %>% 
      pull(1)
  
  d.genomes <- 
    tibble(set = cur_set,
           n_genomes = length(unique(annot)),
           n_genes = length(annot)) %>% 
    bind_rows(d.genomes,.)
  
  
}
# clean up
rm(annot, annot_file)

# arrange data for plot ---------------------------------------------------

#summarize gene number per dataset
d.plot <-
  d.sum_amg %>% 
  filter(!is.na(gene_id)) %>% 
  filter(gene_id %in% enriched_ko) %>%
  pivot_longer(cols = !gene_id, names_to = "set", values_to = "n") 
  
  #normalize by total genes
d.plot <-
  left_join(d.plot, d.genomes, by = "set") %>% 
  mutate(frac_genes = n / n_genes)
  
  # add ecosystem type
  d.plot <-
    left_join(d.plot, d.eco, by = "set")
      

  # total observations per gene
  obs_amg <- 
    d.sum_amg %>% 
    filter(!is.na(gene_id)) %>% 
    filter(gene_id %in% spor_ko) %>%
    rowwise() %>% 
    mutate(obs = sum(across(where(is.numeric)))) %>% 
    arrange(desc(obs)) %>% 
    select(gene_id, obs)




  
  set_order <- 
    d.genomes %>% 
    arrange(n_genes) %>% 
    pull(set)
  


# plot --------------------------------------------------------------------

  

p <-
  d.plot %>% 
  # 0 to na
  filter(frac_genes>0) %>% 
  mutate(ecosystem_type = str_replace(ecosystem_type," ", "\n")) %>% 
    
  # order of Kos
  mutate(gene_id = fct_relevel(gene_id, obs_amg$gene_id)) %>%
  mutate(set = fct_relevel(set, set_order)) %>%
  
  #mark enriched
  mutate(enriched = if_else(gene_id %in% enriched_ko, TRUE, FALSE)) %>% 
  
  #plot
  ggplot(aes(x = gene_id, y =set))+
  geom_tile(aes(fill  = frac_genes), color = "white")+
  #mark enriched
  scale_fill_viridis_c(na.value = "white",
                       name='% of viral genes',
                       breaks = c(0.001,0.002,0.003),
                       labels = scales::percent(c(0.001,0.002,0.003), accuracy = 0.1))+
  scale_color_manual(values = c("transparent", "red"))+
  # separate by ecosystem
  facet_grid(ecosystem_type~., scales = "free_y", space = "free_y", switch = "y")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        axis.ticks = element_blank())

p1 <- obs_amg %>% 
  filter(gene_id %in% enriched_ko) %>% 
  mutate(gene_id = fct_inorder(gene_id)) %>%
  ggplot(aes(gene_id, obs))+
  geom_col()+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

p2 <- d.genomes %>% 
  left_join(., d.eco, by = "set") %>% 
  # mutate(ecosystem_type = str_replace(ecosystem_type," ", "\n")) %>% 
  mutate(set = fct_relevel(set, set_order)) %>%
  ggplot(aes(set, n_genes))+
  geom_col()+
  coord_flip()+
  facet_grid(ecosystem_type~., scales = "free_y", space = "free_y", switch = "y")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank())

ggarrange(p1, ggplot()+ theme_void(),
          p,#+theme(axis.text.x = element_blank())+ xlab("sporulation KOs"),
          p2+ ylab("N genes"), 
          nrow = 2, widths = c(5,1), heights = c(1,5)) %>% 
ggsave(here("dram_metaG/plots/","PA.png"), ., width = 8, height = 6)

