# Enrichment analysis for AMG in phages of sporulating hosts

library(here)
library(tidyverse)
library(cowplot)
library(ggrepel)

data_dir <- ("dram_metaG/data/")

# GVD data ----------------------------------------------------------------

# AMGs detected by dram
amg.gvd <-  read_tsv(here(data_dir, "dram_output","Gregory_gvd", "amg_summary.tsv"))

# host sporulation from A2_gvd_taxa_spor.R
f_spor.gvd <- 
  read_csv(here(data_dir, "enrichment", "gvd_spor_predictions.csv")) %>% 
  select(Contig, f_spor)
  

# combine
amg.gvd <- 
  left_join(amg.gvd, f_spor.gvd, by = c("scaffold" = "Contig"))
               
# GPD data ----------------------------------------------------------------

# AMGs detected by dram
amg.gpd <-  read_tsv(here(data_dir, "dram_output","Camarillo-Guerrero", "amg_summary.tsv"))

# host sporulation from A2_gvd_taxa_spor.R
f_spor.gpd <- 
  read_csv(here(data_dir, "enrichment", "gpd_spor_predictions.csv")) %>% 
  select(GPD_id, f_spor)

# combine
amg.gpd <- 
  left_join(amg.gpd, f_spor.gpd, by = c("scaffold" = "GPD_id"))        


# Combine data-sets --------------------------------------------------------

amg_spor <-
  bind_rows(
    amg.gvd %>% mutate(set = "gvd"),
    amg.gpd %>% mutate(set = "gpd")
  )



# > Add data on sporulation genes -------------------------------------------

#list of sporulation genes
spor_genes <- read_csv(here("spor_gene_list", "data", "dram_spore_genes.csv"))

amg_spor <- 
  amg_spor %>% 
  mutate(spor_gene = if_else(gene_id %in% spor_genes$gene_id.ko,
                           "sporulation_gene", "other"))

# Summarise total number of sporulators and non sporulators ---------------

sum_spor <- 
  bind_rows(f_spor.gpd %>% rename(scaffold = GPD_id) %>% mutate(set = "gpd"),
            f_spor.gvd%>% mutate(set = "gvd"))%>% 
  group_by(set, f_spor) %>% 
  summarise(n=n())



# Enrichment test combined data-sets ---------------------------------------


#significance threshold
p.signif <- 1e-6

k_treshold <- 32

total_sporulators <- 
  sum_spor %>% 
  filter(f_spor) %>% 
  pull(n) %>% sum()

total_nonSporulators <- 
  sum_spor %>% 
  filter(!f_spor) %>% 
  pull(n) %>% sum()

d_enrich <- 
  amg_spor %>% 
  # Summarize KOs by host sporulation
  mutate(f_spor = case_when(
    is.na(f_spor) ~ "unkownSpor_host",
    f_spor ~ "spor_host",
    ! f_spor ~ "nonSpor_host"
  )) %>% 
  filter(!is.na(gene_id)) %>%  
  group_by(gene_id, spor_gene, f_spor) %>% 
  summarise(n=n(), .groups = "drop") %>% 
  pivot_wider(values_from = n, names_from = f_spor, values_fill = 0) %>% 
    
    # q	: the number of white balls drawn without replacement 
    # from an urn which contains both black and white balls.
    # >>> AMG (specific gene) in viruses infecting sporulators
    mutate(q = spor_host) %>% 
    
    # m	:the number of white balls in the urn.
    # >>> number of viruses infecting sporulators in the pool
    mutate(m = total_sporulators) %>% 
    
    # n	:the number of black balls in the urn.
    # >>> number of viruses infecting non-sporulators in the pool
    mutate(n = total_nonSporulators) %>% 
    
    # k: the number of balls drawn from the urn, hence must be in 0,1,., m+n.
    # >>> Total AMG (specific gene) detected in viruses
    mutate(k = spor_host + nonSpor_host) %>% 
    
    mutate(p.val = signif (phyper(q, m, n, k, lower.tail =F),5)) %>% 
  # Adjust Pvalue for multiple testing
    mutate(p.adj = signif (p.adjust(p.val, method = "BH"),5)) %>% 
  # mark significantly enriched genes
   mutate(significant = (p.adj < p.signif) & (k>k_treshold))


# export analysis test results
write_csv(d_enrich %>% arrange(desc(spor_gene), desc(significant)),
          here(data_dir, "enrichment","gvd_gpd_enrich.csv"))


# plot --------------------------------------------------------------------

p <-
  d_enrich %>% 
  ggplot(aes(k,-log10(p.adj+1e-24)))+
  geom_hline(yintercept = -log10(p.signif), linetype=2, color="grey", size = 1)+
  geom_vline(xintercept = k_treshold, linetype=2, color="grey", size = 1)+
  
  # geom_text_repel(daaes(label = lab), max.overlaps=50,color = "grey", size = 3)+
  
  geom_jitter(aes(fill = significant),width = 0.005, height = 0.005,
              shape=21, size=3, stroke = 1, alpha = 0.5, show.legend = F)+
  theme_classic()+
  facet_wrap(~ spor_gene %>% fct_rev(), ncol = 1) +
  panel_border(color = "black")+
  scale_x_continuous(trans = "log2", breaks = (2^(0:11)))+
  xlab("Sample size\nNo. homologs detected (log2)")+
  ylab("Enrichment (-log10 adj. P-value)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(p, filename = here("dram_metaG/plots", "enrich_gvd_gpd.png"), 
       width = 4, height = 6)  

