library(here)
library(tidyverse)
library(cowplot)

data_dir <- ("dram_metaG/data/")

# map enriched genes back to taxa

# Levels of taxonomy


# enriched KOs -------------------------------------------------------------
d.enriched <- read_csv(here(data_dir, "enrichment","spor_enriched.csv"))

# GVD data ----------------------------------------------------------------

# AMGs detected by dram
amg.gvd <-  read_tsv(here(data_dir, "dram_output","Gregory_gvd", "amg_summary.tsv"))

# host sporulation from A2_gvd_taxa_spor.R
f_spor.gvd <- read_csv(here(data_dir, "enrichment", "gvd_spor_predictions.csv"))

# Arrange to bind with gpd
gtdb_tax_gvd <- c("gtdb_d", "gtdb_p","gtdb_f")

f_spor.gvd <- 
  f_spor.gvd %>%
  select(Contig, GTDB_Host, f_spor) %>% 
  separate(GTDB_Host, into = gtdb_tax_gvd, sep = ";", fill = "right") %>% 
  mutate(gtdb_p = str_remove(gtdb_p, "p__"),
         gtdb_f = str_remove(gtdb_f, "f__")) %>% 
  select(-gtdb_d)


# combine
amg.gvd <- 
  left_join(amg.gvd, f_spor.gvd, by = c("scaffold" = "Contig"))

# GPD data ----------------------------------------------------------------

# AMGs detected by dram
amg.gpd <-  read_tsv(here(data_dir, "dram_output","Camarillo-Guerrero", "amg_summary.tsv"))

# host sporulation from A2_gvd_taxa_spor.R
f_spor.gpd <- read_csv(here(data_dir, "enrichment", "gpd_spor_predictions.csv")) %>% 
  select(GPD_id, gtdb_tax_gvd[-1])

# combine
amg.gpd <- 
  left_join(amg.gpd, f_spor.gpd, by = c("scaffold" = "GPD_id"))        


# Combine data-sets --------------------------------------------------------

amg_spor <-
  bind_rows(
    amg.gvd %>% mutate(set = "gvd"),
    amg.gpd %>% mutate(set = "gpd")
  )



# Filter enriched KOs -----------------------------------------------------


amg_enriched <- 
  amg_spor %>% 
  filter(gene_id %in% d.enriched$gene_id)

# clean up
rm(amg.gpd, amg.gvd, amg_spor)
rm(f_spor.gpd, f_spor.gvd)


#adding non-sporulator to all non-Firmicutes
amg_enriched <- 
  amg_enriched %>% 
  mutate(f_spor = case_when(
    str_detect(gtdb_p,"Firmicutes") ~ f_spor,
    is.na(gtdb_p) ~ f_spor,
    TRUE ~ FALSE))



amg_phyla <-
  amg_enriched %>%
  group_by(gene_id,gtdb_p) %>% 
  summarise(n=n()) %>%
  arrange(desc(n)) %>% 
  mutate(firmi = case_when(
    str_detect(gtdb_p, "Firmicutes") ~ "Firmicutes",
    is.na(gtdb_p) ~ "unk. host",
    TRUE ~"Other phyla"))

# order of KO
ko_order <- amg_phyla %>% 
  group_by(gene_id,firmi) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = firmi, values_from = n, values_fill = 0) %>% 
  arrange(`Other phyla`) %>% 
  pull(gene_id)

p <- amg_phyla %>% 
  mutate(gene_id = factor(gene_id, levels = ko_order)) %>% 
  ggplot(aes(gene_id, gtdb_p))+
  geom_tile(aes(fill = log10(n)))+
  facet_grid(. ~ firmi, scales = "free", space = "free",
             labeller = labeller(firmi = label_wrap_gen(width = 8)))+
  scale_fill_viridis_c()+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text( face = "bold"),
        strip.background = element_blank())+
  ylab("Phylum")+
  coord_flip()
# p
ggsave2(here("dram_metaG/plots/amg_phyla.png"),p, width = 6, height = 8)
