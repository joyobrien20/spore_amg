library(here)
library(tidyverse)
library(cowplot)

data_dir <- ("dram_metaG/data/")

# get data from reviewed AMGs ----------------------------------------------

# manual curation results are in "curated" folder

dir.curated <- here(data_dir, "curated")
done <- list.files(dir.curated, pattern = ".csv", full.names = T)

d.curated <- tibble()
for (f in done){
  d.curated <-
    read_csv(f, col_types = cols(
    focal_gene = col_character(),
    set_group = col_character(),
    set = col_character(),
    idx = col_double(),
    scaffold = col_character(),
    v.cat = col_character(),
    skim_plot_isViral = col_character()
  )) %>% 
    bind_rows(d.curated,.)
}


# Focus on Gut viromes ----------------------------------------------------

# summarize scaffolds
d.curated.sum <- d.curated %>% 
  # use gut viromes
  filter(set_group != "CSUsets") %>% 
  group_by(focal_gene, skim_plot_isViral) %>% 
  summarise(n_scaffolds = n()) %>% 
  # to add 0 where needed
  pivot_wider(names_from = skim_plot_isViral, 
              values_from = n_scaffolds, 
              values_fill = 0) %>% 
  pivot_longer(-1, names_to = "skim_plot_isViral", values_to = "n_scaffolds" ) %>% 
  mutate(gene = str_remove(focal_gene, ".*_"))

# make calls based on number of TRUE observations
d.curated.sum <- d.curated.sum %>% 
  filter(skim_plot_isViral=="TRUE") %>% 
  mutate(AMG = case_when(
    n_scaffolds > 5 ~ "Likely viral",
    n_scaffolds < 1 ~ "Unlikely viral",
    TRUE ~ "Possible viral"
  )) %>% 
  # add calls
  select(focal_gene, AMG) %>% 
  left_join(d.curated.sum, by = "focal_gene") 

# plot scaffolds detected ------------------------------------------------------
# arrange plotting order
genes <- d.curated.sum %>% 
  group_by(gene) %>% 
  summarise(total = sum(n_scaffolds)) %>% 
  arrange(total %>% desc()) %>% 
  pull(gene)

#plotting order
d.curated.sum <- d.curated.sum %>% 
  mutate(gene = factor(gene, levels = genes))
# summarize data
d.scafs <- 
  d.curated.sum %>% 
  pivot_wider(names_from = skim_plot_isViral, values_from = n_scaffolds) %>% 
  mutate(inspected = `FALSE`+`TRUE`+MAYBE,
         detected = `FALSE`+`TRUE`+MAYBE+ `NA`)
  

# plot summary
p1 <- d.scafs %>% 
  filter(AMG == "Likely viral") %>%
  ggplot(aes(gene))+
  geom_col(aes(y= detected, fill = "detected"))+
  geom_col(aes(y=inspected, fill = "inspected"))+
  
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "none")+
  scale_fill_manual(values = c("grey20", "red"),name=NULL)+
  scale_y_log10(breaks = 10^(0:4))+
  annotation_logticks(sides = "l")+
  ylab("N scaffolds")

# p1

p2 <- d.scafs %>% 
  filter(AMG != "Likely viral") %>%
  ggplot(aes(gene))+
  geom_col(aes(y= detected, fill = "detected"))+
  geom_col(aes(y=inspected, fill = "inspected"))+
  
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "bottom")+
  scale_fill_manual(values = c("grey20", "red"),name=NULL)+
  scale_y_log10(breaks = 10^(0:4))+
  annotation_logticks(sides = "l")+
  ylab("N scaffolds")

# p2

p_detected <- plot_grid(p1,p2, ncol=1, rel_heights = c(.8,1))


# plot scaffolds inspected ------------------------------------------------------


# plot 
p1 <- d.curated.sum %>% 
  filter(!(skim_plot_isViral == "NA")) %>%
  rename(is_AMG_Viral = skim_plot_isViral) %>% 
  group_by(gene) %>% 
  mutate(perc_scaffolds = 100*n_scaffolds/sum(n_scaffolds)) %>% 
  filter(AMG == "Likely viral") %>%
  ggplot(aes(gene, perc_scaffolds))+
  geom_col(aes(fill = is_AMG_Viral))+
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "none")+
  scale_fill_viridis_d(direction = -1)+
  ylab("% of inspected")

p2 <- d.curated.sum %>% 
  filter(!(skim_plot_isViral == "NA")) %>%
  rename(is_AMG_Viral = skim_plot_isViral) %>% 
  group_by(gene) %>% 
  mutate(perc_scaffolds = 100*n_scaffolds/sum(n_scaffolds)) %>% 
  filter(!(AMG == "Likely viral")) %>%
  ggplot(aes(gene, perc_scaffolds))+
  geom_col(aes(fill = is_AMG_Viral))+
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "bottom")+
  scale_fill_viridis_d(direction = -1)+
  ylab("% of inspected")

p_inspected <- plot_grid(p1,p2, ncol=1, rel_heights = c(.8,1))

p <- plot_grid(NULL, p_detected, NULL ,p_inspected, ncol = 1, 
               rel_heights = c(0.1,1,0.1,1), labels = c("(a)","","(b)"))
ggsave(filename = here("dram_metaG/plots","GUT_scaffold_counts.png"),
       plot = p, width = 8, height = 9, bg = "white")



# summarize counts ---------------------------------------------------

table(d.curated.sum$AMG[d.curated.sum$skim_plot_isViral == "TRUE"])
# Likely viral Possible viral Unlikely viral 
#           30              8             18 

# sum(!is.na(d.curated$skim_plot_isViral))

sum(d.scafs$inspected)
#6784
mean(d.scafs$inspected)
# 119.0175
sd(d.scafs$inspected)
# 56.89825

# CSU sets ----------------------------------------------------------------

# summarize scaffolds
d.csu.sum <- d.curated %>% 
  # use gut viromes
  filter(set_group == "CSUsets") %>% 
  group_by(focal_gene, skim_plot_isViral) %>% 
  summarise(n_scaffolds = n()) %>% 
  # to add 0 where needed
  pivot_wider(names_from = skim_plot_isViral, 
              values_from = n_scaffolds, 
              values_fill = 0) %>% 
  pivot_longer(-1, names_to = "skim_plot_isViral", values_to = "n_scaffolds" ) %>% 
  mutate(gene = str_remove(focal_gene, ".*_"))

# add AMG calls from Guts
d.csu.sum <- d.csu.sum %>% 
  rename(csu.n_scaffolds = n_scaffolds) %>% 
  left_join(d.curated.sum, .)

#plotting order
d.csu.sum <- d.csu.sum %>% 
  mutate(gene = factor(gene, levels = genes))

csu.scafs <-
  d.csu.sum %>% 
  select(-n_scaffolds) %>% 
  rename(n_scaffolds = csu.n_scaffolds ) %>% 
  pivot_wider(names_from = skim_plot_isViral, values_from = n_scaffolds) %>% 
  mutate(inspected = `FALSE`+`TRUE`+MAYBE,
         detected = `FALSE`+`TRUE`+MAYBE+ `NA`)

# plot summary
p1 <- csu.scafs %>% 
  filter(AMG == "Likely viral") %>%
  ggplot(aes(gene))+
  geom_col(aes(y= detected, fill = "detected"))+
  geom_col(aes(y=inspected, fill = "inspected"))+
  
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "none")+
  scale_fill_manual(values = c("grey20", "red"), name=NULL)+
  scale_y_log10(breaks = 10^(0:4))+
  annotation_logticks(sides = "l")+
  ylab("N scaffolds")

# p1

p2 <- csu.scafs %>% 
  filter(AMG != "Likely viral") %>%
  ggplot(aes(gene))+
  geom_col(aes(y= detected, fill = "detected"))+
  geom_col(aes(y=inspected, fill = "inspected"))+
  
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "bottom")+
  scale_fill_manual(values = c("grey20", "red"), name=NULL)+
  scale_y_log10(breaks = 10^(0:4))+
  annotation_logticks(sides = "l")+
  ylab("N scaffolds")

# p2

p_detected <- plot_grid(p1,p2, ncol=1, rel_heights = c(.8,1))


# plot scaffolds inspected ------------------------------------------------------


# plot 
p1 <- d.csu.sum %>% 
  filter(!(skim_plot_isViral == "NA")) %>%
  rename(is_AMG_Viral = skim_plot_isViral) %>% 
  group_by(gene) %>% 
  mutate(perc_scaffolds = 100*csu.n_scaffolds/sum(csu.n_scaffolds)) %>% 
  filter(AMG == "Likely viral") %>%
  ggplot(aes(gene, perc_scaffolds))+
  geom_col(aes(fill = is_AMG_Viral))+
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "none")+
  scale_fill_viridis_d(direction = -1)+
  ylab("% of inspected")

p2 <- d.csu.sum %>% 
  filter(!(skim_plot_isViral == "NA")) %>%
  rename(is_AMG_Viral = skim_plot_isViral) %>% 
  group_by(gene) %>% 
  mutate(perc_scaffolds = 100*csu.n_scaffolds/sum(csu.n_scaffolds)) %>% 
  filter(!(AMG == "Likely viral")) %>%
  ggplot(aes(gene, perc_scaffolds))+
  geom_col(aes(fill = is_AMG_Viral))+
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "bottom")+
  scale_fill_viridis_d(direction = -1)+
  ylab("% of inspected")
# labs(caption =  "panels reflect AMGs examples validated (viral if >5, not viral if 0)")
p_inspected <- plot_grid(p1,p2, ncol=1, rel_heights = c(.8,1))

p <- plot_grid(NULL, p_detected, NULL ,p_inspected, ncol = 1, 
               rel_heights = c(0.1,1,0.1,1), labels = c("(a)","","(b)"))
ggsave(filename = here("dram_metaG/plots","CSU_scaffold_counts.png"),
       plot = p, width = 8, height = 9, bg = "white")


# Map CSU validated genes to environment ------------------------------------
# ecosystem classification
d.eco <- read_csv(here(data_dir, "ecosystem_details.csv"))



csu.eco.sum <- 
  # add ecosytem type data 
  left_join(d.curated, d.eco) %>% 
  mutate(gene = str_remove(focal_gene, ".*_") %>% 
           factor(levels = genes[genes %in% str_remove(likely_viral_genes, ".*_")])) %>% 
  # only CSU
  filter(set_group == "CSUsets") %>% 
  # only likely-viral genes
  filter(focal_gene %in% likely_viral_genes) %>% 
  # only validated by manual inspection
  filter(skim_plot_isViral == "TRUE") %>% 
  # summarize
  group_by(gene, ecosystem_type) %>% 
  summarise(n.scaffolds = n())


p <- csu.eco.sum %>%
  ggplot(aes(gene, ecosystem_type))+
  geom_tile(aes(fill = n.scaffolds), color = "white")+
  geom_text(aes(label = n.scaffolds))+
  facet_grid(.~"Likely-viral")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "bottom")+
  scale_fill_viridis_b(direction = -1)+
  scale_x_discrete(drop=FALSE)

ggsave(filename = here("dram_metaG/plots","CSUsets_ecosystems.png"),
       plot = p, width = 8, height = 4, bg = "white")
# Export validated viral genes --------------------------------------------
# host functions of putative AMGs

# enriched KOs 
d.enriched <- read_csv(here(data_dir, "enrichment","spor_enriched.csv")) 
# naming scheme used in curation
d.enriched <- d.enriched %>% 
  mutate(gene_name = case_when(
    !is.na(bs) ~ str_remove(bs, " \\[.*") %>% paste0("Bs_", .),
    !is.na(cd) ~ str_remove(cd, " \\[.*")%>% paste0("Cd_", .),
    TRUE ~ paste0("NA_",gene_id))) %>% 
  #correction for a gene that has a suffix
  mutate(gene_name = str_remove(gene_name, "-.*"))
# add curation calls
d.enriched <- d.scafs %>% 
  select(gene_name = focal_gene, AMG) %>% 
  distinct() %>% 
  left_join(d.enriched,  by = "gene_name")

# B. subtilis gene functions from subtiwiki
sw.genes <- read_csv(here(data_dir,"TIDY_subtiwiki.gene.export.2022-10-27.csv"))

likely_viral_genes <- d.scafs %>% 
  filter(AMG == "Likely viral") %>% 
  pull(focal_gene) 

likely_viral <-
  d.enriched %>% 
  filter(AMG == "Likely viral") %>% 
  separate(gene_name, into = c("sp", "gene" ), remove = F)

likely_viral <-sw.genes %>% 
  filter(title %in% likely_viral$gene) %>% 
  left_join(likely_viral, ., by = c("gene" = "title"))

# export
likely_viral %>% 
  select(gene_name, gene_id, Bs_locus = locus,
         description, Function = `function`) %>% 
  #clean
  mutate(description = str_remove_all(description,"\\[wiki\\||\\[protein\\|")) %>% 
  mutate(description = str_remove_all(description,"]")) %>% 
  mutate(Function = str_remove_all(Function,"\\[wiki\\||\\[protein\\|")) %>% 
  mutate(Function = str_remove_all(Function,"]")) %>% 
  write_csv(here(data_dir,"viral_sporGenes.csv"))


