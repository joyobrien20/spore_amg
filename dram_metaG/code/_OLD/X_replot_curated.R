library(here)
library(tidyverse)
library(cowplot)



# get data from reviewed AMGs ----------------------------------------------

# manual curation results are in "curated" folder
data_dir <- ("dram_metaG/data/")
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

# summarize scaffolds
d.curated.sum <- d.curated %>% 
  group_by(focal_gene, skim_plot_isViral) %>% 
  summarise(n_scaffolds = n()) %>% 
  # to add 0 where needed
  pivot_wider(names_from = skim_plot_isViral, 
              values_from = n_scaffolds, 
              values_fill = 0) %>% 
  pivot_longer(-1, names_to = "skim_plot_isViral", values_to = "n_scaffolds" )

# make calls based on number of TRUE observations
d.curated.sum <- d.curated.sum %>% 
  filter(skim_plot_isViral=="TRUE") %>% 
  mutate(AMG = case_when(
    n_scaffolds > 9 ~ "viral",
    n_scaffolds < 1 ~ "not viral",
    TRUE ~ "maybe viral"
  )) %>% 
  # add calls
  select(focal_gene, AMG) %>% 
  left_join(d.curated.sum, by = "focal_gene") 

# plot scaffolds inspected ------------------------------------------------------

# arrange plotting order
genes <- d.curated.sum %>% 
  group_by(focal_gene) %>% 
  summarise(total = sum(n_scaffolds)) %>% 
  arrange(total %>% desc()) %>% 
  pull(focal_gene)

d.curated.sum <- d.curated.sum %>% 
  mutate(focal_gene = factor(focal_gene, levels = genes))


# plot 
p1 <- d.curated.sum %>% 
  filter(!(skim_plot_isViral == "NA")) %>%
  rename(is_AMG_Viral = skim_plot_isViral) %>% 
  group_by(focal_gene) %>% 
  mutate(perc_scaffolds = 100*n_scaffolds/sum(n_scaffolds)) %>% 
  filter(AMG == "viral") %>%
  ggplot(aes(focal_gene, perc_scaffolds))+
  geom_col(aes(fill = is_AMG_Viral))+
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "none")+
  scale_fill_viridis_d(direction = -1)+
  ylab("% scaffolds examined")

p2 <- d.curated.sum %>% 
  filter(!(skim_plot_isViral == "NA")) %>%
  rename(is_AMG_Viral = skim_plot_isViral) %>% 
  group_by(focal_gene) %>% 
  mutate(perc_scaffolds = 100*n_scaffolds/sum(n_scaffolds)) %>% 
  filter(!(AMG == "viral")) %>%
  ggplot(aes(focal_gene, perc_scaffolds))+
  geom_col(aes(fill = is_AMG_Viral))+
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "bottom")+
  scale_fill_viridis_d(direction = -1)+
  ylab("% scaffolds examined")+
  labs(caption =  "panels reflect AMGs examples validated (viral if >9, not viral if 0)")
p <- plot_grid(p1,p2, ncol=1, rel_heights = c(.8,1))
ggsave(filename = here("dram_metaG/plots","scaffolds_inspected.png"), plot = p, width = 8, height = 6)


# plot scaffolds detected ------------------------------------------------------

# plot summary
p1 <- d.curated.sum %>% 
  rename(is_AMG_Viral = skim_plot_isViral) %>% 
  filter(AMG == "viral") %>%
  ggplot(aes(focal_gene, n_scaffolds))+
  
  geom_col()+
  geom_col(data=d.curated.sum %>% filter(AMG == "viral") %>% filter (skim_plot_isViral != "NA"), fill = "red")+
  
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "none")+
  scale_fill_viridis_d(direction = -1)+
  ylab("N scaffolds detected")
  # ylim(0,1.1e4)


p2 <- d.curated.sum %>% 
  filter(!(skim_plot_isViral == "NA")) %>%
  rename(is_AMG_Viral = skim_plot_isViral) %>% 
  filter(!(AMG == "viral")) %>%
  ggplot(aes(focal_gene, n_scaffolds))+
  geom_col()+
  geom_col(data=d.curated.sum %>% filter(AMG != "viral") %>% filter (skim_plot_isViral != "NA"), fill = "red")+
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "bottom")+
  scale_fill_viridis_d(direction = -1)+
  ylab("N scaffolds detected")
  # ylim(0,1.1e4)
p <- plot_grid(p1,p2, ncol=1)
ggsave(filename = here("dram_metaG/plots","scaffolds_detected.png"), plot = p, width = 8, height = 6)

# different way of plotting
p <- d.curated.sum %>% 
  group_by(focal_gene,AMG) %>% 
  summarise(n_scaffolds = sum(n_scaffolds)) %>% 
  mutate(AMG = fct_rev(AMG)) %>% 
  ggplot(aes(AMG, n_scaffolds))+
  geom_boxplot(outlier.colour = "white", fill = "grey80", width = 0.3)+
  geom_jitter(height = 0, width = .1, shape = 21, fill = "white", size = 3)+
  scale_y_log10()+
  theme_classic(base_size = 16)+
  ylab("N scaffolds detected")
ggsave(filename = here("dram_metaG/plots","AMG_counts.png"), plot = p, width = 4, height = 4)

# summarize counts ---------------------------------------------------

table(d.curated.sum$AMG[d.curated.sum$skim_plot_isViral == "TRUE"])
sum(!is.na(d.curated$skim_plot_isViral))


# replot enrichment -------------------------------------------------------


# enriched KOs 
d.enriched <- read_csv(here(data_dir, "enrichment","spor_enriched.csv")) 

# naming scheme used in curation
d.enriched <- d.enriched %>% 
  mutate(gene_name = case_when(
    !is.na(bs) ~ str_remove(bs, " \\[.*") %>% paste0("Bs_", .),
    !is.na(cd) ~ str_remove(cd, " \\[.*")%>% paste0("Cd_", .),
    TRUE ~ paste0("NA_",gene_id))) 

# add curation calls
d.enriched <- d.curated.sum %>% 
  select(gene_name = focal_gene, AMG) %>% 
  distinct() %>% 
  left_join(d.enriched,  by = "gene_name")

#significance threshold
p.signif <- 1e-6
k_treshold <- 32

p <-
  d.enriched %>% 
  filter(k>0) %>%
  ggplot(aes(k,-log10(p.adj+1e-120)))+
  # geom_hline(yintercept = -log10(p.signif), linetype=2, color="grey", size = 1)+
  # geom_vline(xintercept = k_treshold, linetype=2, color="grey", size = 1)+
  
  # geom_text_repel(daaes(label = lab), max.overlaps=50,color = "grey", size = 3)+
  
  geom_jitter(aes(fill = AMG),width = 0.005, height = 0.005,
              shape=21, size=2, stroke = 1, alpha = 0.5, show.legend = T)+
  ggrepel::geom_text_repel(aes(label = gene_name, color = AMG), force = 20,
                           size = 3, show.legend = F, max.overlaps = 25)+
  scale_fill_viridis_d(direction = -1, alpha = 0.5)+
  scale_color_viridis_d(direction = -1, alpha = 1)+
  theme_classic()+
  facet_wrap(~ spor_gene %>% fct_rev(), ncol = 2) +
  panel_border(color = "black")+
  scale_x_continuous(trans = "log2", breaks = (2^(0:11)))+
  xlab("Sample size\nNo. homologs detected (log2)")+
  ylab("Enrichment (-log10 adj. P-value)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(filename = here("dram_metaG/plots","scrutinized_enrichment.png"), plot = p, width = 8, height = 8)

# host functions of putative AMGs -----------------------------------------

# B. subtilis gene functions from subtiwiki
sw.genes <- read_csv(here(data_dir,"TIDY_subtiwiki.gene.export.2022-10-27.csv"))

bs.true <- d.enriched %>% 
  filter(AMG == "viral") %>% 
  separate(gene_name, into = c("sp", "gene" ), remove = F)

bs.true <-sw.genes %>% 
  filter(title %in% bs.true$gene) %>% 
  left_join(bs.true, ., by = c("gene" = "title"))


# export
bs.true %>% 
  select(gene_name, gene_id, Bs_locus = locus,
         description, Function = `function`) %>% 
  #clean
  mutate(description = str_remove_all(description,"\\[wiki\\||\\[protein\\|")) %>% 
  mutate(description = str_remove_all(description,"]")) %>% 
  mutate(Function = str_remove_all(Function,"\\[wiki\\||\\[protein\\|")) %>% 
  mutate(Function = str_remove_all(Function,"]")) %>% 
  write_csv(here(data_dir,"viral_sporGenes.csv"))
