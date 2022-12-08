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
    n_scaffolds > 5 ~ "viral",
    n_scaffolds < 1 ~ "not viral",
    TRUE ~ "maybe viral"
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
  filter(AMG == "viral") %>%
  ggplot(aes(gene))+
  geom_col(aes(y= detected, fill = "detected"))+
  geom_col(aes(y=inspected, fill = "inspected"))+
  
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "none")+
  scale_fill_manual(values = c("grey20", "red"))+
  scale_y_log10(breaks = 10^(0:4))+
  annotation_logticks(sides = "l")+
  ylab("N scaffolds")

p1

p2 <- d.scafs %>% 
  filter(AMG != "viral") %>%
  ggplot(aes(gene))+
  geom_col(aes(y= detected, fill = "detected"))+
  geom_col(aes(y=inspected, fill = "inspected"))+
  
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "bottom")+
  scale_fill_manual(values = c("grey20", "red"))+
  scale_y_log10(breaks = 10^(0:4))+
  annotation_logticks(sides = "l")+
  ylab("N scaffolds")

# p2

p_detected <- plot_grid(p1,p2, ncol=1, rel_heights = c(.8,1))

# ggsave(filename = here("dram_metaG/plots","GUT_scaffolds_detected.png"), plot = p, width = 8, height = 6)

# # different way of plotting
# p <- d.scafs %>% 
#   mutate(AMG = fct_rev(AMG)) %>% 
#   ggplot(aes(AMG, detected))+
#   geom_boxplot(outlier.colour = "white", fill = "grey80", width = 0.3)+
#   geom_jitter(height = 0, width = .1, shape = 21, fill = "white", size = 3)+
#   scale_y_log10()+
#   theme_classic(base_size = 16)+
#   ylab("N scaffolds detected")
# 
# # p
# ggsave(filename = here("dram_metaG/plots","GUT_AMG_counts.png"), plot = p, width = 4, height = 4)


# plot scaffolds inspected ------------------------------------------------------


# plot 
p1 <- d.curated.sum %>% 
  filter(!(skim_plot_isViral == "NA")) %>%
  rename(is_AMG_Viral = skim_plot_isViral) %>% 
  group_by(gene) %>% 
  mutate(perc_scaffolds = 100*n_scaffolds/sum(n_scaffolds)) %>% 
  filter(AMG == "viral") %>%
  ggplot(aes(gene, perc_scaffolds))+
  geom_col(aes(fill = is_AMG_Viral))+
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "none")+
  scale_fill_viridis_d(direction = -1)+
  ylab("% scaffolds examined")

p2 <- d.curated.sum %>% 
  filter(!(skim_plot_isViral == "NA")) %>%
  rename(is_AMG_Viral = skim_plot_isViral) %>% 
  group_by(gene) %>% 
  mutate(perc_scaffolds = 100*n_scaffolds/sum(n_scaffolds)) %>% 
  filter(!(AMG == "viral")) %>%
  ggplot(aes(gene, perc_scaffolds))+
  geom_col(aes(fill = is_AMG_Viral))+
  facet_grid(.~fct_rev(AMG), scales = "free_x",space = "free_x")+
  theme_classic()+
  panel_border(color = "black")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust=1),
        legend.position = "bottom")+
  scale_fill_viridis_d(direction = -1)+
  ylab("% scaffolds examined")
  # labs(caption =  "panels reflect AMGs examples validated (viral if >5, not viral if 0)")
p_inspected <- plot_grid(p1,p2, ncol=1, rel_heights = c(.8,1))

p <- plot_grid(NULL, p_detected, NULL ,p_inspected, ncol = 1, 
               rel_heights = c(0.1,1,0.1,1), labels = c("(a)","","(b)"))
ggsave(filename = here("dram_metaG/plots","GUT_scaffold_counts.png"),
       plot = p, width = 8, height = 9, bg = "white")



# summarize counts ---------------------------------------------------

# table(d.curated.sum$AMG[d.curated.sum$skim_plot_isViral == "TRUE"])
# sum(!is.na(d.curated$skim_plot_isViral))



