# Plotting results of enrichment analysis for AMG in phages of sporulating hosts
# combining both genomes and metagenomes

library(here)
library(tidyverse)
library(cowplot)
library(ggrepel)



# import and arrange data -------------------------------------------------

# refseq enrichement
d.refseq <- read_csv(here("enrichment/data","enrichment_stats.csv")) %>% 
  mutate(dataset = "RefSeq")


# # GVD & GPD enrichment
d.gut <- read_csv(here("dram_metaG/data/", "enrichment","gvd_gpd_enrich.csv")) %>%
  mutate(dataset = "Gut Viromes")


# GVD & GPD inspection calls
d.AMG <- read_csv(here("dram_metaG/data/", "Gut_enriched_inspected.csv")) 
  #add to enrichment
d.gut <- d.AMG %>% 
    select(gene_id, gene, AMG) %>% 
    left_join(d.gut, .)

# combine
d <- bind_rows(d.refseq, d.gut) %>% 
  #fix typos
  mutate(spor_gene = str_replace_all(spor_gene, "Spoulation", "Sporulation") %>% 
           str_replace_all( "gene$", "genes")) %>% 
  #rename "other genes"
  mutate(spor_gene = str_replace_all(spor_gene, "Other", "Non sporulation")) %>% 
  # non-inspected genes
  mutate(AMG = case_when(
    !is.na(AMG) ~ AMG,
    dataset == "Gut Viromes" ~ "not inspected",
    dataset == "RefSeq" ~ "phage isolate"
  ))

# mark significantly enriched genes
p.signif <- 1e-6
k_treshold <- 30
d <- d %>% 
  mutate(significant = (p.adj < p.signif) & (k>k_treshold))

# labels
d <-
  d %>% 
  mutate(p_lab = case_when(
    !(significant) ~ "",
    dataset == "RefSeq" ~ p_lab,
    k < 200 | p.adj > 1e-20 ~ "",
    (dataset == "Gut Viromes") & AMG == "Likely viral" ~ gene,
    TRUE ~ ""
  )) 


# plot --------------------------------------------------------------------


p <-
  d %>% 
  filter(k>0) %>%
  ggplot(aes(k,-log10(p.adj+1e-120)))+
  
  
  geom_hline(yintercept = -log10(p.signif), linetype=2, color="grey", size = 1)+
  geom_vline(xintercept = k_treshold, linetype=2, color="grey", size = 1)+
  
  geom_jitter(aes(fill = significant),#, shape = AMG),
              width = 0.01, height = 0,
              shape=21,size=2, stroke = 0.3, alpha = 0.5, show.legend = F)+
  
  geom_text_repel(aes(label = p_lab), max.overlaps = 50,color = "black", size = 3,
                  force = 15, ylim=c(NA,120))+
  
  scale_fill_viridis_d(direction = -1, alpha = 0.5)+
  scale_shape_manual(values = 21:27)+
  theme_classic(base_size = 16)+
  facet_grid(dataset%>% fct_rev() ~ spor_gene %>% fct_rev()) +
  panel_border(color = "black")+
    scale_x_log10()+
    annotation_logticks(base = 10, sides = "b", outside = T, colour = "black",
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm"))+
    coord_cartesian(clip = "off")+
  xlab("Homologs detected")+
  ylab("Enrichment in phages\ninfecting spore-forming hosts")+
  theme(axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(here("dram_metaG/plots/", "enrichment_combined.png"),
       plot = p, width = 8, height = 6)
   
