# library(here)
# library(tidyverse)
# library(Biostrings)
# # BiocManager::install("ggmsa")
# # devtools::install_github("YuLab-SMU/ggmsa")
# 
# library(ggmsa)
# 
# aa_msa <- readAAMultipleAlignment(filepath = here("genes_of_interest/spo0A/data/align-trim-tree/seq2align_MafftEinsi.aln"),
#                                   format = "fasta")
# 
# # Plot all ----------------------------------------------------------------
# 
# d.aa_msa <- tidy_msa(aa_msa) %>% 
#   mutate(g = str_remove(name, "_.*")) %>% 
#   mutate(g = if_else(g == "kegg", "Bacteria", g))
# # mutate(g = if_else(name == "kegg_BSU24220", "Bs", g)) %>% 
# # mutate(g = fct_relevel(g, c("virome", "Bs", "kegg")))
# 
# 
# p <- ggplot() + 
#   geom_msa(data = d.aa_msa2 , 
#            seq_name = F, font = NULL, border = NA, color = 'Zappo_AA')+
#   facet_grid (fct_rev(g)~ ., scales = "free_y", space = "free_y")+
#   theme_classic()+
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title.y = element_blank())
# 
# ggsave(filename = here("genes_of_interest/spo0A", "plots", "msa1.png"), plot = p)
# 
# # Focused Plot ------------------------------------------------------------
# 
# 
# d.aa_msa <- tidy_msa(aa_msa, 250,400) %>% 
#   mutate(g = str_remove(name, "_.*")) %>% 
#   mutate(g = if_else(g == "kegg", "Bacteria", g))
# 
# p1 <- 
#   ggplot()+
#   geom_msa(data = d.aa_msa %>% 
#              filter(g == "virome") %>% droplevels(),
#            seq_name = F, font = NULL, border = NA,color = 'Zappo_AA')+
#   geom_seqlogo(top = T,color = 'Zappo_AA')+
#   facet_grid (fct_rev(g)~ ., scales = "free_y", space = "free_y")+
#   theme_classic()+
#   theme(axis.text = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title = element_blank())
# 
# 
# p2 <-    
#   ggplot()+
#   geom_msa(data = d.aa_msa %>% 
#              filter(g != "virome", ),
#            seq_name = F, font = NULL, border = NA,color = 'Zappo_AA')+
#   geom_seqlogo(top = T,color = 'Zappo_AA')+
#   facet_grid (fct_rev(g)~ ., scales = "free_y", space = "free_y")+
#   theme_classic()+
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title.y = element_blank())
# 
# 
# p3 <- cowplot::plot_grid(p1,p2, ncol = 1, rel_heights = c(1,2))
# 
# ggsave(filename = here("genes_of_interest/spo0A", "plots", "msa2.png"), plot = p3)
# 
# 
# 
# # Focused Plot 2 ------------------------------------------------------------
# 
# 
# d.aa_msa <- tidy_msa(aa_msa, 360,390) %>% 
#   mutate(g = str_remove(name, "_.*")) %>% 
#   mutate(g = if_else(g == "kegg", "Bacteria", g))
# 
# p1 <- 
#   ggplot()+
#   geom_msa(data = d.aa_msa %>% 
#              filter(g == "virome") %>% droplevels(),
#            seq_name = F, font = NULL, border = NA,color = 'Zappo_AA')+
#   geom_seqlogo(top = T,color = 'Zappo_AA')+
#   facet_grid (fct_rev(g)~ ., scales = "free_y", space = "free_y")+
#   theme_classic()+
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title.y = element_blank())
# 
# 
# p2 <-    
#   ggplot()+
#   geom_msa(data = d.aa_msa %>% 
#              filter(g != "virome", ),
#            seq_name = F, font = NULL, border = NA,color = 'Zappo_AA')+
#   geom_seqlogo(top = T,color = 'Zappo_AA')+
#   facet_grid (fct_rev(g)~ ., scales = "free_y", space = "free_y")+
#   theme_classic()+
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title.y = element_blank())
# 
# 
# p3 <- cowplot::plot_grid(p1,p2, ncol = 1)
# 
# ggsave(filename = here("genes_of_interest/spo0A", "plots", "msa3.png"), plot = p3)
# 
# 
# 
# 
