library(here)
library(tidyverse)

# import data from DRAM-v, w. host taxonomy

d.vir <- read_csv( here("data/Viruses/amg_summary_wHost.tsv")) %>% 
  # remove viruses of unknown hosts
  filter(!is.na(order))

# Firmicute or not
# d.vir %>% 
#   group_by(phylum, class) %>% 
#   summarise()
d.vir <- d.vir %>% 
  mutate(firmicute.host = str_detect(phylum, "Firmicutes"))

d.vir.sum <- d.vir %>%
  # filter(str_detect(header, regex("sporulation", ignore_case = T))) %>%
  filter(!is.na(gene_id)) %>% 
  group_by(gene_id, gene_description, category,firmicute.host) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

d.vir.sum %>% 
  ggplot(aes(gene_id, n, color = firmicute.host))+
  geom_point(shape=21, size=2)+
  facet_wrap(~category, scales = "free_y")+
  coord_flip()+
  theme_bw()+
  ggsave(here("plots", "AMG_x_hostTax.png"), width = 14, height = 7)
