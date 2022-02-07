library(here)
library(tidyverse)
library(foreach)



d.amg <- read_csv(here("metag/data/gvd_enriched_fullData.csv"))
coat_ko <- c("K06333","K06334")
# focus on enriched genes
d.amg.coat <- 
  d.amg %>%
  filter(gene_id %in% coat_ko) %>% 
  # parse scaffold-gene
  select(gene_id,scaffold, gene, amg_flags) %>% 
  mutate(n.gene = str_remove(gene, ".*_") %>% parse_number()) %>% 
  select(-gene) 
# %>% 
#   group_by(scaffold) %>% 
#   summarise(genes = str_c(n.gene, collapse = ";"),
#             KOs = str_c(gene_id, collapse = ";"),
#             flags = str_c(amg_flags, collapse = ";"))# %>% 
# # filter(str_detect(genes, ";"))


# gene number per scaffold
# filtering function from chunked reading
f <- function(.data, pos) {
  filter(.data, scaffold %in% d.amg.coat$scaffold) 
}
# read annotations
d.annot <- read_tsv_chunked(here("metaG/data/gvd/annotations.tsv"),
                      DataFrameCallback$new(f), chunk_size = 10000)

d.annot <- d.annot %>% 
  group_by(scaffold) %>% 
  summarise(n.genes = max(gene_position)) %>% 
  arrange(n.genes)

d.amg %>%
  filter(scaffold %in% d.amg.coat$scaffold) %>% 
  # parse scaffold-gene
  select(gene_id,scaffold, gene, amg_flags) %>% 
  mutate(n.gene = str_remove(gene, ".*_") %>% parse_number()) %>% 
  select(-gene)  %>% 
  mutate(scaffold = fct_relevel(scaffold, d.annot$scaffold)) %>% 
  mutate(coat = if_else(gene_id %in% coat_ko, gene_id, "other amg")) %>% 
  ggplot(aes(x = n.gene , y = scaffold))+
  geom_tile(aes(fill = coat))+
  geom_segment(data = d.annot, aes(xend = n.genes, yend = scaffold), x = 1, color = "grey30") +
  scale_fill_manual(values = c("red", "blue", "grey70"))+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        legend.position = c(0.8,0.2))

ggsave(here("metaG/plots/gvd_coat_scaffolds.png"), height = 5, width = 4)
