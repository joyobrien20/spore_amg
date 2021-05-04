library(here)
library(tidyverse)


d.vir <- read_tsv(here("data/Viruses/amg_summary.tsv"))



d.vir.sum <- d.vir %>% 
  # filter(str_detect(header, regex("sporulation", ignore_case = T))) %>% 
  filter(!is.na(gene_id)) %>% 
  group_by(gene_id, gene_description, header) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))

#### match to hosts ####

#import virus-host data
#?# downloaded from ftp://ftp.genome.jp/pub/db/virushostdb/ (24/Nov/2020)
vh.db <- read_tsv(here("data","virushostdb.tsv") ) %>% 
  rename_all(make.names)

d.vir <- d.vir %>% 
  # parse viral ID to match
  mutate(refseq.id = str_remove(gene, ".*sequences.simple_header_")) %>% 
  mutate(refseq.id = str_remove(refseq.id, "\\..*$")) 

vh.db%>%
  filter(str_detect(host.lineage,regex("bacteria",ignore_case = T)))%>%
  group_by(virus.tax.id, virus.name)%>%
  summarise(n=n())%>%
  ggplot(aes(x=n))+
  geom_histogram()+
  # scale_y_log10()+
  ggtitle("VHDB host number for phages" )

#  virus duplicates in VHDB data
# These reflect multiple hosts
duplicated(vh.db$`virus tax id`)%>%sum() # 3461
