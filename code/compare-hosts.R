# setwd("~/GitHub/spore_amg")
library(here)
library(tidyverse)
library(ggrepel)

# import data from DRAM-v, w. host taxonomy

d.vir <- read_csv( here("data/Viruses/amg_summary_wHost.tsv")) %>% 
  # remove viruses of unknown hosts
  filter(!is.na(order))

# Firmicute or not
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

#resolve at family level
firmi <- d.vir %>%
  filter(firmicute.host) %>% 
  separate(family.etc, into = c("family", "genus", "etc"), sep = ";",extra = "merge", remove = F) %>% 
  group_by(phylum, class, family, genus) %>%
  summarise(n=n())

# classification of sporulators by Galperin 2013
spore.fam <- read_csv(here("data/Galperin_2013_MicrobiolSpectrum_table2.csv"))
# add likelihood of sporulation
# Galperin table fraction footnote:
  # "The distribution of sporeformers among the 
  # experimentally characterized members of the respective
  # family is indicated as follows: 
  #   +++, all (or nearly all) characterized members of the family produce spores;
  #   ++, a significant fraction of species are sporeformers; 
  #   +, the family includes some sporeformers; 
  #   -, no known sporeformers in the family."
spore.fam <- spore.fam %>% 
  mutate(spore.likley = case_when(Fraction.spore.forming == "+++" ~ TRUE,
                                  Fraction.spore.forming == "++" ~ TRUE,
                                  Fraction.spore.forming == "+" ~ FALSE,
                                  Fraction.spore.forming == "-"  ~ FALSE)) 

#resolve at family level
d.vir <-
  d.vir %>%
  separate(family.etc, into = c("family", "genus", "etc"), sep = ";",extra = "merge", remove = F) 

#add Galperin sporrulation 
d.vir <- spore.fam %>% 
  select(Family, spore.likley) %>% 
  left_join(d.vir, ., by = c("family" = "Family")) %>% 
  mutate(spore.likley = if_else(is.na(spore.likley), FALSE, spore.likley))
# %>% group_by(firmicute.host, spore.likley) %>% summarise(n=n())

d.vir.sum <- d.vir %>%
  # filter(str_detect(header, regex("sporulation", ignore_case = T))) %>%
  filter(!is.na(gene_id)) %>% 
  group_by(gene_id, gene_description, category,spore.likley) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# spore.fam$Family %in% fams

#### Enrichment test ####

n.phages <- d.vir %>% 
  # keep one record of each phage
  filter(!duplicated(virus.tax.id)) %>% 
  group_by(spore.likley) %>% 
  summarise(n=n())

d.stat <- d.vir.sum %>% 
  mutate(spore.likley = if_else(spore.likley, "sporulator.host", "nonsporulator.host")) %>% 
  pivot_wider(names_from = "spore.likley", values_from = "n", values_fill = 0)

d.vir %>% 
   filter(!is.na(gene_id)) %>%
     filter(!duplicated(virus.tax.id)) %>% 
     group_by(spore.likley) %>% 
     summarise(n=n())


# q	: the number of white balls drawn without replacement from an urn which contains both black and white balls.
  # AMG in viruses infecting firmicutes
# m	:the number of white balls in the urn.
  # number of viruses infecting firmicutes in the pool
# n	:the number of black balls in the urn.
  # number of viruses infecting non-firmicutes in the pool
# k: the number of balls drawn from the urn, hence must be in 0,1,., m+n.
  # Total AMG detected in viruses

d.stat <- d.stat %>% 
  mutate(M = n.phages$n[n.phages$spore.likley],
         N = n.phages$n[!n.phages$spore.likley],
         K = nonsporulator.host + sporulator.host) %>% 
  mutate(p.val = phyper(q = sporulator.host,
                        m = M, n = N, k = K, lower.tail = F, log.p = F))

d.stat$adj.p <- p.adjust(d.stat$p.val, method = "BH")
d.stat$sig <- d.stat$adj.p < 0.05         

# d.stat %>% 
#   filter(adj.p > 0) %>% 
#   pull(adj.p) %>% 
#   min()

# add labels to genes significantly enriched with K>10
d.stat <- d.stat%>% 
  mutate(print.lab = adj.p<0.05 & K > 10) %>% 
  mutate(lab = if_else(print.lab, gene_description,"")) %>% 
  mutate(lab = str_replace(lab, "dUTP pyrophosphatase", "dUTP pyrophosphatase;")) %>% 
  mutate(lab = str_extract(lab, pattern = ".*;")) %>% 
  mutate(lab = str_remove(lab, ";"))

# str_extract(d.stat$lab, pattern = regex(".*;|.*["))
d.stat %>% 
  ggplot(aes(K,adj.p+1e-24))+
  geom_hline(yintercept = 0.05, linetype=3, color="grey")+
  geom_vline(xintercept = 8, linetype=3, color="grey")+
  geom_jitter(width = 0.05, height = 0.05,
              shape=21,fill="grey", size=3, stroke = 2, alpha = 0.5)+
  geom_text_repel(aes(label = lab), na.rm = T )+
  scale_y_log10()+
  # scale_x_log10()+
  scale_x_continuous(trans = "log2", breaks = (2^(1:8)))+
  theme_classic()+
  xlab("No. homologs detected (log2)")+
  ylab("P-value (log10 adj. BH)")+
  labs(title = "Enrichment of phages with hosts likely to sporulate",
       subtitle = "calculated by hypergeometric distribution",
       caption = paste(d.stat$M[1],"phages with hosts likely to sporulate\n",
                       d.stat$N[1],"phages with hosts NOT likely to sporulate,"))
  ggsave(here("plots/enrichmnent_likely_sporultor.png"), width = 6, height = 6)
  
##########
  # stopped here
##########

# x <- 0:100
# y <- dhyper(x,303,1044,280)
# qplot(x,y)+geom_line()+geom_vline(xintercept = 54)
# phyper(q = 9,m = 303,n = 1044,k = 10, lower.tail = F, log.p = T)


#focus on spoIIIE
d.vir %>% 
  filter(str_detect(gene_description, "spoIIIE")) %>% 
  separate(family.etc, into = c("family", "genus", "etc"), extra = "merge", fill = "right") %>% 
  group_by(phylum,class, order, family, genus) %>% 
  summarize(n=n())

d.vir %>% 
  filter(str_detect(gene_description, "spoIIIE")) %>% 
  pull(virus.name)
