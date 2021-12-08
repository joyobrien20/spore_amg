library(here)
library(tidyverse)

# import data on hosts in refseq
d.vir <- read_tsv(here("data/Viruses/vMAG_stats.tsv"))


#### match to hosts ####

# #-----------------------------#
# # Downloading virus-host database data
# #-----------------------------#
# # downloaded on 4/May/2021
# library(curl)
# url <- "ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv"
# vh.db <- read_tsv(curl(url))
# 
# vh.db <- rename_all(vhdb, make.names)
# 
# write_tsv(vh.db, here("data","virushostdb.tsv"))
# #-----------------------------#

#import virus-host data
vh.db <- read_tsv(here("data","virushostdb.tsv") ) 

# It will be easier to match data by taxid
# extracting from vh.db a mapping-table between refseq and taxid
# there are many viruses that have multiple refseq ids (=multiple sequences)
# and there are viruses that have multiple hosts (= multiple rows)
map_taxid_refseq <- vh.db %>% 
  select(virus.tax.id, refseq.id) %>% 
  # get all sequence IDs
  separate(refseq.id, into = as.character(1:200), sep = ",",fill = "right") %>% 
  pivot_longer(-virus.tax.id, names_to = "num", values_to = "refseq.id") %>% 
  filter(! is.na(refseq.id)) %>% 
  select(-num) %>% 
  # discard replicated ID pairs
  distinct() %>% 
  mutate(refseq.id = trimws(refseq.id))


d.vir <- d.vir %>% 
  # parse viral ID to match
  mutate(refseq.id = str_remove(refseq.id, "\\..*$")) %>% 
  mutate(refseq.id = trimws(refseq.id))

# # Are we missing any viruses in vh.db
# sum(!d.vir$refseq.id %in% map_taxid_refseq$refseq.id) # Only 3! 
# d.vir$refseq.id[which(!d.vir$refseq.id %in% map_taxid_refseq$refseq.id)] #"NC_029050"
# "NC_029050" "NC_029072" "NC_042059"
#add taxon ID to viral amg data
d.vir <-left_join(d.vir, map_taxid_refseq, by = "refseq.id") 

#### dealing with multiple rows(-hosts) per virus ####
# keep only vh.db data relavent to amg data
vh.db <- semi_join(vh.db, d.vir, by = c("virus.tax.id"))
# how much duplication is there?
vh.db%>%
  group_by(virus.tax.id, virus.name)%>%
  summarise(N=n())%>%
  group_by(N) %>% 
  summarise(nN=n()) %>% 
  arrange(desc(N))
# 233 phages with 2 hosts and 20 phages with 3-6 hosts

# how much variation is there in multiple host taxonomy ?
# Typically phages will infect very closely related hosts ( sibling strains or species)

# First we will break up the host lineage data so we can check if  hosts are of simmilar taxonomy

# the host lineage is not uniform having 3-28 levels
# str_count(vh.db$host.lineage,";")%>%range(na.rm = TRUE) #0 28 => 1-29 levels
# str_count(vh.db$host.lineage,";")%>%hist()
# The 29 level is a single case of a crASS-like phage (IAS virus) that has humans listed as host
# should be excluded
vh.db <-  vh.db %>% filter(!str_detect(host.lineage, " Homininae; Homo"))
# # this led me to find also other viruses that are not phages
# i.e. host is not bacteria 
vh.db <-  vh.db %>% filter(str_detect(host.lineage, "Bacteria"))

# going back to host lineage:
str_count(vh.db$host.lineage,";")%>%range(na.rm = TRUE) #0 9 => 1-10 levels
# 3 have only bacteria listed as host:
# Halocynthia phage JM-2012
# Hamiltonella virus APSE1
# uncultured crAssphage
# removing these as well
vh.db <-  vh.db %>% filter(!str_detect(host.lineage, "Bacteria$"))

# going back to host lineage levels:
str_count(vh.db$host.lineage,";")%>%range(na.rm = TRUE) #1 9 => 2-10 levels

# the first looks like this: 
# Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Alcaligenaceae; Achromobacter
# Which corresponds to 
# Domain; Phylum; class; order, family, genus.

# but some have an intermediate rank between domain and phylum
# specifically FCB group AND Terrabacteria group 
# I will remove those for consistency

vh.db <- vh.db%>%
  mutate(host.lineage=str_remove(host.lineage," FCB group;"))%>%
  mutate(host.lineage=str_remove(host.lineage," Terrabacteria group;"))

vh.db <- vh.db%>%
  separate(host.lineage, sep="; ",
           into = c("domain", "phylum", "class", "order", "family.etc"), 
           extra = "merge", remove = F)

# Are there any viruses with hosts of different classes?
vh.db%>%
  group_by(virus.tax.id, virus.name, class)%>%
  summarise(n=n())%>%
  arrange(desc(n)) %>% 
  group_by(virus.tax.id, virus.name)%>%
  summarise(n=n())%>%
  filter(n>1) %>% 
  arrange(desc(n))

# There is only one such virus,  Thermus phage phi OH2
# PMID 23950135: this is a prophage of a Gebacillus strain
# vhdb comment: This virus is now called "Bacteriophage Phi OH2". No longer any mention of "Thermus" host mentioned in the sequence record AB823818. 
# dropping that host record

vh.db <- vh.db%>%
  filter(! ((virus.name == "Thermus phage phi OH2") & 
           str_detect(host.lineage, "Thermus")))


# all hosts for each phage are the same down to class
# some phages have refseq in evidence for one of their hosts
# I think this is the isolation host 

#remove host duplicate
vh.db <- vh.db %>% 
  # mark refseq evidence as frist priority
  mutate(priority1 = grepl("RefSeq", .$evidence)) %>% 
  # count host taxonomy levels, for second priority
  mutate(priority2 = str_count(host.lineage,";")) %>% 
  arrange(desc(priority1), desc(priority2)) %>% 
  # remove duplicated
  mutate(dup = duplicated(virus.tax.id)) %>% 
  filter(!dup) %>% 
  select(-priority1, -priority2, -dup)

# duplicated(vh.db$virus.tax.id)%>%sum() # 0


d.vir <- left_join(d.vir, vh.db, by = "virus.tax.id")




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

#add Galperin sporulation 
d.vir <- spore.fam %>% 
  select(Family, spore.likley) %>% 
  left_join(d.vir, ., by = c("family" = "Family")) %>% 
  mutate(spore.likley = if_else(is.na(spore.likley), FALSE, spore.likley))


# d.vir.sum <- d.vir %>%
#   mutate(firmicute = str_detect(phylum, "Firmicutes")) %>% 
#   filter(!is.na(firmicute)) %>% 
#   group_by(firmicute,spore.likley) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))

write_csv(d.vir, here("data/Viruses/vMAG_host_sporul.csv"))
