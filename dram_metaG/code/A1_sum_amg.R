# summarize genes by host

library(here)
library(tidyverse)

# Data sets -------------------------
data_dir <-"dram_metaG/data/dram_output"
sets <- list.dirs(here(data_dir), full.names = F, recursive = F)

# remove extra data
sets <- sets[sets!="extra_data"]

# LOOP starts here ---------

set <- sets[2]

# > host sporulation assignment  ----------------
  # execute separate code that does this
  # loads into the environment spor_mags df
  #relies on the "set" variable
source(here("dram_metaG/code/A_assign_sporFam.R"))


# > read in AMG table ----------------

d.amg <- read_tsv(here(data_dir, set, "amg_summary.tsv"))

  # add host sporulation
d.amg_host <- 
  spor_mags %>% 
  select(scaffold = `Virus ID`, Host_taxonomy = gtdb_dpf,
         spore_likely = f_spor) %>% 
  left_join(d.amg, . , by = "scaffold") 

spor_mags %>% 
  
  ggplot(aes(`Virus ID`, Score))+
  geom_point(aes(fill=f_spor),shape=21)+coord_flip()+theme_bw()


  # remove genes for which host is unknown
  filter(!str_detect(Host_taxonomy,"Host Not Assigned")) %>% 
  # Remove genes dropped from Gregory
  filter(!is.na(Host_taxonomy)) %>% 
  #remove genes for which there is uncertainty on host sporuation
  filter(!is.na(spore_likely)) 

# # summarize
# spor_mags_sum <- 
#   spor_mags %>%
#   group_by(f_spor) %>%
#   summarise(n=n())
