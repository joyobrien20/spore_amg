# summarize genes by host

library(here)
library(tidyverse)

# Data sets -------------------------
data_dir <-"dram_metaG/data/dram_output"
sets <- list.dirs(here(data_dir), full.names = F, recursive = F)

# remove extra data
sets <- sets[sets!="extra_data"]

# LOOP starts here ---------

set <- sets[12]

# > host sporulation assignment  ----------------
  # execute separate code that does this
  # loads into the environment spor_mags df
  #relies on the "set" variable
source(here("dram_metaG/code/A2_assign_sporFam.R"))


# > read in AMG table ----------------

d.amg <- read_tsv(here(data_dir, set, "amg_summary.tsv"))

  # add host sporulation
d.amg_host <- 
  spor_mags %>% 
  select(scaffold = `Virus ID`, Host_taxonomy = gtdb_dpf,
         spore_likely = f_spor) %>% 
  left_join(d.amg, . , by = "scaffold") 

p <- spor_mags %>% 
  
  ggplot(aes(`Virus ID`, Score))+
  geom_point(aes(fill=f_spor),shape=21)+coord_flip()+theme_bw()+
  theme(axis.text.y = element_blank())
  ggsave(here("dram_metaG/plots/","methyl_host.png"), p, width = 8, height = 20)
  
  top <- spor_mags %>% 
    group_by(`Virus ID`) %>% 
    summarise(n=n()) %>%
    filter(n>=80) %>% pull(`Virus ID`)

  spor_mags %>% 
    filter(`Virus ID` %in% top) %>% 
    ggplot(aes(Score))+
    geom_histogram(aes(fill = f_spor))+
    facet_wrap(~ `Virus ID`)+
    theme_classic()+
    # theme(#strip.background = element_blank(),
    #       strip.text.x = element_blank())+
    cowplot::panel_border(color = "black")
  
    
  mutate(`Virus ID` = fct_inorder(`Virus ID`)) %>% 
    ggplot(aes(`Virus ID`, n)) +
    geom_col() +
    coord_flip()+
    theme_classic()+
    theme(axis.text.y = element_blank())


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
