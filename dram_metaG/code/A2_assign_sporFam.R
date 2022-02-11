# accept list of MAG species and return sporulation assignment
library(here)
library(tidyverse)

# Curated list of sporulation in families of Firmicutes
fam_spor <- read_csv(here("gtdb_spor/data/gtdb_families_sporulation.csv")) %>% 
  distinct() %>% 
  mutate(gtdb_dpf = str_c(gtdb_d,gtdb_p, gtdb_f, sep = ";"))


#load list of species
mags <- read_csv(here("dram_metaG/data/dram_output", set, "mag_vmag_info.csv"))


mags <- mags %>% 
  mutate(gtdb_dpf = str_c(
    str_extract(Taxonomy_String, "d__.*?;") %>% str_remove(";"),
    str_extract(Taxonomy_String, "p__.*?;") %>% str_remove(";"),
    str_extract(Taxonomy_String, "f__.*?;") %>% str_remove(";"),
    sep = ";"))



spor_mags <- 
  fam_spor %>% 
  select(gtdb_dpf, f_spor) %>% 
  left_join(mags, .)

# assign non-Firmicutes as non-sporulators
spor_mags <-
  spor_mags %>% 
  mutate(f_spor = if_else(
    str_detect(gtdb_dpf,"p__Firmicutes"),
    f_spor, FALSE))

rm(mags, fam_spor)