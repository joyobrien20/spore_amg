library(here)
library(tidyverse)

# downloaded gpd scaffold meta with predicted host taxonomy
# download.file(
#   url = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/GPD_metadata.tsv",
#   destfile = here("dram_metaG/data/dram_output/Camarillo-Guerrero/GPD_metadata.tsv")
# )

# from README
# * GPD_metadata.tsv: Metadata of the GPD. 
# 
# It contains the following columns:
#   
# 1) GPD_id: Genome identifier
# 2) Source: Indicates whether the genome was predicted from a metagenome or a bacterial isolate
# 3) GPD_VC: Viral cluster to which the genome belongs
# 4) Size: Genome length in bp
# 5) Predicted_phage_taxon: Predicted viral taxon of the genome
# 6) Host_range_isolates: Host assignment by CRISPR exact matching and prophage matching
# 7) Host_range_taxon: GTDB taxonomy assignment of Host_range_isolates
# 8) Metagenomic_runs_detected: Human gut metagenomes in which the genome was detected by read mapping
# 9) Continents_detected: Continent of origin of Metagenomic_runs_detected
# 10) Countries_detected: Country of origin of Metagenomic_runs_detected
# 11) CheckV_MIUVIG: MIUVIG quality predicted by checkV
# 12) CheckV_completion: Genome completion predicted by checkV
# 13) CheckV_viral_region: Viral regions predicted by checkV
# 14) CheckV_host_region: Host regions predicted by checkV
# 15) CheckV_prophage: Indicates whether the genome is predicted to be a prophage or not
# 16) CheckV_termini: Indicates whether the genome contains a direct terminal repeat
# 17) CheckV_copies: Number of genome copies detected by checkV
# 18) Novel: Indicates whether the genome is exclusive to GPD or not

gpd_meta <- read_tsv(here("dram_metaG/data/dram_output/Camarillo-Guerrero/GPD_metadata.tsv"))

gpd_meta %>% 
  filter(!is.na(Host_range_taxon)) %>% 
  mutate(n_pred = str_count(Host_range_taxon, ",")+1) %>% 
  # ggplot(aes(n_pred))+geom_histogram()
  group_by(n_pred) %>% 
  summarise(n=n()) %>% 
  arrange(n_pred) %>% 
  mutate(cum_n = cumsum(n)) %>% 
  mutate(cum_perc = 100 * cum_n/sum(n)) #%>% 
  # ggplot(aes(n_pred, cum_perc))+
  # geom_line()+
  # # scale_y_log10()+
  # theme_bw()

# range of predictions 1-46, with 90% having 8 or less.


# parse host taxonomy -----------------------------------------------------

# Levels of taxonomy
gtdb_tax <- c("gtdb_p","gtdb_c","gtdb_o","gtdb_f","gtdb_g","gtdb_s")

gpd_hosts <- 
  gpd_meta %>% 
  #keep only phages with predicted hosts
  filter(!is.na(Host_range_taxon)) %>% 
  #transform to long format
  select(GPD_id, Host_range_taxon) %>% 
  separate(Host_range_taxon, into = paste0("host", 1:50), sep = ",", fill = "right") %>% 
  pivot_longer(-GPD_id, names_to = "host_n", values_to = "Host_range_taxon") %>% 
  filter(!is.na(Host_range_taxon)) %>% 
  #parse taxonomy
  separate(Host_range_taxon, into = gtdb_tax, sep = "/")



# Host families per phage -------------------------------------------------

# Since I make calls of host sporulation by family, I wonder how many 
# of the Firmicutes phages have predictions that span different host families
gpd_hosts %>% 
  mutate(fam = str_c(gtdb_p, gtdb_c, gtdb_o, gtdb_f, sep = "-")) %>% 
  group_by(GPD_id, gtdb_p, fam) %>% 
  summarise() %>% 
  ungroup() %>% 
  group_by(GPD_id, gtdb_p) %>% 
  summarise(n_fams=n()) %>% 
  ungroup() %>% 
  group_by(gtdb_p, n_fams) %>% 
  summarise(n = n())
# only 75 of 25651 Firmicutes phages (40932 total phages)
# those 75 have only 2 families, never more
# (in other phyla some have 3 families)


gpd_hosts %>%
  filter(gtdb_p == "Firmicutes") %>% 
  mutate(fam = str_c(gtdb_p, gtdb_c, gtdb_o, gtdb_f, sep = "-")) %>% 
  group_by(fam) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))
  pull(fam) %>% 
  anyDuplicated()

# ===================================
# old code (unmodified) below this
  
d.gvd <- read_csv(here("metaG/data/gvd/gvd_all_hosts.csv"))
d.f_spor <- read_csv(here("gtdb_spor/data/gtdb_families_sporulation.csv"))


# make gvd-like taxonomy
d.f_spor <- d.f_spor %>% 
  mutate(GTDB_Host = str_c(str_remove(gtdb_d, "d__"), gtdb_p,gtdb_f, sep = ";"))

#clean gvd
d.gvd <- d.gvd %>% 
  filter(Eukaryotic_or_Prokaryotic_Virus == "Bacteriophage") %>% 
  filter(!str_detect(GTDB_Host, "Host Not Assigned" )) %>% 
  filter(str_detect(GTDB_Host,"p__Firmicutes"))

# compatability check
d.gvd %>% 
  filter(!str_detect(GTDB_Host, "Host Not Assigned" )) %>% 
  pull(GTDB_Host)
           
table(d.gvd$GTDB_Host %in% d.f_spor$GTDB_Host, useNA = "a")
# FALSE  TRUE  <NA> 
#   301  6623     0 

d.gvd %>% 
  filter(!d.gvd$GTDB_Host %in% d.f_spor$GTDB_Host) %>% 
  pull(GTDB_Host) %>% 
  unique() %>% 
  length() #18


 # the other way
table(d.f_spor$GTDB_Host %in% d.gvd$GTDB_Host,useNA = "a")
# FALSE  TRUE  <NA> 
#   336    66     0

d.f_spor %>% 
  filter(d.f_spor$GTDB_Host %in% d.gvd$GTDB_Host) %>% 
  pull(f_spor) %>% 
  table(useNA = "a")
# FALSE  TRUE  <NA> 
#   37    22     7 


# There are 18 Firmicutes families missing from my list,
# accounting for 301 phages.
# my list accounts for 66 Firmicutes families in GVD 
# for which I have sporulation prediction for all but 7.

#reload gvd and join
d.gvd <- read_csv(here("metaG/data/gvd/gvd_all_hosts.csv"))

d.gvd <- d.f_spor %>% 
  select(GTDB_Host, f_spor) %>% 
  left_join(d.gvd, ., by = "GTDB_Host" )



#
  d.gvd %>% 
  filter(Eukaryotic_or_Prokaryotic_Virus == "Bacteriophage") %>% 
  filter(!str_detect(GTDB_Host, "Host Not Assigned" )) %>% 
    mutate(f_spor = if_else(
      str_detect(GTDB_Host,regex("p__Firmicutes", ignore_case = T)),
      f_spor, FALSE)) %>% 
    group_by(f_spor) %>% 
    summarise(n=n())
  # f_spor     n
  # <lgl>  <int>
  
  # 1 FALSE   9383
  # 2 TRUE    4224
  # 3 NA       326

#save gvd with sporulation prediction 
gvd_spor <- d.gvd; rm(d.gvd) 
save(gvd_spor,file = here("metaG/data/gvd/gvd_spor.Rdata"))

# delete gvd download 
unlink(here("metaG/data/gvd/gvd_all_hosts.csv"))
