library(here)
library(tidyverse)

# downloaded gpd scaffold meta with predicted host taxonomy
# download.file(
#   url = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/gut_phage_database/GPD_metadata.tsv",
#   destfile = here("dram_metaG/data/dram_output/Camarillo-Guerrero")
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
  
  
  

# gpd Vs GTDB sporulator list ---------------------------------------------

f_spor <- read_csv(here("gtdb_spor/data/gtdb_families_sporulation.csv"))  

# are Firmicutes family names unique?
  dup_fam <- duplicated(f_spor$gtdb_f) %>%
    which() %>% 
    f_spor$gtdb_f[.]
  # only one "f__Sporanaerobacteraceae"
  
  f_spor %>% 
    filter(gtdb_f %in% dup_fam) %>% view()
  #these are fully duplicated rows
  # removing duplicate
  f_spor <- distinct(f_spor)
  
  duplicated(f_spor$gtdb_f) %>% which() 
  # no duplicated families left
  
  # make sporulator list and gpd comparable  
  gpd_firmi_fam <- 
    gpd_hosts %>%
    filter(gtdb_p == "Firmicutes") %>% 
    mutate(fam = str_c("f__",gtdb_f, sep = "")) %>% 
    pull(fam) %>% 
    unique()

  # any gpd families missing from f_spor?
  gpd_firmi_fam [!gpd_firmi_fam %in% f_spor$gtdb_f]
  # "f__DTU089" "f__Helcococcaceae" "f__Bacillaceae_A" "f__NA"    
  
  #using GTDB to track these outdated taxon names:
  tax_update <- tibble(gpd_f = rep(NA,4), gtdb_f=rep(NA,4))
  
  # https://gtdb.ecogenomic.org/taxon_history/?from=R80&to=R202&query=f__DTU089
  tax_update[1,] <- list("f__DTU089", "f__Acutalibacteraceae")
  
  # https://gtdb.ecogenomic.org/taxon_history/?from=R80&to=R202&query=f__Helcococcaceae
  tax_update[2,] <- list("f__Helcococcaceae", "f__Peptoniphilaceae")

  # https://gtdb.ecogenomic.org/taxon_history/?from=R80&to=R202&query=f__Bacillaceae_A
  tax_update[3,] <- list("f__Bacillaceae_A", "f__DSM-18226")
  tax_update[4,] <- list("f__Bacillaceae_A", "f__DSM-1321")
    #this family has been split
  
  # Are the updated families listed in f_spor?
  f_spor %>% 
    select(gtdb_f, f_spor) %>% 
    left_join(tax_update, .)
  # yes!!
  
  

# Assign sporulation prediction to GPD hosts ------------------------------

  gpd_hosts <-
    f_spor %>% 
    # join update to f_spor
    select(gtdb_f, f_spor) %>% 
    mutate(gpd_f = gtdb_f) %>% 
    bind_rows(., tax_update) %>% 
    # join f_spor to gpd_hosts
    select (-gtdb_f) %>% 
    mutate(gpd_f = str_remove(gpd_f, "f__")) %>% 
    left_join(gpd_hosts, ., by = c("gtdb_f" = "gpd_f"))
  
  # assign non-sporulating to non-Firmicutes
  gpd_hosts <- 
    gpd_hosts %>% 
    mutate(f_spor = if_else(gtdb_p =="Firmicutes", f_spor, FALSE))
  
  

# Summarize sporulation for each scaffold ---------------------------------------------
  gpd_hosts %>% 
    group_by(GPD_id, f_spor) %>% 
    summarise() %>% 
    ungroup() %>% 
    group_by(GPD_id) %>% 
    summarise(n_spor=n()) %>% 
    pull(n_spor) %>% table()
  
  #     1     2 
  # 40891    41
  # 41 scaffolds have conflicting host sporulation predictions

  # which scaffolds?
  conflict_id <- 
    gpd_hosts %>% 
    group_by(GPD_id, f_spor) %>% 
    summarise() %>% 
    ungroup() %>%
    filter(duplicated(GPD_id)) %>% 
    pull(GPD_id) %>% 
    unique()
  
  gpd_hosts %>% 
    filter(GPD_id %in% conflict_id) %>% 
    mutate(GPD_id = fct_infreq(GPD_id)) %>% 
    group_by(GPD_id, f_spor) %>% 
    summarise(n=n()) %>% 
    arrange(n) %>% 
    ggplot(aes(GPD_id, n)) + 
    geom_col(aes(fill = f_spor)) + 
    coord_flip()
  
  # all have at least one predicted host that is likely a sporulator
  # so they have a potential to infect a sporulating cell
  # therefore calling them all phags of sporulators
  
  # finalize host sporulatin predictions for GPD ids
  gpd_hosts <- 
    gpd_hosts %>%
    select(GPD_id, f_spor) %>% 
    filter(!duplicated(GPD_id)) %>% 
    mutate(f_spor = if_else(GPD_id %in% conflict_id, TRUE, f_spor))
  
  #plot
  gpd_hosts %>%
    group_by(f_spor) %>%
    summarise(n=n()) %>%
    ggplot(aes(f_spor,n)) +
    geom_col(aes(fill = f_spor), color = "black")+
    theme_classic()+
    scale_fill_viridis_d(na.value = "grey50")+
    scale_color_viridis_d(na.value = "grey50")
  

# Export GPD host sporulation predictions ---------------------------------


  write_csv(gpd_hosts,
            here("dram_metaG/data/enrichment/gpd_spor_predictions.csv"))  
  