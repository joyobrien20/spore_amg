library(here)
library(tidyverse)

# download GTDB files -----------
setwd(here())
system(paste("wsl mkdir -p", "gtdb_spor/data/gtdb_downloads"))
setwd(here("gtdb_spor/data/gtdb_downloads"))

url_latest <- "https://data.gtdb.ecogenomic.org/releases/latest/"

# From FILE DESCRIPTIONS
# * bac120_metadata_<release>.tar.gz
# Metadata for all bacterial genomes including GTDB and NCBI taxonomies, completeness and
# contamination estimates, assembly statistics, and genomic properties.

system(paste0("wsl wget ", url_latest,"VERSION"))
system(paste0("wsl wget ", url_latest,"bac120_metadata.tar.gz"))
f_meta <- system(paste0("wsl tar -xzvf ","bac120_metadata.tar.gz"), intern = T)

# find columns with taxonmy
taxon_cols <- read_tsv(here("gtdb_spor/data/gtdb_downloads", f_meta), n_max = 1) %>% 
  colnames() %>% 
  grep("taxonomy", .)

# read in only Firmicutes
  
  # function to filter firmicutes by chunks
  f <- function(x, pos) {
    x %>% filter(str_detect(gtdb_taxonomy,
                            regex("Firmicutes", ignore_case = T))) # %>% 
      # select(taxon_cols) 
      
  }
  
  d_meta <-
    read_tsv_chunked(here("gtdb_spor/data/gtdb_downloads", f_meta),
                     DataFrameCallback$new(f))
  
  d_meta_comp <- d_meta %>% 
    select(accession, gtdb_taxonomy, ncbi_taxonomy) %>% 
    separate(gtdb_taxonomy, sep = ";", into = paste0("gtdb_",c("d","p","c","o","f","g","s"))) %>% 
    separate(ncbi_taxonomy, sep = ";", into = paste0("ncbi_",c("d","p","c","o","f","g","s")))
  
# list of GTDB firmicutes families
  gtdb_fams <- d_meta_comp %>% 
    select(gtdb_d,gtdb_p,gtdb_c,gtdb_o,gtdb_f) %>% 
    group_by(gtdb_d,gtdb_p,gtdb_c,gtdb_o,gtdb_f) %>%
    summarise(n=n())

# Classify by Galperin ----------------------------------------------------
  # classification of sporulators by Galperin 2013
  spore.fam <- read_csv(here("enrichment/data/Galperin_2013_MicrobiolSpectrum_table2.csv"))
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
    mutate(spore.likely = case_when(Fraction.spore.forming == "+++" ~ TRUE,
                                    Fraction.spore.forming == "++" ~ TRUE,
                                    Fraction.spore.forming == "+" ~ FALSE,
                                    Fraction.spore.forming == "-"  ~ FALSE)) 
  
  
  
  # how many can I classify using this list
  ncbi_fam <- 
    d_meta_comp %>% 
    pull(ncbi_f) %>% 
    str_remove("f__")
  
  (ncbi_fam %in% spore.fam$Family %>% sum())/length(ncbi_fam)
  # 92% of GTDB genomes (=species?)
  
  #add Galperin sporulation classification
  d_meta_comp <- d_meta_comp %>% 
    # adjust NCBI family
    mutate(ncbi_fam = str_remove(ncbi_f , "f__")) %>% 
    left_join(., spore.fam %>% select(Family, spore.likely),
              by = c("ncbi_fam" = "Family"))
  
  # # is the classification consistent in GTDB taxonomy?
  # d_meta_comp %>% 
  #   group_by(gtdb_f, spore.likely) %>% 
  #   summarise(n=n()) %>% 
  #   group_by(gtdb_f) %>% 
  #   summarise(count = n()) %>% 
  #   filter(count>1) %>% 
  #   view()
  
  p <- d_meta_comp %>% 
    filter(! gtdb_c %in% c("c__Bacilli", "c__Clostridia")) %>% 
    group_by(gtdb_c,gtdb_f, spore.likely) %>% 
    summarise(n=n()) %>% 
    group_by(gtdb_c,gtdb_f) %>% 
    mutate(rel_abun = n/sum(n)) %>% 
    ggplot(aes(gtdb_f, rel_abun))+
    geom_col(aes(fill = spore.likely))+
    facet_wrap(~gtdb_c, scales = "free")+
    theme_classic()+
    theme(axis.text = element_blank())

  ggsave(here("gtdb_spor/plots/","Galperin_prediction.png"), p,
         width = 11, height = 8 )  
  
  p <- d_meta_comp %>% 
    filter( gtdb_c %in% c("c__Bacilli", "c__Clostridia")) %>% 
    group_by(gtdb_c,gtdb_f, spore.likely) %>% 
    summarise(n=n()) %>% 
    group_by(gtdb_c,gtdb_f) %>% 
    mutate(rel_abun = n/sum(n)) %>% 
    ggplot(aes(gtdb_f, rel_abun))+
    geom_col(aes(fill = spore.likely))+
    facet_wrap(~gtdb_c, scales = "free", ncol = 1)+
    theme_classic()+
    theme(axis.text = element_blank())
  
  ggsave(here("gtdb_spor/plots/","Galperin_prediction_2.png"), p,
         width = 11, height = 8 )  
  
  # call by majority
  gtdb_fams <- 
    d_meta_comp %>% 
      group_by(gtdb_f, spore.likely) %>% 
      summarise(n = n() ) %>% 
      pivot_wider(names_from = "spore.likely", values_from = "n", values_fill = 0) %>% 
      mutate(total = `TRUE` + `FALSE`) %>% # total wxcluding NAs
      mutate(perc_spor = 100*`TRUE`/total) %>% 
      mutate(glp_spor = perc_spor > 50) %>% 
      select(gtdb_f, glp_spor) %>% 
      left_join(gtdb_fams, .)

# Browne classification ----------------------------------------------------
  # Browne, H.P., Almeida, A., Kumar, N. et al. 
  # Host adaptation in gut Firmicutes is associated with sporulation loss and 
  # altered transmission cycle. Genome Biol 22, 204 (2021). 
  # https://doi.org/10.1186/s13059-021-02428-6
  
setwd(here())
d_browne <- read_csv(here("gtdb_spor/data/Browne_2021_tbl-S1.csv"))

#filter out non-firmicutes (have species name in column #2)
d_browne <- d_browne %>% 
  filter(! str_detect(`major taxonomic family/ species name for non-Firmicutes`,
                    " "))


# matching Browne to GTDB by assembly
d_browne %>% 
  mutate(in_meta = 
           d_browne$Accession %>% str_remove(".*_") %>% str_remove("\\..") #remove version decimal
         %in% 
           (d_meta$accession %>% str_remove(".*_") %>% str_remove("\\..")))  %>% 
  group_by(in_meta,`spore-former characterisation`) %>% 
  summarise(n=n())
# not all there, but there are ~1000 shared genomes, 
# evenly split between sporulators and non=sporulators




# add Browne classification
d_meta_comp <- d_browne %>% 
  mutate(join_acc = Accession %>% str_remove(".*_") %>% str_remove("\\..")) %>% 
  select(join_acc ,brwn_spor_stat= `spore-former characterisation`) %>% 
  left_join(d_meta_comp %>% 
              mutate(join_acc = accession %>% str_remove(".*_") %>% 
                       str_remove("\\..")), .)

# call by majority
gtdb_fams <- 
  d_meta_comp %>% 
    mutate(brwn_spor_stat = as.logical(brwn_spor_stat)) %>% 
    group_by(gtdb_f, brwn_spor_stat) %>% 
    summarise(n = n() ) %>% 
    pivot_wider(names_from = "brwn_spor_stat", values_from = "n", values_fill = 0) %>% 
    mutate(total = `TRUE` + `FALSE`) %>% # total wxcluding NAs
    mutate(perc_spor = 100*`TRUE`/total) %>% 
    mutate(brwn_spor = perc_spor > 50) %>% 
    select(gtdb_f, brwn_spor) %>% 
    left_join(gtdb_fams, .) %>% 
    mutate(final_spor = case_when(
      brwn_spor == TRUE ~ TRUE,
      brwn_spor == FALSE ~ FALSE,
      glp_spor == TRUE ~ TRUE,
      glp_spor == FALSE ~ FALSE,
      TRUE ~ NA
    )) 

gtdb_fams %>% 
  group_by(final_spor) %>% 
  summarise(fams =n(), sn = sum(n))

# d_meta_comp %>% 
#   filter(! is.na(brwn_spor_stat)) %>% 
#   group_by(gtdb_c,gtdb_f, brwn_spor_stat) %>% 
#   summarise(n=n()) %>% 
#   filter(n>3) %>% 
#   group_by(gtdb_c,gtdb_f) %>% 
#   mutate(rel_abun = n/sum(n)) %>% 
#   ggplot(aes(gtdb_f, rel_abun))+
#   geom_col(aes(fill = as.character(brwn_spor_stat)))+
#   facet_wrap(~gtdb_c, scales = "free")+
#   theme_classic()+
#   theme(axis.text = element_blank(),
#         legend.position = "bottom")
  


  
  


  