# Make plots and tables of scaffolds with putative AMGs
# for manual curation
library(here)
library(tidyverse)
library(doParallel)
registerDoParallel()
library(foreach)


# list of sporulation genes enriched in gvd/gpd ------------------------------
d.enriched <- read_csv(here("dram_metaG/data/", "enrichment","spor_enriched.csv"))

#load plotting function
source(here("dram_metaG/code/enrichment_metaG/Z_function_plot_amg_scaffolds.R"))

# directory with the DRAM raw output
data_dir <- here("dram_metaG/data/dram_output")
sets <-  list.dirs(data_dir, full.names = F, recursive = F)


# function to group data sets
choose_sets <-  function(sg, v.set){
  if(sg == "GVD") return("Gregory_gvd")
  if(sg == "GPD") return("Camarillo-Guerrero")
  if(sg == "CSUsets") return(v.set[!v.set %in% c("Gregory_gvd","Camarillo-Guerrero")])
  
  
}
# start loop ------------------------------------------------
#For each putative AMG (by KO)
# extract scaffold annotation data and make plots
set_groups = c("CSUsets", "GVD", "GPD")
# loop over gene table
foreach (i = 1:nrow(d.enriched), .packages = "tidyverse") %:%
  foreach(j = seq(set_groups)) %dopar% {
  # folders of data to plot
  sets_2_plot <- choose_sets(set_groups[j], sets)
  
  # short name of gene to write on plots
  gene_name <-
    case_when(
      !is.na(d.enriched$bs[i]) ~
        str_remove(d.enriched$bs[i], " \\[.*") %>% paste0("Bs_", .),
      !is.na(d.enriched$cd[i]) ~
        str_remove(d.enriched$cd[i], " \\[.*")%>% paste0("Cd_", .),
      TRUE ~ paste0("NA_", i))
  
  plot_amg_scaffolds(data_dir = data_dir,
                     sets = sets_2_plot,
                     sets_name = set_groups[j],
                     amg_name = gene_name,
                     ko_amg = d.enriched$gene_id[i])


}
