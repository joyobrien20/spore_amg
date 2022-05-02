# test plotting function
library(here)

source(source(here("dram_metaG/code/plot_amg_scaffolds.R")))

cur_dir = "dram_metaG/data/dram_output"

# sets from CSU (excluding GVD and GPD)
sets_name <- "CSUsets"
sets <- list.dirs(here(cur_dir), full.names = F, recursive = F)
sets <- sets[!sets %in% "Camarillo-Guerrero"]
sets <- sets[!sets %in% "Gregory_gvd"]


plot_amg_scaffolds(data_dir = cur_dir,
                   sets = sets,
                   sets_name = "CSUsets",
                   amg_name = "ispF",
                   ko_amg = "K01770")



# Data for building the function -----------------
data_dir <- "dram_metaG/data/dram_output"




# AMG looking at currently
amg_name <- "cotJ"
ko_amg <- c("K06333", "K06334")