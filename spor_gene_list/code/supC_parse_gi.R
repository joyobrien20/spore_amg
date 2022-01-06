library(here)
library(tidyverse)

# convert GI to more useful locus tags
# using eutils installed in  Windows Subsystem for Linux (wsl)

# Ramos-Silva's list of C.diff sporulation genes
d.ramos <- read_csv(here("spor_gene_list/data", "cdif_RamosSilva_2019_S1CD.csv" ))
gi <- parse_number(d.ramos$GI)

d <- tibble()

for (g in gi){
  if(is.na(g)) next
  wsl <- paste("wsl esearch -db protein -query ",g," | efetch -format gb")
  
  gb <- system(wsl, intern = T)
  gb2 <- gb[c(grep("locus", gb),grep("product", gb),grep("ACCESSION", gb))]
  
  d <- as_tibble(gsub("^ */","",x = gb2)) %>% 
    mutate(value = str_remove_all(value,"\"")) %>% 
    mutate(value =str_replace(value, "ACCESSION *", "acc=")) %>% 
    separate(value, into = c("name","val"), sep = "[=]") %>% 
    pivot_wider(names_from = "name", values_from = "val") %>% 
    mutate(gi = g) %>% 
    bind_rows(d,.)
  
  Sys.sleep(1)
}

write_csv(d, here("spor_gene_list/data/giLocus_RamosSilva_.csv"))


