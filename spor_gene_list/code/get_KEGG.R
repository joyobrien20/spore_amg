# BiocManager::install("KEGGREST")
library(KEGGREST)
library(tidyverse)

keggs <-  keggLink("bsu", "ko")

          
keggs <- enframe(keggs, name = "ko", value = "bsu")

keggs$bsu <-  str_replace(keggs$bsu,pattern = "bsu:BSU",replacement = "BSU_")
keggs$ko <-  str_replace(keggs$ko,pattern = "ko:",replacement = "")
write_csv(keggs,"kegg_bsu.csv")

dspore[!(dspore$gene%in%keggs$bsu),]
