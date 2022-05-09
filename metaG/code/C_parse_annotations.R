library(here)
library(tidyverse)
library(foreach)

#### read in data ------------------------------
enriched <- read_csv(here("metaG/data/","gvd_enriched_fullData.csv"))
scaffolds_2get <- unique(enriched$scaffold)



# filtering function from chunked reading
f <- function(.data, pos) {
  filter(.data, scaffold %in% scaffolds_2get) 
}
# read annotations
y <- read_tsv_chunked(here("metaG/data/gvd/annotations.tsv"),
                 DataFrameCallback$new(f), chunk_size = 10000)


#### write annotations by scaffold ----------------

# make directory
if(!dir.exists(here("metaG/data/gvd/annotations"))){
  dir.create(here("metaG/data/gvd/annotations"))
}

foreach (scaf = scaffolds_2get, .packages = "readr") %do%
  {
    write_csv(y[y$scaffold==scaf,],
              here("metaG/data/gvd/annotations", paste0(scaf,".csv")))
  }


rm(y)
