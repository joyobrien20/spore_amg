# get NCBI taxid matched with nucleotide accesions
library(here)
library(tidyverse)


#Downlad the mapping file from NCI FTP
download.file(url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
              destfile = here("data", "nucl_wgs.accession2taxid.gz"))

download.file(url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz.md5",
              destfile = here("data", "nucl_wgs.accession2taxid.gz.md5"))

