#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
#
# Get AA sequences for a COG ----
#
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

# load packaes ------------------------------------------------------------
library(here)
library(tidyverse)
library(rentrez)
library(XML)
library(seqinr)

# get data for genes ----------------------------------------------------------------
cog <- "COG5830" # spoIID
cDir <- "genes_of_interest/spoIIID_cog"

# download fasta
download.file(paste0("https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fasta/",cog,".fa.gz"),
              destfile = here(cDir, "data", "cog.fa.gz"))

#read in fasta
faa <- read.fasta(here(cDir, "data", "cog.fa.gz"),as.string = F, seqtype = "AA")

# downlaod metadata
download.file(paste0("https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fasta/",cog,".tsv.gz"),
              destfile = here(cDir, "data", "cog.tsv.gz"))

# read in metadata
d <- read_tsv(here(cDir, "data", "cog.tsv.gz")) %>% 
  rename(taxid = axid)
  

# get taxonomy ------------------------------------------------------------

#initialize empty tibble
d.tax <- tibble()
tax_ids <- unique(d$taxid)

while(length(tax_ids)>0){
  for(i in tax_ids){
    
    #to catch errors (https://stackoverflow.com/a/59462794)
    skip_to_next <- FALSE
    tryCatch(
      # fetch taxonomy data from NCBI
      tax_rec <- entrez_fetch(db="taxonomy", id=i, 
                              rettype="xml", parsed=TRUE),
      error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next } 
    # if no error remove taxid from vector
    tax_ids <- tax_ids[-which(tax_ids==i)]
    
      # parse NCBI taxonomy
      #https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
      tax_list <- XML::xmlToList(tax_rec)
      tax <- do.call(rbind, tax_list$Taxon$LineageEx)
      
      # convert lineage to tibble
      tmp.tax <- 
        unlist(tax) %>%
        matrix(nrow = nrow(tax), ncol = ncol(tax)) %>% 
        as_tibble(.name_repair = "unique") %>% 
        set_names(dimnames(tax)[[2]]) 
      
      # wide format and bind to other listings
      d.tax <- 
        bind_rows(deframe(tmp.tax[,3:2])) %>% 
        # lots of different thigs are arke "no rank"
        # get rid of them
        select(-contains("no rank")) %>% 
        mutate(taxid = i) %>% 
        bind_rows(d.tax,.)   
    
  }
}
# add taxonomy to metadata
d <- left_join(d, d.tax, by = "taxid")




# get sequence  -----------------------------------------------------------

# convert to tibble
d.faa <- tibble(protein_id = getName(faa),
                # annot = getAnnot(faa),
                seq = getSequence(faa,as.string = T) %>% unlist())
# add sequence to metadata
d <- left_join(d, d.faa, by = "protein_id")


# export data -------------------------------------------------------------
write_csv(d, here(cDir, "data", "cog_taxa_seq.csv"))
  