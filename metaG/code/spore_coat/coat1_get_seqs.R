library(here)
library(tidyverse)
library(seqinr)


# pull spore coat proteins detected in GVD 

#### data not on GitHub
dram.annotations <- "/N/slate/danschw/GVD/annotations.tsv"
gvd <- "/N/slate/danschw/GVD/GVDv1_viralpopulations.fna"
gvd_headers <- "/N/slate/danschw/GVD/gvd_headers.txt"
  #made by grep:
  #grep -n ">" GVDv1_viralpopulations.fna > gvd_headers.txt

#### read in data ------------------------------
enriched <- read_csv(here("metaG/data/","gvd_enriched_fullData.csv"))

enriched <- enriched %>% 
  filter(gene_id=="K06333")
scaffolds_2get <- unique(enriched$scaffold)


# filtering function from chunked reading
f <- function(.data, pos) {
  filter(.data, ko_id=="K06333") 
}
# read annotations
d.annot <- read_tsv_chunked(dram.annotations, DataFrameCallback$new(f), chunk_size = 10000)


#extract and translate CDS
headers <- read_delim(gvd_headers, delim = ":", col_names = c("n.line", "scaffold")) %>% 
  mutate(scaffold = str_remove(scaffold, ">")) %>% 
  filter(scaffold %in% scaffolds_2get)

for(i in 1:nrow(d.annot)){
  cur_scaffold <- d.annot$scaffold[i]
  cur_line <- headers %>% filter(scaffold==cur_scaffold) %>% pull(n.line)
  dna <- read_lines(gvd, n_max = 1,
                    skip = cur_line) %>% as.SeqFastadna(name = cur_scaffold)
  
  getFrag(dna, begin = d.annot$start_position[i], end = d.annot$end_position[i]) %>% 
    translate(sens = if_else(d.annot$strandedness[i]==1, "F", "R")) %>% 
    # remove stop codon asterisk
    getFrag(begin = 1, end = getLength(.)-1) %>% 
    write.fasta(file.out = here("metaG/data/coat/cotJB_gvd.faa"), 
                names = d.annot$X1[i], open = if_else(i==1, "w","a"))
}



##### combine sequences -------

# sequences from 3 sources
# 1. GVD DRAM output (above)
seq_gvd <- here("metaG/data/coat","cotJB_gvd.faa")
# 2. BLASTp of B. subtilis cotJB with taxonomy limited to viruses
  # results from gut virome
seq_blast <- here("metaG/data/coat","cotJB_blastp_virus.faa")
# 3. bacteria sequences clusteres at 50% ID (uniref50).
seq_uniref <- here("metaG/data/coat","uniref-cotjb-filtered-identity 0.5.fasta")

fa.gvd <- read.fasta(seq_gvd)
fa.blast <- read.fasta(seq_blast)
fa.uniref <- read.fasta(seq_uniref)

table_fasta <- function(.seq, set){
  tibble(dataset= set,
         header = getName(.seq),
         seq = getSequence.character(.seq),
         len = getLength(.seq)) %>% return(.)
}

all.seq <- 
  bind_rows(table_fasta(fa.gvd,"gvd"), 
            table_fasta(fa.blast,"blastp"), 
            table_fasta(fa.uniref,"uniref"))



all.seq <- 
  all.seq %>% 
  mutate(vb = if_else(dataset=="uniref", "bacteria", "phage")) %>% 
  group_by(dataset) %>% 
  mutate(i=1, n=cumsum(i)) %>% 
  ungroup() %>% 
  mutate(new_header = str_c(vb, dataset, n, sep = "_"))


write.fasta(sequences = all.seq$seq, 
            names = all.seq$new_header, 
            file.out = here("metaG/data/coat","cotJB_to_align.faa"))

as.SeqFastaAA(all.seq$str_seq)


# save data
save(all.seq, file = here("metaG/data/coat","CoatSeqData.Rdata"))
#######__________________--------------
# # extract DNA with surrounding
# for(i in 1:nrow(d.annot)){
#   cur_scaffold <- d.annot$scaffold[i]
#   cur_line <- headers %>% filter(scaffold==cur_scaffold) %>% pull(n.line)
#   dna <- read_lines(gvd, n_max = 1,
#                     skip = cur_line) %>% as.SeqFastadna(name = cur_scaffold)
#   
#   getFrag(dna, begin = d.annot$start_position[i]-1000, end = d.annot$end_position[i]+1000) %>% 
#     write.fasta(file.out = here("metaG/data/coat/cotJB_region.fna"), 
#                 names = d.annot$X1[i], open = if_else(i==1, "w","a"))
# }





# # filtering function from chunked reading
# f <- function(.data, pos) {
#   filter(.data, scaffold %in% scaffolds_2get) 
# }
# # read annotations
# y <- read_tsv_chunked(dram.annotations,
#                       DataFrameCallback$new(f), chunk_size = 10000)
