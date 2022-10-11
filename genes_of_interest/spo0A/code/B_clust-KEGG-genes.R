library(here)
library(tidyverse)
library(seqinr)
library(foreach)

# load KEGG genes table
d <-read_csv(here("genes_of_interest","spo0A","data","spo0A_KEGG_sequences.csv"))


# cluster with cd-clust at different levels ---------------------------------------------
data_dir <- here("genes_of_interest","spo0A","data/")
if(!dir.exists(data_dir)){
  dir.create(data_dir)
}

# > make multi fasta ------------------------

faa_file = paste0(data_dir, "all_KEGG.faa")

d %>% 
  # prepare as fasta
  mutate(fasta = paste0(">",kegg_gene,"\n", seqAA)) %>% 
  pull(fasta) %>% 
  # write to file
  write_lines(., file = faa_file)
  
cdhit.cmd <- "module load cd-hit; cd-hit"

# cd-hit parameters 

# -aS	alignment coverage for the shorter sequence, default 0.0
# if set to 0.9, the alignment must covers 90% of the sequence

# c	sequence identity threshold, default 0.9
# this is the default cd-hit's "global sequence identity" calculated as:
#  	number of identical amino acids in alignment
#  	divided by the full length of the shorter sequence

v_aS = 0.9 #seq(0.5,1,0.05) # aS does not matter here
v_c = seq(0.4,1,0.05)

d.clusts <- tibble(aS = rep(v_aS, length(v_c)),
                   C = rep(v_c, each= length(v_aS)),
                   n.clusters = NA)

foreach(i = seq(1:nrow(d.clusts))) %do%
  {
    
    cmd <- 
      paste(cdhit.cmd, 
            "-aS", d.clusts$aS[i],
            "-c", d.clusts$C[i],
            "-d 0 -g 1 -sc 1", 
            "-i", faa_file, 
            "-o", paste0(data_dir,"clust_aS",d.clusts$aS[i],"_c",d.clusts$C[i])
            
      )
    
    cd.x <- system(cmd, intern = T)
    
    # extract number of clusters
    d.clusts$n.clusters[i] <- 
      cd.x %>% str_detect("clusters") %>% cd.x[.] %>% 
      str_remove(.,".*finished *") %>% parse_number()
    
    # -aS  alignment coverage for the shorter sequence, default 0.0
    #      if set to 0.9, the alignment must covers 90% of the sequence
    # -c   sequence identity threshold, default 0.9
    #      this is the default cd-hit's "global sequence identity" calculated as:
    #      number of identical amino acids in alignment
    #      divided by the full length of the shorter sequence
    # -d	length of description in .clstr file, default 20
    # 	if set to 0, it takes the fasta defline and stops at first space
    # -g	1 or 0, default 0
    #      by cd-hit's default algorithm, a sequence is clustered to the first 
    #      cluster that meet the threshold (fast cluster). If set to 1, the program
    #      will cluster it into the most similar cluster that meet the threshold
    #      (accurate but slow mode)
    #      but either 1 or 0 won't change the representatives of final clusters     
    # -sc	sort clusters by size (number of sequences), default 0, output clusters by decreasing length
    # 	if set to 1, output clusters by decreasing size
    
  }



d.clusts %>% 
  filter(!is.na(n.clusters)) %>% 
  mutate(aS = as_factor(aS)) %>% 
  ggplot(aes(x = C, y = n.clusters, color = aS))+
  geom_line()+
  geom_point(size=2, shape=21, fill = "white")
  # geom_text(aes(label = aS), nudge_y = 10, size =3 )

################
#
# Using C = 0.65
#
################

# parse cd-hit-est clstr ----
clstr <- read_lines(paste0(data_dir,"clust_aS",v_aS,"_c",0.65,".clstr"))

# initialize
d.clstr <- tibble()
n.clstr = 0

for(i in 1:length(clstr)){
  
  # assign clstr number
  n.clstr <- if_else(str_detect(clstr[i],"^>Cluster",), 
                     parse_number(clstr[i]),n.clstr)
  #row with clstr number has no more data
  if(str_detect(clstr[i],"^>Cluster")) next
  
  n.member <- str_remove(clstr[i], "\t.*")
  acc <- str_remove(clstr[i],".*, >") %>% str_remove("\\.\\.\\..*")
  aa <- str_extract(clstr[i],"\t.*aa") %>% parse_number()
  is.rep <- str_detect(clstr[i], "\\*$")
  perc.id <- if_else(is.rep, 100,
                     str_extract(clstr[i], "at.*%") %>% parse_number())
  
  
  d.clstr <- tibble(protClstr = n.clstr,
                    protein = acc,
                    aa = aa,
                    perc.id_protClstr = perc.id ,
                    protClstr_rep = is.rep) %>% 
    bind_rows(d.clstr, .)
  
}

# add cluster data to fasta table
d <- left_join(d, d.clstr, by = c("kegg_gene" = "protein") )

# replace representatives if needed -----------

# replacment function
swap_cluster_rep <- function(d, organism, cluster){
  current_rep <- which(d$protClstr == cluster &
                         d$protClstr_rep)
  new_rep <- which(str_detect(d$organism, organism) &
                     d$protClstr == cluster)
  d$protClstr_rep[current_rep] <- FALSE
  d$protClstr_rep[new_rep] <- TRUE
  return(d)
}

# > B. subtilis 168 -------
Bs168.clst <- d %>% 
  filter(str_detect(organism, "168")) %>% pull(protClstr)
d %>% 
  filter(protClstr %in% Bs168.clst) %>% view
#large cluster of 242 sequences, 168 is 81.24% simmilar to cluster rep.

d <- swap_cluster_rep(d, "168", Bs168.clst)
#check
d %>% 
  filter(protClstr %in% Bs168.clst) %>% view

# > C. difficile 630Derm (lab strain) -----------
# There are two 630 strains beside Derm, all in same cluster
# Clostridioides difficile 630Derm is already cluster rep
# I will "swap" it for reproducibility (in case he is not cluster rep in another run.)
Cd630.clst <- d %>% 
  filter(str_detect(organism, "630Derm")) %>% pull(protClstr)
d %>% 
  filter(protClstr %in% Cd630.clst) %>% view

d <- swap_cluster_rep(d, "630Derm", Cd630.clst)
#check
d %>% 
  filter(protClstr %in% Cd630.clst) %>% view

# Save results
write_csv(d, file = paste0(data_dir,"faa_KEGG_clustered.csv"))

# delete clustering files
to_delete <- list.files(data_dir, pattern = "clust_aS", full.names = T)
unlink(to_delete)
