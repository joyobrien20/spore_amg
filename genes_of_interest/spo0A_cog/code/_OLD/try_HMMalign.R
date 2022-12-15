library(here)
library(tidyverse)
library(seqinr)

cDir <- "genes_of_interest/spo0A_cog"
oDir <- here(cDir, "data/hmm")
hmms <- list.files(oDir, pattern = ".hmm|.HMM", full.names = F)
input <- ("virome_spo0A.faa")

setwd(oDir)


wsl <- paste("wsl hmmalign -o",
             "try_spo0A_hmm.aln",
             hmms[2],input)
system(wsl)
