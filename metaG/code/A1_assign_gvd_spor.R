library(here)
library(tidyverse)

# downloaded gvd scaffold meta with predicted host taxonomy
# https://www.cell.com/cms/10.1016/j.chom.2020.08.003/attachment/f714d07e-7b44-4710-a4ca-2cdf1c027082/mmc3.xlsx
# deleted frst line and saved as csv(direct download from R got bad file)   

d.gvd <- read_csv(here("metaG/data/gvd/gvd_all_hosts.csv"))
d.f_spor <- read_csv(here("gtdb_spor/data/gtdb_families_sporulation.csv"))


# make gvd-like taxonomy
d.f_spor <- d.f_spor %>% 
  mutate(GTDB_Host = str_c(str_remove(gtdb_d, "d__"), gtdb_p,gtdb_f, sep = ";"))

#clean gvd
d.gvd <- d.gvd %>% 
  filter(Eukaryotic_or_Prokaryotic_Virus == "Bacteriophage") %>% 
  filter(!str_detect(GTDB_Host, "Host Not Assigned" )) %>% 
  filter(str_detect(GTDB_Host,"p__Firmicutes"))

# compatability check
d.gvd %>% 
  filter(!str_detect(GTDB_Host, "Host Not Assigned" )) %>% 
  pull(GTDB_Host)
           
table(d.gvd$GTDB_Host %in% d.f_spor$GTDB_Host, useNA = "a")
# FALSE  TRUE  <NA> 
#   301  6623     0 

d.gvd %>% 
  filter(!d.gvd$GTDB_Host %in% d.f_spor$GTDB_Host) %>% 
  pull(GTDB_Host) %>% 
  unique() %>% 
  length() #18


 # the other way
table(d.f_spor$GTDB_Host %in% d.gvd$GTDB_Host,useNA = "a")
# FALSE  TRUE  <NA> 
#   336    66     0

d.f_spor %>% 
  filter(d.f_spor$GTDB_Host %in% d.gvd$GTDB_Host) %>% 
  pull(f_spor) %>% 
  table(useNA = "a")
# FALSE  TRUE  <NA> 
#   40    19     7 


# There are 18 Firmicutes families missing from my list,
# accounting for 301 phages.
# my list accounts for 66 Firmicutes families in GVD 
# for which I have sporulation prediction for all but 7.

#reload gvd and join
d.gvd <- read_csv(here("metaG/data/gvd/gvd_all_hosts.csv"))

d.gvd <- d.f_spor %>% 
  select(GTDB_Host, f_spor) %>% 
  left_join(d.gvd, ., by = "GTDB_Host" )



#
  d.gvd %>% 
  filter(Eukaryotic_or_Prokaryotic_Virus == "Bacteriophage") %>% 
  filter(!str_detect(GTDB_Host, "Host Not Assigned" )) %>% 
    mutate(f_spor = if_else(
      str_detect(GTDB_Host,regex("p__Firmicutes", ignore_case = T)),
      f_spor, FALSE)) %>% 
    group_by(f_spor) %>% 
    summarise(n=n())
  # f_spor     n
  # <lgl>  <int>
  
  # 1 FALSE  10571
  # 2 TRUE    3036
  # 3 NA       326

#save gvd with sporulation prediction 
gvd_spor <- d.gvd; rm(d.gvd) 
save(gvd_spor,file = here("metaG/data/gvd/gvd_spor.Rdata"))

# delete gvd download 
unlink(here("metaG/data/gvd/gvd_all_hosts.csv"))
