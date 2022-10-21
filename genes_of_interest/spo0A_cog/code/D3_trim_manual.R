library(here)
library(tidyverse)
library(Biostrings)
library(seqinr)
# BiocManager::install("ggmsa")
library(ggmsa)
cDir <- "genes_of_interest/spo0A_cog"

aa_msa <- readAAMultipleAlignment(
  filepath = here(cDir,"/data/align-trim-tree/maxcc.afa"),
  format = "fasta")



# Plot all ----------------------------------------------------------------

d.aa_msa <- tidy_msa(aa_msa) %>%
  mutate(g = str_remove(name, "_.*")) %>%
  mutate(g = if_else(g == "cog", "Bacteria", g))



p <- ggplot() +
  geom_msa(data = d.aa_msa,# %>% filter(str_detect(g, "virome")) ,
           seq_name = T, font = NULL, border = NA, color = 'Zappo_AA')+
  facet_grid (fct_rev(g)~ ., scales = "free_y", space = "free_y")+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

ggsave(filename = here(cDir, "plots", "msa1.png"), plot = p, width = 8, height = 4)



# sum columns -------------------------------------------------------------

p.gaps <- d.aa_msa %>% 
  group_by(g,position) %>% 
  summarise(gaps = sum(character == "-")/n(),
            non_gaps = 1-gaps) %>% 
  ggplot(aes(position, non_gaps))+
  geom_line(aes(color = g), size = 0.3, show.legend = F)+
  geom_area(aes(fill = g), alpha = 0.7, show.legend = F)+
  theme_classic()+
  facet_grid(fct_rev(g)~.)+
  scale_x_continuous(breaks = seq(0,450, 50))

ggsave(filename = here(cDir, "plots", "align_coverage.png"), 
       plot = p.gaps, width = 8, height = 3)


# Focused Plot ------------------------------------------------------------
# based on column gap summary I will trim the right-hand-side domain
# where viral genes are aligning

#zoom in to etermine trimmin limits
zoom.start <- 220
zoom.end <- 380

p1 <-
  ggplot()+
  geom_msa(data = d.aa_msa %>% 
             filter(position >zoom.start &
                      position<zoom.end),
           seq_name = T, font = NULL, border = NA,color = 'Zappo_AA')+
  facet_grid (fct_rev(g)~ ., scales = "free_y", space = "free_y")+
  theme_classic()+
  scale_x_continuous(breaks = seq(zoom.start,zoom.end, 10))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())
p1

# based on zoom-in plot I will trim at these positons
trim.start <- 230
trim.end <- 375
p1 <-
  ggplot()+
  geom_msa(data = d.aa_msa,# %>% 
             # filter(position >pos.start &
             #          position<pos.end),
           seq_name = T, font = NULL, border = NA,color = 'Zappo_AA')+
  geom_vline(xintercept = c(trim.start,trim.end), size = 0.5)+
  facet_grid (fct_rev(g)~ ., scales = "free_y", space = "free_y")+
  theme_classic()+
  scale_x_continuous(breaks = c(seq(0,400, 100),trim.start,trim.end))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
# p1



ggsave(filename = here(cDir, "plots", "msa_trim.png"), 
       plot = p1, width = 8, height = 4)


# remove short virome sequences -------------------------------------------

# There are three virome sequences that are missing the left side of the domain
# I want to exclude these from phylogeny

to_exclude <- 
  d.aa_msa %>% 
  filter(g == "virome") %>% 
  filter(position %in% 250:270) %>% 
  group_by(name) %>% 
  summarise(gaps = sum(character == "-")/n()) %>% 
  arrange(desc(gaps)) %>% 
  filter(gaps>0.5) %>% 
  pull(name) %>% 
  as.character()



# Trim and remove ---------------------------------------------------------


trimmed_aa <- 
  d.aa_msa %>% 
  # trim
  filter(position %in% seq(trim.start,trim.end)) %>% 
  # remove gappy virome sequences
  filter(!name %in% to_exclude) %>% 
  # collapse back to aligned sequence
  arrange(position) %>% 
  group_by(g, name) %>% 
  summarise(seq = str_c(character, collapse = ""))


# export ------------------------------------------------------------------



write.fasta(names = trimmed_aa$name,
            sequences = str_split(trimmed_aa$seq,pattern = ""),
            file.out = here(cDir, "data", "/align-trim-tree/maxcc.afa.trim"))


