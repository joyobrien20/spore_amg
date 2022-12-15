library(here)
library(tidyverse)
library(Biostrings)
library(seqinr)
library(cowplot)
# BiocManager::install("ggmsa")
library(ggmsa)
cDir <- "genes_of_interest/spo0A_cog"

aa_msa <- readAAMultipleAlignment(
  filepath = here(cDir,"/data/align-trim-tree/maxcc.afa"),
  format = "fasta")



# Plot all ----------------------------------------------------------------

d.aa_msa <- tidy_msa(aa_msa) %>%
  mutate(g = str_remove(name, "_.*")) %>%
  mutate(g = if_else(g == "cog", "bacteria", g))



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

# ggsave(filename = here(cDir, "plots", "align_coverage.png"), 
#        plot = p.gaps, width = 8, height = 3)


# Focused Plot ------------------------------------------------------------
# based on column gap summary I will trim the right-hand-side domain
# where viral genes are aligning

#zoom in to determine trimming limits
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
  geom_msa(data = d.aa_msa %>% 
             filter(position >trim.start &
                      position<trim.end),
           seq_name = T, font = NULL, border = NA,color = 'Zappo_AA')+
  facet_grid (fct_rev(g)~ ., scales = "free_y", space = "free_y")+
  theme_classic()+
  scale_x_continuous(breaks = seq(trim.start,trim.end, 20))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())

p1

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


# redo coverage plot
# sum columns -------------------------------------------------------------

p.gaps <- d.aa_msa %>% 
  # remove gappy virome sequences
  filter(!name %in% to_exclude) %>% 
  group_by(g,position) %>% 
  summarise(gaps = sum(character == "-")/n(),
            non_gaps = 1-gaps) %>% 
  ggplot(aes(position, non_gaps))+
  # geom_line(aes(color = g), size = 0.3, show.legend = F)+
  # geom_area(aes(fill = g), alpha = 0.7, show.legend = F)+
  geom_line( size = 0.5)+
  geom_area( fill = "grey70")+
  theme_classic()+
  facet_grid(fct_rev(g)~.)+
  panel_border(color = "black", size=1)+
  scale_x_continuous(breaks = c(1,seq(50,450, 50)))+
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0,0.5,1),limits = c(0,1.05))+
  coord_cartesian(expand = FALSE)+
  ylab("alignment coverage")+
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "bold"))
  
# p.gaps  

ggsave(filename = here(cDir, "plots", "align_coverage.png"), 
       plot = p.gaps, width = 8, height = 4)


# Focus on model bacteria ---------------------------------------------------------

model_bacteria <- c(
  # Bacillus subtilis 168
  "cog_NP_390302.1",
  # Clostridioides_difficile_630
  "cog_YP_001087707.1"
)


d.plot <- 
  d.aa_msa %>% 
  # trim
  filter(position %in% seq(trim.start,trim.end)) %>% 
  # remove gappy virome sequences
  filter(!name %in% to_exclude) %>% 
  #remove non model bacteria
  filter(name %in% model_bacteria | str_detect(g, "virome")) %>% 
  droplevels() %>% 
  mutate(name = fct_relevel(name, model_bacteria, after = 0) )


p1 <- ggplot() +
  geom_msa(data =  d.plot  ,
           seq_name =T, border = NA, color = "Zappo_AA")+
  scale_x_continuous(breaks = seq(250,400, 50))+
  theme_classic()+
  coord_cartesian(expand = FALSE)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())
#p1
ggsave(filename = here(cDir, "plots", "msa_focused.png"),
       plot = p1, width = 11.5, height = 4.2)


# combined plot -------------------------

p.bottom = plot_grid(NULL,p1, nrow = 1, rel_widths = c(1,9))
p.top = plot_grid(NULL,p.gaps, nrow = 1, rel_widths = c(1,20))
ggsave(filename = here(cDir, "plots", "msa_both.png"), 
       plot = plot_grid(p.top, NULL ,p.bottom, ncol = 1, rel_heights = c(1,0.2,2),
                        labels = c("(a)","","(b)")),
       width = 8, height = 6, bg = "white", dpi = 600)

# numbers --------------------------------
d.aa_msa %>%
  # remove gappy virome sequences
  filter(!name %in% to_exclude) %>% 
  select(name, g) %>% 
  distinct() %>% 
  group_by(g) %>% 
  summarise(n = n())
  

