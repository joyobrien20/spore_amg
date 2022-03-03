# plot scaffold genes
library(here)
library(tidyverse)
library(gggenes)
library(cowplot)
# source(here("dram_metaG/code/ggwrap.R"))

# input variables --------------------------------------------------------
data_dir <- "dram_metaG/data/dram_output"

# set <- "Anaerobic_Digestors"
# 
# scaffolds_2get <- c(
#   "VIRSorter_AD_CERC_C7A_k121_18704_flag_0_multi_51_5194_len_112620_gene_9_gene_106-4790-95673-cat_4",
#   "VIRSorter_AD_CERC_C7A_k121_267218_flag_0_multi_67_9926_len_109514-cat_2"
# )
# 
# ko_amg <- "K07699" #spo0A
ko_amg <- c("K06333", "K06334") #cotJ
cur_set <-  "Huttenhower"
# Sporulation genes -------------------------------------------------------

spor_genes <- read_csv(here("spor_gene_list/data/dram_spore_genes.csv")) %>% 
  # genes with KO
  filter(!is.na(gene_id.ko)) %>% 
  separate(gene_description, into = c("symbol", "description"), sep = ";") %>% 
  rename(ko = gene_id.ko)

spor_ko <- 
  unique(spor_genes$ko)


# extract annot with amg --------------------------------------------------

d.amg <- read_tsv(here(data_dir, cur_set, "amg_summary.tsv"))

d.amg <- d.amg %>% 
  filter(gene_id %in% ko_amg)

# # host taxonomy from A1_assign_gvd_spor.R
# load(file = here("metaG/data/gvd/gvd_spor.Rdata"))
# gvd_spor <- 
#   gvd_spor %>% filter(Contig %in% d.amg$scaffold) %>% 
#   filter(f_spor)

# get scaffold annotation -------------------------------------------------





  # filtering function from chunked reading
  f <- function(.data, pos) {
    filter(.data, scaffold %in% d.amg$scaffold) 
  }
  
  annot_file <- 
    list.files(here(data_dir, cur_set), pattern = "annotation", full.names = T)
  
  # read annotations
  d.annot <- read_tsv_chunked(annot_file,
                              DataFrameCallback$new(f), 
                              chunk_size = 10000) 



# add index and set
d.annot$idx <- 
  d.annot %>% 
  group_by(scaffold) %>% 
  group_indices()

d.annot$set <- cur_set



# add informative properties for AMG --------------------------------------

# kegg_hit & pfam_hits
  # Phages will typically have many genes with no annotation ("NA") 
  # or hypothetical annotations. Bacteria typically have well annotated stretches.
d.annot <- d.annot %>% 
  mutate(has_annot = !(is.na(kegg_hit) & is.na(pfam_hits)))


# Hallmark viral gene terms in annotation
v_hallmark <- 
  c( "virion", "capsid", "tail", "terminase", "Baseplate",
     "phage", "virus", "Reverse transcriptase", "head")

d.annot <- d.annot %>% 
  mutate(viral_hallmark = 
           str_detect(kegg_hit, str_c(v_hallmark, collapse = "|"))|
           str_detect(pfam_hits, str_c(v_hallmark, collapse = "|"))) %>% 
  mutate(viral_hallmark = if_else(is.na(viral_hallmark), FALSE, viral_hallmark)) %>% 
  # false positive
  mutate(viral_hallmark = 
           if_else(
             str_detect(kegg_hit, "Minor_tail_Z Laminin_I")|
               str_detect(pfam_hits, "Minor_tail_Z Laminin_I"),
             FALSE, viral_hallmark))

# hypothetical genes are a special class
d.annot <-
  d.annot %>% 
  mutate(one_NA = is.na(kegg_hit) | is.na(pfam_hits)) %>% 
  mutate(hypothetical = 
           case_when(
             #both hypothetical
             str_detect(kegg_hit, regex("hypothetical", ignore_case = T))&
               str_detect(pfam_hits, regex("hypothetical", ignore_case = T)) ~ TRUE,
             
             # one hypothetical an other NA
             str_detect(kegg_hit, regex("hypothetical", ignore_case = T))&
               one_NA ~ TRUE,
             str_detect(pfam_hits, regex("hypothetical", ignore_case = T))&
               one_NA ~ TRUE,
             TRUE ~ FALSE
           )
         )
  
  
           
###
# combine
d.annot <- d.annot %>% 
  mutate(gene_type = case_when(
    viral_hallmark ~ "viral",
    hypothetical ~ "hypothetical",
    has_annot ~ "other annotation",
    TRUE ~ "unannotated"
  ) %>% as_factor() %>% fct_relevel("other annotation"))


# plot --------------------------------------------------------------------
# https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html

virsort_temparate <- 
  d.annot %>% 
  filter(str_detect(scaffold, "cat_4|cat_5")) %>% 
  pull(idx) %>% unique()
        
virsort_lytic <- 
  d.annot %>% 
  filter(str_detect(scaffold, "cat_1|cat_2")) %>% 
  pull(idx) %>% unique()



n_genomes <- virsort_lytic[seq(1,22)] %>% sort()
n_genomes <- virsort_lytic[seq(33,48)] %>% sort()
n_genomes <- virsort_temparate[seq(1,22)] %>% sort()
n_genomes <- virsort_temparate[seq(22,32)] %>% sort()
# n_genomes <- seq(22, 22+21)
p <- d.annot %>% 
  filter(idx %in% n_genomes) %>% 
  ggplot(aes(y = strandedness/3,
             xmin = start_position, 
             xmax = end_position,
             fill = gene_type,
             forward = strandedness
  )) +
  geom_hline(yintercept = 0, color = "black", size=0.1)+
  geom_rect(data = d.annot %>% 
              filter(str_detect(kegg_hit, ko_amg) | (ko_id %in% ko_amg)) %>% 
              filter(idx %in% n_genomes),
            aes(xmin = start_position, xmax = end_position),
            ymin = -Inf, ymax = Inf, fill = "pink") +
  geom_gene_arrow(color = "black") +
  # geom_gene_arrow(data = filter(d.annot, 
  #                               str_detect(kegg_hit, spor_ko) |(kegg_id %in% spor_ko)),
  #                 fill = "purple") +
  geom_gene_arrow(data = d.annot %>% 
                    filter(str_detect(kegg_hit, ko_amg) | (ko_id %in% ko_amg)) %>% 
                    filter(idx %in% n_genomes),
                  aes(y=strandedness/3),
                  size=0.5, fill = "red") +
  geom_gene_label(aes(label = gene_position),
                  align = "centre") +
  facet_wrap(set ~ idx, scales = "free_x", strip.position = "left", ncol = 2,
             labeller = label_wrap_gen()) +
  theme_classic()+
  panel_border(color = "black", size=0.5)+
  theme(axis.text.y = element_blank(),
        legend.position = "bottom",
        # strip.background = element_blank(),
        strip.text = element_text(size = 5),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  scale_x_continuous(expand = c(0, 0))+
  scale_fill_viridis_d()+
  ylim(-1,1)+
  guides(fill = guide_legend(title = "gene\ncolor"))

ggsave(here("dram_metaG/plots", "cotJ_multi_genomes.png"),
       plot = p,
       width = 12, height = 8)

# write_csv(scaf_2plot, here("dram_metaG/data/scaffolds_index.csv"))
# write_csv(d.annot, here("dram_metaG/data/scaffolds_annotationa.csv"))

d.annot %>% 
  filter(idx == 56) %>% 
  select(scaffold,gene_position, strandedness, ko_id, kegg_hit, pfam_hits,
         hypothetical, viral_hallmark, has_annot) %>% view()
