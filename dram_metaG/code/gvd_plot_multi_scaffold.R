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
# Sporulation genes -------------------------------------------------------

spor_genes <- read_csv(here("spor_gene_list/data/dram_spore_genes.csv")) %>% 
  # genes with KO
  filter(!is.na(gene_id.ko)) %>% 
  separate(gene_description, into = c("symbol", "description"), sep = ";") %>% 
  rename(ko = gene_id.ko)

spor_ko <- 
  unique(spor_genes$ko)


# extract annot with amg --------------------------------------------------

d.amg <- read_tsv(here(data_dir, "extra_data/gvd", "annotations.tsv"))

d.amg <- d.amg %>% 
  filter(ko_id %in% ko_amg)

# host taxonomy from A1_assign_gvd_spor.R
load(file = here("metaG/data/gvd/gvd_spor.Rdata"))
gvd_spor <- 
  gvd_spor %>% filter(Contig %in% d.amg$scaffold) %>% 
  filter(f_spor)

# get scaffold annotation -------------------------------------------------





  # filtering function from chunked reading
  f <- function(.data, pos) {
    filter(.data, scaffold %in% gvd_spor$Contig) 
  }
  
  annot_file <- 
    list.files(here(data_dir, "extra_data/gvd"), pattern = "annotation", full.names = T)
  
  # read annotations
  d.annot <- read_tsv_chunked(annot_file,
                              DataFrameCallback$new(f), 
                              chunk_size = 10000) 



# add index and set
d.annot$idx <- 
  d.annot %>% 
  group_by(scaffold) %>% 
  group_indices()

d.annot$set <- "gvd"

# plot --------------------------------------------------------------------
# https://cran.r-project.org/web/packages/gggenes/vignettes/introduction-to-gggenes.html

# gggenes::example_genes
# molecule  gene  start    end  strand orientation
# Genome5  genA 405113 407035 forward          -1


# ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
#   geom_gene_arrow() +
#   facet_wrap(~ molecule, scales = "free", ncol = 1) +
#   scale_fill_brewer(palette = "Set3")
n_genomes <- 30
p <- d.annot %>% 
  filter(idx<=n_genomes) %>% 
  mutate(vir_hit = if_else(is.na(viral_hit), NA_character_, "viral hit")) %>% 
  ggplot(aes(y = strandedness/3,
             xmin = start_position, 
             xmax = end_position,
             fill = vir_hit,
             forward = strandedness
  )) +
  geom_hline(yintercept = 0, color = "black", size=0.1)+
  geom_gene_arrow(color = "black") +
  # geom_gene_arrow(data = filter(d.annot, 
  #                               str_detect(kegg_hit, spor_ko) |(kegg_id %in% spor_ko)),
  #                 fill = "purple") +
  geom_gene_arrow(data = d.annot %>% 
                    filter(str_detect(kegg_hit, ko_amg) | (ko_id %in% ko_amg)) %>% 
                    filter(idx<=n_genomes),
                  aes(y=strandedness/2, fill = ko_id),
                  size=0.5, color = "red") +
  geom_gene_label(aes(label = gene_position),
                  align = "centre") +
  facet_wrap(set ~ idx, scales = "free_x", strip.position = "left", ncol = 3,
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
  ylim(-1,1)+
  guides(fill = guide_legend(title = "gene\ncolor"))

ggsave(here("dram_metaG/plots", "gvd_cotJ_multi_genomes.png"),
       plot = p,
       width = 12, height = 8)

# write_csv(scaf_2plot, here("dram_metaG/data/scaffolds_index.csv"))
# write_csv(d.annot, here("dram_metaG/data/scaffolds_annotationa.csv"))

