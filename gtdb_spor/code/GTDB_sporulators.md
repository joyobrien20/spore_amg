GTDB sporulators
================
Daniel Schwartz
January/2022

# The Goal

We would like to test if sporulation genes found in phages are enriched
among phages that infect hosts that can sporulate. To this end we need
to define by taxonomy who are the sporulators among the *Firmicutes*.
The taxonomic level we would like to use fro prediction is that of
Family. The taxonomic scheme used in the metagenomic analysis is that of
the [Genome Taxonomy Database (GTDB)](https://gtdb.ecogenomic.org/).

## How do we identify a sporulator?

The most obvious way is to get experimental evidence from culture
analysis through microscopy and/or the presence of cells resistant to
heat, chemicals, etc. [Galperin
2013](https://doi.org/10.1128/microbiolspectrum.TBS-0015-2012) gives a
table with a summary of such data across the *Firmicutes*. See data
below. However this data is limited.

A more generalizable approach is the prediction of sporulation capacity
from genome-wide analysis, with a focus on the presence of predictive
sporulation genes. Such an approach was implemented by [Weller and Wu
2015](https://doi.org/10.1111/evo.12597) and on a more broad scale by
[Browne et al. 2021](https://doi.org/10.1186/s13059-021-02428-6).

The issue with both of these approaches with respect to our goal, is
that sporulation appears to have been lost multiple times within the
*Firmicutes*. This result is that within any taxa (in our case *family*)
of sporulators there are commonly also some non-sporulating members.
Therefore we will try and establish some prediction at the family level
of the likelihood of sporulation (likelihood used here to mean best
guess, not in its strict statstical meaning).

# Firmicutes in GTDB

Files will be obtained from the [latest version of
GTDB](https://data.gtdb.ecogenomic.org/releases/latest/).

We will use the genomes listed in the GTDB as Firmicurtes and their
taxonomy. The file which contains this data is the GTDB’s bacterial
metadata file. This is described in **FILE DESCRIPTIONS**:

> Metadata for all bacterial genomes including GTDB and NCBI taxonomies,
> completeness and contamination estimates, assembly statistics, and
> genomic properties.

To download the file I use the *Windows Linux Subsystem* (WSL). This
data cannot be synced by github so it needs to be downloaded at every
machine.

``` r
# make directory
setwd(here())
system(paste("wsl mkdir -p", "gtdb_spor/data/gtdb_downloads"))
setwd(here("gtdb_spor/data/gtdb_downloads"))

# file source
url_latest <- "https://data.gtdb.ecogenomic.org/releases/latest/"

# version
system(paste0("wsl wget ", url_latest,"VERSION"))


# download and extract
system(paste0("wsl wget ", url_latest,"bac120_metadata.tar.gz"))
f_meta <- system(paste0("wsl tar -xzvf ","bac120_metadata.tar.gz"), intern = T)

#delete original file
unlink(here("gtdb_spor/data/gtdb_downloads","bac120_metadata.tar.gz"))

setwd(here())
```

**GTDB version used is v202 Released April 27, 2021.**

### Firmicutes in the GTDB

``` r
f_meta <- list.files(here("gtdb_spor/data/gtdb_downloads"), pattern = "bac120.*tsv")
# read in only Firmicutes

  # function to filter firmicutes by chunks
  f <- function(x, pos) {
    x %>% filter(str_detect(gtdb_taxonomy,
                            regex("Firmicutes", ignore_case = T))) # %>% 
      # select(taxon_cols) 
      
  }
  
  d_meta <-
    read_tsv_chunked(here("gtdb_spor/data/gtdb_downloads", f_meta),
                     DataFrameCallback$new(f))
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   .default = col_character(),
    ##   ambiguous_bases = col_double(),
    ##   checkm_completeness = col_double(),
    ##   checkm_contamination = col_double(),
    ##   checkm_marker_count = col_double(),
    ##   checkm_marker_set_count = col_double(),
    ##   checkm_strain_heterogeneity = col_double(),
    ##   coding_bases = col_double(),
    ##   coding_density = col_double(),
    ##   contig_count = col_double(),
    ##   gc_count = col_double(),
    ##   gc_percentage = col_double(),
    ##   genome_size = col_double(),
    ##   gtdb_representative = col_logical(),
    ##   gtdb_type_species_of_genus = col_logical(),
    ##   l50_contigs = col_double(),
    ##   l50_scaffolds = col_double(),
    ##   longest_contig = col_double(),
    ##   longest_scaffold = col_double(),
    ##   lsu_23s_count = col_double(),
    ##   lsu_5s_count = col_double()
    ##   # ... with 24 more columns
    ## )
    ## i Use `spec()` for the full column specifications.

``` r
# separate taxa levels by GTDB and NCBI for each genome  
d_meta_comp <- d_meta %>% 
    select(accession, gtdb_taxonomy, ncbi_taxonomy) %>% 
    separate(gtdb_taxonomy, sep = ";", into = paste0("gtdb_",c("d","p","c","o","f","g","s"))) %>% 
    separate(ncbi_taxonomy, sep = ";", into = paste0("ncbi_",c("d","p","c","o","f","g","s")))

# list of GTDB firmicutes families
  gtdb_fams <- d_meta_comp %>% 
    select(gtdb_d,gtdb_p,gtdb_c,gtdb_o,gtdb_f) %>% 
    group_by(gtdb_d,gtdb_p,gtdb_c,gtdb_o,gtdb_f) %>%
    summarise(.groups = "drop") 
  

p.gtdb <- gtdb_fams %>% 
    group_by(gtdb_d,gtdb_p,gtdb_c,gtdb_o) %>% 
    summarise(n_families=n(), .groups = "drop") %>% 
  mutate(gtdb_o = fct_reorder(gtdb_o, desc(n_families))) %>% 
    ggplot(aes(gtdb_o, n_families))+
    geom_bar(stat="identity")+
    theme_classic()+
    theme(axis.text.x = element_blank())+
  xlab("Order")
  

# list of NCBI firmicutes families
  ncbi_fams <- d_meta_comp %>% 
    select(ncbi_d,ncbi_p,ncbi_c,ncbi_o,ncbi_f) %>% 
    group_by(ncbi_d,ncbi_p,ncbi_c,ncbi_o,ncbi_f) %>%
    summarise(.groups = "drop") 
  

p.ncbi <- ncbi_fams %>% 
      group_by(ncbi_d,ncbi_p,ncbi_c,ncbi_o) %>% 
    summarise(n_families=n(), .groups = "drop") %>% 
  mutate(tax=paste0(ncbi_p,ncbi_c,ncbi_o)) %>% 
  mutate(tax = fct_reorder(tax, desc(n_families))) %>% 
    ggplot(aes(tax, n_families))+
    geom_bar(stat="identity")+
    theme_classic()+
    theme(axis.text.x = element_blank())+
  xlab("Order")

plot_grid(p.gtdb, p.ncbi, labels = c("GTDB", "NCBI"), ncol = 1, label_x = 0.5, label_y = 0.9)
```

![](GTDB_sporulators_files/figure-gfm/firmi%20families-1.png)<!-- -->

GTDB has much more taxa levels. Overall GTDB lists 0 families in the
Firmicutes.

## Sporulator status prediticion from Browne et al. 

Data from table S1 from paper that lists prediction of sporulation
ability fro &gt;1000 Firmicutes genomes. Predictions are based on the
distribution of 66 sporulation genes within families, with special
weight to spo0A (See Browne et al. methods).

``` r
d_browne <- read_csv(here("gtdb_spor/data/Browne_2021_tbl-S1.csv"))
```

    ## Warning: Missing column names filled in: 'X13' [13], 'X14' [14], 'X15' [15],
    ## 'X16' [16]

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   `tree order` = col_double(),
    ##   `major taxonomic family/ species name for non-Firmicutes` = col_character(),
    ##   `Firmicutes taxonomic order/ non-Firmicutes phylum` = col_character(),
    ##   Accession = col_character(),
    ##   `Sanger accession` = col_character(),
    ##   environment = col_character(),
    ##   `sporulation signature score %` = col_double(),
    ##   `presence spo0A` = col_double(),
    ##   `spore-former characterisation` = col_double(),
    ##   `genome size` = col_double(),
    ##   `ethanol resistance` = col_character(),
    ##   `Erysipelotrichaceae experiments` = col_character(),
    ##   X13 = col_logical(),
    ##   X14 = col_logical(),
    ##   X15 = col_logical(),
    ##   X16 = col_logical()
    ## )

``` r
#filter out non-firmicutes (have species name in column #2)
d_browne <- d_browne %>% 
  filter(! str_detect(`major taxonomic family/ species name for non-Firmicutes`,
                    " "))


  browne_fams <- d_browne %>% 
    select(order = `Firmicutes taxonomic order/ non-Firmicutes phylum`,
           family = `major taxonomic family/ species name for non-Firmicutes`,
           spore_score = `sporulation signature score %`,
           sporulator = `spore-former characterisation` )

  
  browne_fams %>% 
    mutate(sporulator = as.logical(sporulator)) %>% 
     mutate(family = fct_infreq(family) ) %>% 
    ggplot(aes(family,spore_score ))+
    geom_jitter(aes(fill = sporulator), width = 0.1, shape=21, alpha=0.3)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
```

![](GTDB_sporulators_files/figure-gfm/get%20Browne-1.png)<!-- -->

I will use the family majority to predict sporulation likelihood of a
family.

``` r
browne_fams <- browne_fams %>% 
  select(-spore_score) %>% 
  mutate(sporulator = if_else(sporulator==1, "sporulators","non")) %>% 
  group_by(order, family, sporulator) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = sporulator, values_from = n, values_fill = 0) %>% 
  mutate(total = sporulators + non) %>% # total wxcluding NAs
    mutate(perc_spor = 100*sporulators/total) %>% 
    mutate(fam_spor = perc_spor > 50)
```

    ## `summarise()` has grouped output by 'order', 'family'. You can override using the `.groups` argument.

``` r
browne_fams %>% 
  mutate(family = fct_infreq(family) ) %>%
  ggplot(aes(family, perc_spor+1))+
  geom_col(aes(fill = fam_spor))+
  # geom_point(aes(y=total))+
  facet_grid(~order, scales = "free_x", space = "free_x")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        legend.position = "bottom")
```

![](GTDB_sporulators_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Add Browne predictions to GTDB families

``` r
gtdb_fams <- 
  browne_fams %>%
    ungroup() %>% 
    select(family, browne_spor=fam_spor) %>% 
    mutate(family = paste0("f__", family)) %>% 
    left_join(gtdb_fams, ., by =  c("gtdb_f" = "family"))

browne_n_pred <- sum(!is.na(gtdb_fams$browne_spor)) %>% as.numeric()
```

Using Browne et al. data we can assign sporulation prediction to 75 of
402 Firmicutes families in GTDB.

## predictions by Galperin

Galperin lists a 4 level system to describe fraction of sporulators in a
family. The table footnote says this:

> The distribution of sporeformers among the experimentally
> characterized members of the respective family is indicated as
> follows:  
> +++, all (or nearly all) characterized members of the family produce
> spores;  
> ++, a significant fraction of species are sporeformers;  
> +, the family includes some sporeformers;  
> -, no known sporeformers in the family.
