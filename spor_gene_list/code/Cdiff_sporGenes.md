C. difficile sporulation genes
================
Daniel Schwartz
Dec/2021

The goal of this analysis is to generate a list of sporulation genes
from *Clostridioides difficle* to use in DRAM.

# Lists of C. diff sporulation genes

1.  [Fimlaid et al. 2013](https://doi.org/10.1371/journal.pgen.1003660)

> we generated loss-of-function mutations in genes encoding these
> sporulation sigma factors and performed RNA-Sequencing to identify
> specific sigma factor-dependent genes.

Analysis done in strain 630, genes listed with locus\_tag
(*CD630\_NNNN*).

``` r
d.fim <- read_csv(here("spor_gene_list/data", "cdif_Fimlaid_2013_S9.csv" ))
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   Name = col_character(),
    ##   locus_tag = col_character(),
    ##   description = col_character()
    ## )

2.  [Dembek et al. 2015](https://doi.org/10.1128/mBio.02383-14)

These authors made a transposon mutant library in C. difficile epidemic
strain R20291. They then grew library cells in sporulation media and
purified spores. Sporulation genes were those Tn-mutants that were
missing in the spores. In a similar way they also isentifies germination
genes as those absent in a culture grown from purified spores.

Analysis done in strain R20291, genes listed with locus\_tag
(*CDR20291\_NNNN*). They also list for each gene the ortholog in strain
630, if available. However, they do not mention how orthology was
determined.

``` r
d.dem <- read_csv(here("spor_gene_list/data", "cdif_Dembek_2015_S2.csv" ))
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   locus_tag = col_character(),
    ##   `function` = col_character(),
    ##   logFC = col_double(),
    ##   q.value = col_double(),
    ##   `Cd630 orthologue locus_tag` = col_character()
    ## )

3.  [Ramos-Silva et al. 2019](https://doi.org/10.1101/473793)

> All known sporulation genes from B. subtilis strain 168 and C.
> difficile strain 630, the two species that are better characterized in
> terms of their sporulation machinery, were collected from the
> literature…)

Genes of strain 630 are listed bu GI number. In a separate script
(parse\_gi.R) I converted the GI numbers to C. diff 630 locus tags. That
is the list I will use here.

``` r
# d.ramos <- read_csv(here("spor_gene_list/data", "cdif_RamosSilva_2019_S1CD.csv" ))
d.ramos <- read_csv(here("spor_gene_list/data", "giLocus_RamosSilva_.csv" )) %>% 
  # Remove columns created without need
  select(locus_tag, old_locus_tag, product, acc, gi)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   locus_tag = col_character(),
    ##   old_locus_tag = col_character(),
    ##   product = col_character(),
    ##   acc = col_character(),
    ##   gi = col_double(),
    ##   `substrate for methyltransfer, creating the product` = col_logical(),
    ##   note = col_character(),
    ##   `(DIM6/NTAB) family [Energy production and conversion];` = col_logical()
    ## )

4.  [Saujet et al. 2013](https://doi.org/10.1371/journal.pgen.1003756)

These autors constructes mutants in the major regulators of sporulation
(sigma factors EFGK and spoIIID) and compared gene expression between
each mutant abd a WT strain during sporulation (time selected by
preliminary test to maximize differential expression) using microarrays.

Analysis done in strain 630, genes listed with old locus\_taga
(*CDNNNN*). To match these I use a gene data table from [Petit et
al. 2014](https://doi.org/10.1186/1471-2164-15-160).

``` r
d.sau <- read_csv(here("spor_gene_list/data", "cdif_Saujet_2013_sup.csv"),
                   trim_ws = T)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   Gene = col_character(),
    ##   symbol = col_character(),
    ##   Function = col_character(),
    ##   FC = col_double(),
    ##   header = col_character(),
    ##   mutant.strain = col_character()
    ## )

``` r
d.pet <- read_csv(here("spor_gene_list/data", "cdif_Pettit_2014_spo0A.csv" ),
                  trim_ws = T)
```

    ## 
    ## -- Column specification --------------------------------------------------------
    ## cols(
    ##   old_locus_tag = col_character(),
    ##   locus_tag = col_character(),
    ##   `Gene product` = col_character(),
    ##   `gene name` = col_character(),
    ##   `Functional class` = col_character(),
    ##   `Functional colour` = col_double(),
    ##   `Transcriptome p-adj value` = col_character(),
    ##   `Transcriptome log2fold change in 630?erm?spo0A` = col_character(),
    ##   `Proteome log2fold change in 630?erm?spo0A` = col_character(),
    ##   `Mature spore proteome (PMID:19542279)` = col_character()
    ## )

``` r
d.sau <- d.sau %>% 
  mutate(sub_locus = str_extract(Gene, "\\..*") %>% str_remove("\\.")) 

d.sau <- d.sau %>% 
  mutate(pet_sub_locus =  LETTERS[as.numeric(d.sau$sub_locus)]) %>% 
  mutate(pet_locus = str_replace(Gene, "\\..*", pet_sub_locus)) %>% 
  select(-sub_locus, -pet_sub_locus) 

d.sau <- left_join(d.sau,
          d.pet %>% select(1:2),
          by = c("pet_locus" =  "old_locus_tag"))
```

``` r
library(gplots)
```

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
vn <- venn(list(DEM=d.dem$`Cd630 orthologue locus_tag` %>% unique(),
          FIM=d.fim$locus_tag %>% unique(),
          RAM=d.ramos$locus_tag %>% unique(),
          SAU=d.ramos$locus_tag %>% unique()))
```

![](Cdiff_sporGenes_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

I will include in the list of sporultion genes those that occur in
multiple data sets.

``` r
intersections <- attributes(vn)$intersections
keepers <- names(intersections)
keepers <- keepers[grep(":", keepers)]

spor_genes <- intersections[keepers] %>% unlist()
```

currently 329 genes. I will look at the genes left out to see if any
should be kept.

Which genes occur only in the Dembek data set?

``` r
look <- d.dem %>% 
  filter(`Cd630 orthologue locus_tag` %in% intersections$DEM)

look
```

    ## # A tibble: 656 x 5
    ##    locus_tag   `function`                     logFC  q.value `Cd630 orthologue ~
    ##    <chr>       <chr>                          <dbl>    <dbl> <chr>              
    ##  1 CDR20291_0~ seryl-tRNA synthetase          -3.53 2.97e- 2 CD630_00140        
    ##  2 CDR20291_0~ putative cytosine/adenosine ~ -10.8  3.18e- 9 CD630_00150        
    ##  3 CDR20291_0~ recombination protein          -9.64 3.48e- 6 CD630_00180        
    ##  4 CDR20291_0~ putative membrane protein      -9.58 1.39e- 6 CD630_00290        
    ##  5 CDR20291_0~ AraC-family transcriptional ~  -2.82 1.31e-14 CD630_00310        
    ##  6 CDR20291_0~ acetoin:2,6-dichlorophenolin~  -7.55 4.01e-13 CD630_00360        
    ##  7 CDR20291_0~ E2 component of acetoin dehy~  -2.30 1.17e- 7 CD630_00380        
    ##  8 CDR20291_0~ E3 component of acetoin dehy~  -2.12 4.80e-15 CD630_00390        
    ##  9 CDR20291_0~ putative dual-specificity pr~  -2.68 1.38e-10 CD630_00500        
    ## 10 CDR20291_0~ RNA polymerase sigma-H factor  -9.89 1.73e-39 CD630_00570        
    ## # ... with 646 more rows

There are two genes that by description are involved in sporulation. I
will keep those in the list.

``` r
look <- look %>% 
  filter(str_detect(`function`, "spore") |
           str_detect(`function`, "sporulation"))

look
```

    ## # A tibble: 2 x 5
    ##   locus_tag    `function`                logFC  q.value `Cd630 orthologue locus~
    ##   <chr>        <chr>                     <dbl>    <dbl> <chr>                   
    ## 1 CDR20291_33~ putative spore protein    -3.81 4.21e- 4 CD630_34940             
    ## 2 CDR20291_33~ stage V sporulation pro~ -11.4  4.38e-11 CD630_35160

``` r
spor_genes <- c(spor_genes, look$`Cd630 orthologue locus_tag`)
```

Which genes occur only in the Fimlaid data set?

``` r
look <- d.fim %>% 
  filter(locus_tag %in% intersections$FIM)

look
```

    ## # A tibble: 66 x 3
    ##    Name   locus_tag   description                                               
    ##    <chr>  <chr>       <chr>                                                     
    ##  1 CD2373 CD630_23730 CstA-like carbon starvation protein                       
    ##  2 murG   CD630_26510 UDP-NAG-NAM-(pentapeptide) pyrophosphoryl-undecaprenol NA~
    ##  3 drm    CD630_12230 phosphopentomutase                                        
    ##  4 CD1824 CD630_18240 P-type calcium transport ATPase                           
    ##  5 CD1404 CD630_14040 oligopeptide transporter                                  
    ##  6 CD3290 CD630_32900 hypothetical protein                                      
    ##  7 CD0145 CD630_01450 S1 RNA-binding domain-containing protein                  
    ##  8 CD0622 CD630_06220 hypothetical protein (DUF1629)                            
    ##  9 CD1222 CD630_12220 integrase site-specific recombinase XerD-like             
    ## 10 prfB   CD630_01440 peptide chain release factor 2 (RF-2)                     
    ## # ... with 56 more rows

Here too there are two genes that by description are involved in
sporulation. I will keep those in the list as well.

``` r
look <- look %>% 
  filter(str_detect(description, "spore") |
           str_detect(description, "sporulation"))

look
```

    ## # A tibble: 2 x 3
    ##   Name   locus_tag   description                          
    ##   <chr>  <chr>       <chr>                                
    ## 1 sigK   CD630_12300 sporulation factor ?K                
    ## 2 CD1045 CD630_10450 sporulation integral membrane protein

``` r
spor_genes <- c(spor_genes, look$locus_tag)
```

Currently I have 333 sporulation genes.

``` r
look <- d.pet %>% 
  filter(! locus_tag %in% spor_genes)

look %>% 
  filter(str_detect(`Gene product`, "spore") |
           str_detect(`Gene product`, "sporulation")) 
```

    ## # A tibble: 17 x 10
    ##    old_locus_tag locus_tag  `Gene product`          `gene name` `Functional cla~
    ##    <chr>         <chr>      <chr>                   <chr>       <chr>           
    ##  1 CD0007        CD630_000~ putative spore protein  <NA>        1.8.1           
    ##  2 CD1021        CD630_102~ putative spore coat pr~ <NA>        1.8.1           
    ##  3 CD1492        CD630_149~ Two-component sensor h~ <NA>        6.1.2           
    ##  4 CD1579        CD630_157~ Two-component sensor h~ <NA>        6.1.2           
    ##  5 CD1935        CD630_193~ Stage V sporulation pr~ spoVS       1.8.1           
    ##  6 CD2035        CD630_203~ putative sporulation i~ <NA>        1.8.1           
    ##  7 CD2273        CD630_227~ putative sporulation i~ <NA>        1.8.1           
    ##  8 CD2492        CD630_249~ Two-component sensor h~ <NA>        1.8.1           
    ##  9 CD2498        CD630_249~ D-alanyl-D-alanine car~ dacF1       1.8.1           
    ## 10 CD2681        CD630_268~ putative sporulation p~ <NA>        1.8.1           
    ## 11 CD2717        CD630_271~ putative sporulation p~ <NA>        1.8.1           
    ## 12 CD3397        CD630_339~ putative sporulation t~ whiA        1.8.1           
    ## 13 CD3498        CD630_349~ Stage V sporulation pr~ spoVB       1.8.1           
    ## 14 CD3548        CD630_354~ putative YaaT-like pro~ <NA>        1.8.1           
    ## 15 CD3671        CD630_367~ Stage 0 sporulation pr~ spo0J       1.8.1           
    ## 16 CD3672        CD630_367~ Transcriptional regula~ soj         1.8.1           
    ## 17 CD3673        CD630_367~ putative stage 0 sporu~ <NA>        1.8.1           
    ## # ... with 5 more variables: Functional colour <dbl>,
    ## #   Transcriptome p-adj value <chr>,
    ## #   Transcriptome log2fold change in 630?erm?spo0A <chr>,
    ## #   Proteome log2fold change in 630?erm?spo0A <chr>,
    ## #   Mature spore proteome (PMID:19542279) <chr>

``` r
look %>% 
  filter(str_detect(`gene name`, "spo")) 
```

    ## # A tibble: 3 x 10
    ##   old_locus_tag locus_tag  `Gene product`           `gene name` `Functional cla~
    ##   <chr>         <chr>      <chr>                    <chr>       <chr>           
    ## 1 CD1935        CD630_193~ Stage V sporulation pro~ spoVS       1.8.1           
    ## 2 CD3498        CD630_349~ Stage V sporulation pro~ spoVB       1.8.1           
    ## 3 CD3671        CD630_367~ Stage 0 sporulation pro~ spo0J       1.8.1           
    ## # ... with 5 more variables: Functional colour <dbl>,
    ## #   Transcriptome p-adj value <chr>,
    ## #   Transcriptome log2fold change in 630?erm?spo0A <chr>,
    ## #   Proteome log2fold change in 630?erm?spo0A <chr>,
    ## #   Mature spore proteome (PMID:19542279) <chr>

### Check if any other genes were left out

I downloade the reannotated C. Diff 630
([Genbank:CP016318.1](https://www.ncbi.nlm.nih.gov/nuccore/CP016318.1/))
genome as GFF. I will check if any genes that I donot have are noted as
related to sporulation.

``` r
library(ape)
gff <- here("spor_gene_list/data","/Cdiff.gff3")
d <- read.gff(gff)

d.cds <- d %>% filter(type == "CDS") 

d.parse.cds <- tibble()
for(i  in 1:nrow(d.cds[])){
  att <- d.cds$attributes[i]
  d.parse.cds <- 
    tibble(locus_tag = str_extract(att, "locus_tag=.*?;") %>% str_remove(";") %>% str_remove("locus_tag="),
           old_locus_tag = str_extract(att, "corresponds to CD630.*?;") %>% str_remove(";") %>% str_remove("corresponds to"),
           acc = str_extract(att, "Name=.*?;") %>% str_remove(";") %>% str_remove("Name="),
           gene = str_extract(att, "gene=.*?;") %>% str_remove(";") %>% str_remove("gene="),
           product = str_extract(att, "product=.*?;") %>% str_remove(";") %>% str_remove("product=")) %>% 
    bind_rows(d.parse.cds, .)
  
}

d.parse.cds$old_locus_tag <- trimws(d.parse.cds$old_locus_tag)

look <- d.parse.cds %>% 
  filter(! old_locus_tag %in% spor_genes)

look %>% 
  filter(str_detect(product, "spore") |
           str_detect(product, "sporulation")) 
```

    ## # A tibble: 17 x 5
    ##    locus_tag    old_locus_tag acc     gene  product                             
    ##    <chr>        <chr>         <chr>   <chr> <chr>                               
    ##  1 CDIF630erm_~ CD630_00070   ARE609~ <NA>  putative spore protein              
    ##  2 CDIF630erm_~ CD630_02750   ARE611~ splB  spore photoproduct (thymine dimer) ~
    ##  3 CDIF630erm_~ CD630_10210   ARE619~ <NA>  putative spore coat protein         
    ##  4 CDIF630erm_~ CD630_14920   ARE624~ <NA>  two-component sensor histidine kina~
    ##  5 CDIF630erm_~ CD630_15790   ARE625~ <NA>  two-component sensor histidine kina~
    ##  6 CDIF630erm_~ CD630_19350   ARE628~ spoVS stage V sporulation protein S       
    ##  7 CDIF630erm_~ CD630_20350   ARE629~ <NA>  putative sporulation integral membr~
    ##  8 CDIF630erm_~ CD630_22730   ARE631~ <NA>  putative sporulation integral membr~
    ##  9 CDIF630erm_~ CD630_24920   ARE634~ <NA>  two-component sensor histidine kina~
    ## 10 CDIF630erm_~ CD630_26810   ARE636~ <NA>  putative sporulation protein        
    ## 11 CDIF630erm_~ CD630_27170   ARE636~ <NA>  putative sporulation protein        
    ## 12 CDIF630erm_~ CD630_32710   ARE642~ <NA>  Spo0E-like sporulation regulatory p~
    ## 13 CDIF630erm_~ CD630_33970   ARE643~ whiA  putative sporulation transcription ~
    ## 14 CDIF630erm_~ CD630_34980   ARE644~ spoVB stage V sporulation protein B       
    ## 15 CDIF630erm_~ CD630_35480   ARE645~ <NA>  putative YaaT-like protein involved~
    ## 16 CDIF630erm_~ CD630_36710   ARE646~ spo0J stage 0 sporulation protein J       
    ## 17 CDIF630erm_~ CD630_36720   ARE646~ soj   sporulation initiation inhibitor pr~

There are 17 such genes! I’ll add them to the list.

``` r
spor_genes <- look %>% 
  filter(str_detect(product, "spore") |
           str_detect(product, "sporulation")) %>% 
  pull(old_locus_tag) %>% 
  c(spor_genes,.)
```

This brings me to 350 sporulation genes.

# get KEGG data on C. diff

## Strain 630

In KEGG this strain has the code *cdf*, and taxon number
[*T00487*](https://www.genome.jp/entry/T00487).

``` r
# all kegg genes
raw.kegg.cdf <- keggFind("genes", "cdf:CD630") 
d.kegg.cdf <- raw.kegg.cdf %>% 
  enframe(name = "kegg", value = "kegg.txt") %>% 
  separate(kegg, into = c("strain", "locus_tag"), sep = ":") %>% 
  separate(kegg.txt, into = c("symbol", "description"), sep = ";", fill = "left", extra = "merge") 


# KOs
ko <- keggLink("cdf", "ko")
ko <- enframe(ko, name = "ko", value = "cdf")
ko$cdf <-  str_replace(ko$cdf,pattern = "cdf:",replacement = "")
ko$ko <-  str_replace(ko$ko,pattern = "ko:",replacement = "")

# join
d.kegg.cdf <- left_join(d.kegg.cdf, ko, c("locus_tag" = "cdf"))
```

# sporulation genes with KO

``` r
spore_ko <- d.kegg.cdf %>% 
  filter (locus_tag %in% spor_genes) %>% 
  filter (! is.na(ko))

n.ko <- spore_ko$ko %>% unique() %>% length()
spore_ko
```

    ## # A tibble: 152 x 5
    ##    strain locus_tag   symbol  description                                  ko   
    ##    <chr>  <chr>       <chr>   <chr>                                        <chr>
    ##  1 cdf    CD630_00160 dnaX    " DNA polymerase III subunits gamma and tau" K023~
    ##  2 cdf    CD630_00170 <NA>    "DNA binding protein"                        K097~
    ##  3 cdf    CD630_00470 ispD    " 2-C-methyl-D-erythritol 4-phosphate cytid~ K009~
    ##  4 cdf    CD630_00480 ispF    " 2-C-methyl-D-erythritol 2,4-cyclodiphosph~ K017~
    ##  5 cdf    CD630_01060 cwlD    " Germination-specific N-acetylmuramoyl-L-a~ K014~
    ##  6 cdf    CD630_01190 glmM    " Phosphoglucosamine mutase"                 K034~
    ##  7 cdf    CD630_01200 glmS    " Glucosamine--fructose-6-phosphate aminotr~ K008~
    ##  8 cdf    CD630_01240 spoIID  " Stage II sporulation protein D"            K063~
    ##  9 cdf    CD630_01260 spoIIID " Stage III sporulation protein D"           K062~
    ## 10 cdf    CD630_02510 fliI    " ATP synthase subunit beta FliI"            K024~
    ## # ... with 142 more rows

Finally we are left with 152 Cdiff sporulation genes that correspond to
140 uniqe KOs.

``` r
write.csv(spore_ko, here("spor_gene_list/data", "cdif_spor_KOs.csv"))

# spor_genes[!spor_genes %in% d.kegg.cdf$locus_tag]
d.kegg.cdf %>% filter(locus_tag %in% spor_genes) %>% 
write.csv(here("spor_gene_list/data", "cdif_spor_KEGSs.csv"))
```

<!-- ## Strain R20291  -->
<!-- In KEGG this strain has the code *cdl*, and taxon number [*T00998*](https://www.genome.jp/entry/T00998). -->
<!-- ```{r} -->
<!-- # all kegg genes -->
<!-- raw.kegg.cdl <- keggFind("genes", "cdl:CDR20291")  -->
<!-- d.kegg.cdl <- raw.kegg.cdl %>%  -->
<!--   enframe(name = "kegg", value = "kegg.txt") %>%  -->
<!--   separate(kegg, into = c("strain", "locus_tag"), sep = ":") %>%  -->
<!--   separate(kegg.txt, into = c("symbol", "description"), sep = ";", fill = "left", extra = "merge")  -->
<!-- # KOs -->
<!-- ko <- keggLink("cdl", "ko") -->
<!-- ko <- enframe(ko, name = "ko", value = "cdl") -->
<!-- ko$cdl <-  str_replace(ko$cdl,pattern = "cdl:",replacement = "") -->
<!-- ko$ko <-  str_replace(ko$ko,pattern = "ko:",replacement = "") -->
<!-- # join -->
<!-- d.kegg.cdl <- left_join(d.kegg.cdl, ko, c("locus_tag" = "cdl")) -->
<!-- ``` -->
