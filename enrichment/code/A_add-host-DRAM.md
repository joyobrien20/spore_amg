Untitled
================
Daniel Schwartz

# Goal

To match Refeq phages to hosts. This will enable to assign a phage
infecting a sporulator or a non-sporulator host.

``` r
# list of all phages processed by DRAM
d.vir <- read_tsv(here("enrichment","data/Viruses/vMAG_stats.tsv"))
```

    ## Rows: 3689 Columns: 15
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr  (1): refseq.id
    ## dbl (12): VIRSorter category, Prophage, Gene count, Strand switches, potenti...
    ## lgl  (2): Circular, Transposase present
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#import virus-host data
vh.db <- read_tsv(here("enrichment","data","virushostdb.tsv") ) 
```

    ## Rows: 16617 Columns: 14
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\t"
    ## chr (11): virus.name, virus.lineage, refseq.id, KEGG.GENOME, KEGG.DISEASE, D...
    ## dbl  (3): virus.tax.id, host.tax.id, source.organism
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# #-----------------------------#
# # Downloading virus-host database data
# #-----------------------------#
# # downloaded on 4/May/2021
# library(curl)
# url <- "ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv"
# vh.db <- read_tsv(curl(url))
# 
# vh.db <- rename_all(vhdb, make.names)
# 
# write_tsv(vh.db, here("data","virushostdb.tsv"))
# #-----------------------------#
```

## match data by taxid

extracting from vh.db a mapping-table between refseq and taxid. There
are many viruses that have multiple refseq ids (=multiple sequences),
and there are viruses that have multiple hosts (= multiple rows).

``` r
map_taxid_refseq <- vh.db %>% 
  select(virus.tax.id, refseq.id) %>% 
  # get all sequence IDs
  separate(refseq.id, into = as.character(1:200), sep = ",",fill = "right") %>% 
  pivot_longer(-virus.tax.id, names_to = "num", values_to = "refseq.id") %>% 
  filter(! is.na(refseq.id)) %>% 
  select(-num) %>% 
  # discard replicated ID pairs
  distinct() %>% 
  mutate(refseq.id = trimws(refseq.id))
```

## organize refseq taxID’s

``` r
d.vir <- d.vir %>% 
  # parse viral ID to match vhdb
  mutate(refseq.id = str_remove(refseq.id, "\\..*$")) %>% 
  mutate(refseq.id = trimws(refseq.id))
```

Are we missing any viruses in vh.db?

``` r
sum(!d.vir$refseq.id %in% map_taxid_refseq$refseq.id)
```

    ## [1] 3

``` r
d.vir$refseq.id[which(!d.vir$refseq.id %in% map_taxid_refseq$refseq.id)]
```

    ## [1] "NC_029050" "NC_029072" "NC_042059"

Only 3! We will do without them.

``` r
#add taxon ID to viral amg data
d.vir <-left_join(d.vir, map_taxid_refseq, by = "refseq.id") 
```

# Dealing with multiple rows(=hosts) per virus in vhdb

-   We only need to keep virus-host data that is relavent to refseq
    viruses
-   We will keep only viruses of bacteria (phages)

``` r
# keep only vh.db data relavent to refseq
vh.db <- semi_join(vh.db, d.vir, by = c("virus.tax.id"))

# keep only phages of bacteria
vh.db <- vh.db %>% 
  filter(str_detect(host.lineage, "^Bacteria"))
```

How much duplication is there?

``` r
vh.db%>%
  group_by(virus.tax.id, virus.name)%>%
  summarise(N_hosts=n())%>%
  group_by(N_hosts) %>% 
  summarise(n=n()) %>% 
  arrange(desc(N_hosts))
```

    ## `summarise()` has grouped output by 'virus.tax.id'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 6 x 2
    ##   N_hosts     n
    ##     <int> <int>
    ## 1       6     2
    ## 2       5     1
    ## 3       4     2
    ## 4       3    15
    ## 5       2   232
    ## 6       1  3366

Roughly 10% of phages have multiple hosts listed.

How much variation is there in multiple host taxonomy ?

Typically phages will infect very closely related hosts (sibling strains
or species).

First we will break up the host lineage data so we can check if hosts
are of similar taxonomy. We have made sporulation predictions at the
family level so that is the level of interest. However, there is some
cleaning to do first, since the the host lineage specification is not
uniform in levels specified.

Some have an intermediate rank between domain and phylum specifically
FCB group AND Terrabacteria group. I will remove those for consistency
and then separate out the taxonomy.

``` r
vh.db <- vh.db%>%
  mutate(host.lineage=str_remove(host.lineage," FCB group;"))%>%
  mutate(host.lineage=str_remove(host.lineage," Terrabacteria group;"))

vh.db <- vh.db%>%
  separate(host.lineage, sep="; ",
           into = c("domain", "phylum", "class", "order", "family","genus.etc"), 
           extra = "merge", remove = F)
```

    ## Warning: Expected 6 pieces. Missing pieces filled with `NA` in 240 rows [16,
    ## 118, 198, 304, 391, 407, 409, 486, 497, 521, 523, 524, 552, 575, 598, 608, 610,
    ## 654, 696, 730, ...].

Are there any viruses with hosts of different families?

``` r
vh.db%>%
  group_by(virus.tax.id, virus.name, family)%>%
  summarise(n=n())%>%
  arrange(desc(n)) %>% 
  group_by(virus.tax.id, virus.name)%>%
  summarise(n=n())%>%
  filter(n>1) %>% 
  arrange(desc(n))
```

    ## `summarise()` has grouped output by 'virus.tax.id', 'virus.name'. You can
    ## override using the `.groups` argument.
    ## `summarise()` has grouped output by 'virus.tax.id'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 12 x 3
    ## # Groups:   virus.tax.id [12]
    ##    virus.tax.id virus.name                         n
    ##           <dbl> <chr>                          <int>
    ##  1        10678 Escherichia virus P1               2
    ##  2        10679 Escherichia virus P2               2
    ##  3        10693 Enterobacteria phage RB51          2
    ##  4       373407 Mycobacterium virus Halo           2
    ##  5       373415 Mycobacterium phage Wildcat        2
    ##  6       440250 Phormidium virus WMP3              2
    ##  7       663560 Mycobacterium phage UncleHowie     2
    ##  8       981323 Gordonia phage GTE2                2
    ##  9      1041524 Escherichia phage K30              2
    ## 10      1100814 Gordonia phage GTE7                2
    ## 11      1352230 Thermus phage phi OH2              2
    ## 12      2681603 Escherichia phage Mu               2

``` r
#keep IDs of viruses with multi-family hosts 
multi_fam_id <- vh.db%>%
  group_by(virus.tax.id, virus.name, family)%>%
  summarise(n=n())%>%
  arrange(desc(n)) %>% 
  group_by(virus.tax.id, virus.name)%>%
  summarise(n=n())%>%
  filter(n>1) %>% 
  arrange(desc(n)) %>% 
  pull(virus.tax.id)
```

    ## `summarise()` has grouped output by 'virus.tax.id', 'virus.name'. You can
    ## override using the `.groups` argument.
    ## `summarise()` has grouped output by 'virus.tax.id'. You can override using the
    ## `.groups` argument.

There are only twelve such viruses. We will get back to those. But
first, ee will keep a single host for all the phages with hosts that
converge at the family level. Some phages have refseq in evidence for
one of their hosts, I think this is the isolation host, so we will
prioritize that row.

``` r
#remove host duplicate
vh.db_clean <- vh.db %>% 
  # remove multi-family, deal with below
  filter(!virus.tax.id %in% multi_fam_id) %>% 
  # mark refseq evidence as frist priority
  mutate(priority1 = grepl("RefSeq", .$evidence)) %>% 
  # count host taxonomy levels, for second priority
  mutate(priority2 = str_count(host.lineage,";")) %>% 
  arrange(desc(priority1), desc(priority2)) %>% 
  # remove duplicated
  mutate(dup = duplicated(virus.tax.id)) %>% 
  filter(!dup) %>% 
  select(-priority1, -priority2, -dup)

# check duplication
duplicated(vh.db_clean$virus.tax.id)%>%sum() 
```

    ## [1] 0

Now we deal with viruses that have more than one family iin host
assignment.

``` r
vh.db%>%
  filter(virus.tax.id %in% multi_fam_id) %>% 
  select(virus.tax.id, virus.name, family)%>%
  distinct() %>%
  arrange(virus.tax.id)
```

    ## # A tibble: 24 x 3
    ##    virus.tax.id virus.name                  family            
    ##           <dbl> <chr>                       <chr>             
    ##  1        10678 Escherichia virus P1        <NA>              
    ##  2        10678 Escherichia virus P1        Enterobacteriaceae
    ##  3        10679 Escherichia virus P2        <NA>              
    ##  4        10679 Escherichia virus P2        Enterobacteriaceae
    ##  5        10693 Enterobacteria phage RB51   <NA>              
    ##  6        10693 Enterobacteria phage RB51   Enterobacteriaceae
    ##  7       373407 Mycobacterium virus Halo    Mycobacteriaceae  
    ##  8       373407 Mycobacterium virus Halo    <NA>              
    ##  9       373415 Mycobacterium phage Wildcat Mycobacteriaceae  
    ## 10       373415 Mycobacterium phage Wildcat <NA>              
    ## # ... with 14 more rows

Most have a family assigned only once and NA in the other.

For those with a host having NA in family we will keep the other row,
with family assigned.

``` r
keep_id <- vh.db%>%
  filter(virus.tax.id %in% multi_fam_id) %>% 
  select(virus.tax.id, virus.name, family)%>%
  distinct() %>%
  filter(is.na(family)) %>% 
  pull(virus.tax.id)

# add to clean list
vh.db_clean <- 
vh.db%>%
  filter(virus.tax.id %in% keep_id) %>% 
   distinct() %>% 
  filter(!is.na(family)) %>% 
  # remove duplicated
  mutate(dup = duplicated(virus.tax.id)) %>% 
  filter(!dup) %>% 
  select(-dup) %>% 
  bind_rows(vh.db_clean,.)
 
# check duplication
duplicated(vh.db_clean$virus.tax.id)%>%sum() 
```

    ## [1] 0

What is left?

``` r
multi_fam_id <- vh.db%>%
  group_by(virus.tax.id, virus.name, family)%>%
  filter(!is.na(family)) %>% 
  summarise(n=n())%>%
  arrange(desc(n)) %>% 
  group_by(virus.tax.id, virus.name)%>%
  summarise(n=n())%>%
  filter(n>1) %>% 
  arrange(desc(n)) %>% 
  pull(virus.tax.id)
```

    ## `summarise()` has grouped output by 'virus.tax.id', 'virus.name'. You can
    ## override using the `.groups` argument.
    ## `summarise()` has grouped output by 'virus.tax.id'. You can override using the
    ## `.groups` argument.

``` r
vh.db%>%
  filter(virus.tax.id %in% multi_fam_id) %>% 
  select(virus.tax.id, virus.name, phylum,family)%>%
  distinct() %>%
  arrange(virus.tax.id)
```

    ## # A tibble: 8 x 4
    ##   virus.tax.id virus.name            phylum                              family 
    ##          <dbl> <chr>                 <chr>                               <chr>  
    ## 1       440250 Phormidium virus WMP3 Cyanobacteria/Melainabacteria group Leptol~
    ## 2       440250 Phormidium virus WMP3 Cyanobacteria/Melainabacteria group Oscill~
    ## 3       981323 Gordonia phage GTE2   Actinobacteria                      Nocard~
    ## 4       981323 Gordonia phage GTE2   Actinobacteria                      Gordon~
    ## 5      1100814 Gordonia phage GTE7   Actinobacteria                      Nocard~
    ## 6      1100814 Gordonia phage GTE7   Actinobacteria                      Gordon~
    ## 7      1352230 Thermus phage phi OH2 Deinococcus-Thermus                 Therma~
    ## 8      1352230 Thermus phage phi OH2 Firmicutes                          Bacill~

We are left with four phages, but only one has the potential to infect a
sporulator, *Thermus phage phi OH2*. According to
[VHDB](https://www.genome.jp/virushostdb/1352230) the evidence for
*Thermaceae* infection is “*Unpublished, from submitter provided text in
sequence record AB823818.*”, For Geobacillus infection there is a
sequenced prophage in the host genome. I will keep the latter. For the
other three phages I will discard dupplicates as above.

``` r
# add to clean list
vh.db_clean <-
vh.db%>%
  filter(virus.tax.id %in% multi_fam_id) %>%   
  # remove Thermus host
  filter(! family =="Thermaceae") %>% 
  # mark refseq evidence as frist priority
  mutate(priority1 = grepl("RefSeq", .$evidence)) %>% 
  # count host taxonomy levels, for second priority
  mutate(priority2 = str_count(host.lineage,";")) %>% 
  arrange(desc(priority1), desc(priority2)) %>% 
  # remove duplicated
  mutate(dup = duplicated(virus.tax.id)) %>% 
  filter(!dup) %>% 
  select(-priority1, -priority2, -dup) %>% 
    bind_rows(vh.db_clean,.)
 
# check duplication
duplicated(vh.db_clean$virus.tax.id)%>%sum() 
```

    ## [1] 0

# add hosts to refseq phages

``` r
d.vir <- left_join(d.vir, vh.db_clean, by = "virus.tax.id")

write_csv(d.vir, here("enrichment","data/Viruses/refseq_phages_wHost.csv"))
```
