# spore_amg
Sporulation gene in viromes

(1) code/add-host-data.R

This code uses data from the virus-host database (https://www.genome.jp/virushostdb/) to add host taxonomy data to potential AMGs detected by DRAM-v. Joint data table written to 'data/Viruses/amg_summary_wHost.tsv' 

(2) code/compare-hosts.R

Using hypergeometric enrichment analysis to detct genes found in phages of Firmicute hosts, more than expected by random
