###############################################################################
## Alignment & variant statistics

library(tidyverse)
library(magrittr)

source('scripts/demographic_history.R')

###############################################################################
## Parameters
dir <- '01_bwa'

###############################################################################
## Data

# Coverage table
cov <- list.files(path = dir, pattern = '.tsv', full.names = TRUE, recursive = TRUE) %>%
  map(read_tsv, col_names = TRUE, col_types = cols()) %>%
  bind_rows() %>%
  mutate(reference = case_when(
    reference == 'GCF_900518725.1_TS' ~ 'Notechis scutatus',
    reference == 'GCF_900518735.1_EBS' ~ 'Pseudonaja textilis',
    reference == 'GCA_004320005.1_hydMel' ~ 'Hydrophis melanocephalus',
    reference == 'aipysurusLaevis' ~ 'Aipysurus laevis',
    reference == 'GCA_009733165.1_Nana' ~ 'Naja naja',
    reference == 'GCA_004320025.1_latLat' ~ 'Laticauda laticaudata',
    reference == 'GCA_004320045.1_latCor' ~ 'Laticauda colubrina'
  ),
  sample = case_when(
    sample == 'DRR144984_DRR144985' ~ 'Laticauda laticaudata',
    sample == 'DRR147552' ~ 'Laticauda colubrina',
    sample == 'DRR147394' ~ 'Hydrophis melanocephalus',
    sample == 'ERR2714264_ts' ~ 'Notechis scutatus',
    sample == 'ERR2714265' ~ 'Pseudonaja textilis',
    sample == 'KLS0691' ~ 'Aipysurus laevis',
    sample == 'SRR10428161' ~ 'Naja naja'
  ))

# % Reads aligned
aln <- list.files(path = dir, pattern = 'flagstat', full.names = TRUE, recursive = TRUE) %>%
  set_names(sub('.+/(.*)/(.*).flagstat', '\\1::\\2', .)) %>%
  map(ngsReports::importNgsLogs, type = 'flagstat') %>%
  bind_rows(.id = 'reference') %>%
  rename(sample = Filename) %>%
  mutate(reference = sub('::.*', '', reference),
         reference = case_when(
           reference == 'GCF_900518725.1_TS' ~ 'Notechis scutatus',
           reference == 'GCF_900518735.1_EBS' ~ 'Pseudonaja textilis',
           reference == 'GCA_004320005.1_hydMel' ~ 'Hydrophis melanocephalus',
           reference == 'aipysurusLaevis' ~ 'Aipysurus laevis',
           reference == 'GCA_009733165.1_Nana' ~ 'Naja naja',
           reference == 'GCA_004320025.1_latLat' ~ 'Laticauda laticaudata',
           reference == 'GCA_004320045.1_latCor' ~ 'Laticauda colubrina'
         ),
         sample = sub('.flagstat', '', sample),
         sample = case_when(
           sample == 'DRR144984_DRR144985' ~ 'Laticauda laticaudata',
           sample == 'DRR147552' ~ 'Laticauda colubrina',
           sample == 'DRR147394' ~ 'Hydrophis melanocephalus',
           sample == 'ERR2714264_ts' ~ 'Notechis scutatus',
           sample == 'ERR2714265' ~ 'Pseudonaja textilis',
           sample == 'KLS0691' ~ 'Aipysurus laevis',
           sample == 'SRR10428161' ~ 'Naja naja'
         )) %>%
  filter(flag %in% c('in total','mapped', 'properly paired')) %>%
  select(-`QC-failed`) %>%
  pivot_wider(names_from = flag, values_from = `QC-passed`) %>%
  mutate(prop_aln = (mapped/`in total`) * 100,
         prop_pp = (`properly paired`/`in total`) * 100) %>%
  rename(total_reads = `in total`, properly_paired = `properly paired`)

# VCF statistics
vcf <- list.files(path = dir, pattern = '.vcf.stats', full.names = TRUE, recursive = TRUE) %>%
  set_names(sub('.+/(.*)/(.*).vcf.stats', '\\1::\\2', .)) %>%
  readBcftoolsStats() %>%
  extract2('SN') %>%
  mutate(reference = case_when(
    reference == 'GCF_900518725.1_TS' ~ 'Notechis scutatus',
    reference == 'GCF_900518735.1_EBS' ~ 'Pseudonaja textilis',
    reference == 'GCA_004320005.1_hydMel' ~ 'Hydrophis melanocephalus',
    reference == 'aipysurusLaevis' ~ 'Aipysurus laevis',
    reference == 'GCA_009733165.1_Nana' ~ 'Naja naja',
    reference == 'GCA_004320025.1_latLat' ~ 'Laticauda laticaudata',
    reference == 'GCA_004320045.1_latCor' ~ 'Laticauda colubrina'
  ),
  sample = case_when(
    sample == 'DRR144984_DRR144985' ~ 'Laticauda laticaudata',
    sample == 'DRR147552' ~ 'Laticauda colubrina',
    sample == 'DRR147394' ~ 'Hydrophis melanocephalus',
    sample == 'ERR2714264_ts' ~ 'Notechis scutatus',
    sample == 'ERR2714265' ~ 'Pseudonaja textilis',
    sample == 'KLS0691' ~ 'Aipysurus laevis',
    sample == 'SRR10428161' ~ 'Naja naja'
  )) %>%
  filter(key %in% c('number of records:', 'number of SNPs:')) %>%
  pivot_wider(names_from = key, values_from = value) %>%
  select(-c(condition, id)) %>%
  rename(n_sites = `number of records:`, n_SNPs = `number of SNPs:`)

# Joining data frames
tbl <- reduce(list(cov, aln, vcf), left_join)

###############################################################################
## Coverage tables

# All samples to each reference
all_cov <- tbl %>%
  select(reference, sample, total_reads, mapped,avg_depth, breadth_depth, prop_aln, prop_pp) %>%
  mutate(double_avg = avg_depth * 2,
         third_avg = avg_depth / 3) %>%
  select(reference, sample, total_reads, mapped, avg_depth, double_avg, third_avg, breadth_cov = breadth_depth, prop_aln, prop_properPair = prop_pp)

# All samples max/min coverage/alignment by reference
all_cov_max_min <- all_cov %>%
  filter(reference != sample) %>%
  group_by(reference) %>%
  mutate(avg_aln = sum(mapped)/sum(total_reads) * 100,
         max_avg = max(avg_depth),
         min_avg = min(avg_depth),
         max_breadth = max(breadth_cov),
         min_breadth = min(breadth_cov),
         min_aln = min(prop_aln),
         max_aln = max(prop_aln),
         min_aln_properPair = min(prop_properPair),
         max_aln_properPair = max(prop_properPair)) %>%
  ungroup() %>%
  select(reference, contains(match = 'max'), contains('min'), avg_aln) %>%
  distinct()

# Samples to respective reference genomes
ref_cov <- tbl %>%
  select(reference, sample, avg_depth, breadth_depth, prop_aln, prop_pp) %>%
  filter(reference == sample) %>%
  mutate(double_avg = avg_depth * 2,
         third_avg = avg_depth / 3) %>%
  select(reference_sample = reference, avg_depth, double_avg, third_avg, breadth_cov = breadth_depth, prop_aln, prop_pp)

# Samples to respective reference genomes: max/min coverage comparison between samples
ref_cov_max_min <- ref_cov %>%
  mutate(max_avg = max(avg_depth),
         min_avg = min(avg_depth),
         max_breadth = max(breadth_cov),
         min_breadth = min(breadth_cov)) %>%
  select(contains(match = 'max'), contains('min')) %>%
  distinct()

###############################################################################
## Variant tables

# All samples to each reference
all_var <- tbl %>%
  select(reference, sample, n_sites, n_SNPs) %>%
  rename()

# Samples to respective reference genomes
ref_var <- tbl %>%
  select(reference, sample, n_sites, n_SNPs) %>%
  filter(reference == sample)

###############################################################################
## Plotting - Coverage
ref_cov %>%
  select(-c(double_avg, third_avg)) %>%
  pivot_longer(-reference_sample, names_to = 'desc', values_to = 'values') %>%
  mutate(desc = case_when(
    desc == 'avg_depth' ~ 'Avg. depth',
    desc == 'breadth_cov' ~ 'Avg. breadth',
    desc == 'prop_aln' ~ '% reads aligned',
    desc == 'prop_pp' ~ '% read properly paired'
  ),
  desc = factor(x = desc, levels = c('Avg. depth',
                                     'Avg. breadth',
                                     '% reads aligned',
                                     '% read properly paired'))) %>%
  ggplot(aes(x = reference_sample, 
             y = values, 
             fill = reference_sample)) +
  geom_bar(stat = 'identity', 
           position = 'dodge',
           colour = 'black') +
  scale_y_continuous(breaks = seq(0, max(ref_cov$avg_depth), 20),
                     expand = c(0,2.5)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 11, name = 'Spectral')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, 
                                   size = 12),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +
  facet_wrap(.~desc)
