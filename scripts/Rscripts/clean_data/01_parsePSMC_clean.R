###############################################################################
## PSMC - parse data
library(tidyverse)
library(magrittr)

source('scripts/Rscripts/demographic_history.R')

###############################################################################
## Parameters
c_lvls <- factor(c("4+5*3+4", "4+30*2+4+6+10"))
mu <- c(1.25e-8, 7.2e-09)  ## Mutation 1.25e-8 7.2e-09
g <- c(10, 3)        ## Generation 10 3
s <- 100       ## Window size

## All combinations
params <- list(
  mutation = mu,
  generation = g,
  windowSize = s
) %>%
  cross_df() %>%
  slice(c(1,4))

###############################################################################
## Input data - main + bootstrap combined files
files_psmcCombined <- list.files(path = '03_analysis/psmc_50', 
                                 pattern = '-combined.psmc', 
                                 full.names = TRUE, 
                                 recursive = TRUE) %>%
  set_names(sub('.+/(.*)/(.*)-combined.psmc', '\\1::\\2', .))

###############################################################################
## Data frame of PSMC results for plotting
tibble_psmc <- pmap(params, function(mutation, generation, windowSize){
  
  map(files_psmcCombined, psmc.result, mu = mutation, s = windowSize, g = generation) %>%
    bind_rows(.id = 'temp') %>%
    separate(col = temp, into = c('reference', 'sample'), sep = '::') %>%
    separate(col = sample, into = c('sample', 'clock'), sep = '-') %>%
    mutate(clock = gsub('_', '+', clock),
           clock = gsub('x', '*', clock),
           mu = mutation,
           gen = generation) %>%
    select(mu, gen, reference, sample, clock, bootstrap, YearsAgo, Ne)
  
}) %>%
  bind_rows() %>%
  mutate(clock = factor(x = clock, levels = c_lvls),
         bootstrap = as.numeric(bootstrap),
         lt = ifelse(bootstrap == 0, 'Standard', 'Bootstrapped'),
         reference = case_when(
           reference == 'GCF_900518725.1_TS' ~ 'Notechis scutatus',
           reference == 'GCF_900518735.1_EBS' ~ 'Pseudonaja textilis',
           reference == 'GCA_004320005.1_hydMel' ~ 'Hydrophis melanocephalus',
           reference == 'aipysurusLaevis' ~ 'Aipysurus laevis',
           reference == 'GCA_009733165.1_Nana' ~ 'Naja naja',
           reference == 'GCA_004320025.1_latLat' ~ 'Laticauda laticaudata',
           reference == 'Hcur1.v1.1' ~ 'Hydrophis curtus'
         ),
         sample = case_when(
           sample == 'DRR144984_DRR144985' ~ 'Laticauda laticaudata',
           sample == 'DRR147394' ~ 'Hydrophis melanocephalus',
           sample == 'ERR2714264_ts' ~ 'Notechis scutatus',
           sample == 'ERR2714265' ~ 'Pseudonaja textilis',
           sample == 'KLS0691' ~ 'Aipysurus laevis',
           sample == 'SRR10428161' ~ 'Naja naja',
           sample == 'SRR10861675_SRR10861676' ~ 'Hydrophis curtus'
         )) %>%
  group_by(reference, gen, mu) %>%
  nest()

## Save object to prevent re-running
write_rds(x = tibble_psmc, 
          file = 'results/psmc_z-clean.Rds', 
          compress = 'gz')
