###############################################################################
## PSMC - plot psmc
library(tidyverse)
library(magrittr)

source('scripts/Rscripts/demographic_history.R')

dir <- 'plots'

###############################################################################
## parameters
c_lvls <- factor(c("4+5*3+4", "4+30*2+4+6+10"))
s <- 100       ## Window size

cp <- c(RColorBrewer::brewer.pal(n = 11, 'RdYlBu')[c(9,10,11)],
        RColorBrewer::brewer.pal(n = 9, 'YlGnBu')[c(4)],
        RColorBrewer::brewer.pal(n = 9, 'YlOrRd')[5],
        RColorBrewer::brewer.pal(n = 9, 'YlOrRd')[3],
        RColorBrewer::brewer.pal(n = 9, 'YlOrRd')[8])

names(cp) <- c('Aipysurus laevis', 
               'Hydrophis melanocephalus',
               'Hydrophis curtus',
               'Laticauda laticaudata', 
               'Naja naja', 
               'Notechis scutatus', 
               'Pseudonaja textilis')
pc <- c('#000000', '#FF0000')
zc <- c('#74ADD1', '#FED976', '#BD0026')

###############################################################################
## Import data
tibble_psmc <- read_rds('results/psmc.Rds')
tibble_psmc_50 <- read_rds('results/psmc_z-clean.Rds')

###############################################################################
## Data preparation
tibble_bestRef <- pmap(tibble_psmc, function(gen, mu, reference, data){
  
  filter(data, sample == reference) %>%
    mutate(gen = gen,
           mu = mu,
           reference = reference) %>%
    select(reference, sample, gen, mu, clock, bootstrap, everything())
  
}) %>%
  bind_rows() %>%
  mutate(zchr = 'Z-chromosome included') %>%
  filter(clock == '4+30*2+4+6+10')

## Z-chromosome removed
tibble_bestRef_50 <- pmap(tibble_psmc_50, function(gen, mu, reference, data){
  
  filter(data, sample == reference) %>%
    mutate(gen = gen,
           mu = mu,
           reference = reference) %>%
    select(reference, sample, gen, mu, clock, bootstrap, everything())
  
}) %>%
  bind_rows() %>%
  mutate(zchr = 'Z-chromosome filtered') %>%
  filter(clock == '4+30*2+4+6+10')

## Subset for comparison
tibble_zchr <- bind_rows(tibble_bestRef, tibble_bestRef_50) %>%
  filter(reference %in% c('Notechis scutatus', 'Pseudonaja textilis', 'Hydrophis curtus')) %>%
  mutate(zchr = factor(x = zchr, levels = c('Z-chromosome included',
                                            'Z-chromosome filtered')))

###############################################################################
## Plotting

## Z-chromosome effect
plot_zchr_effect <- plotPSMC(tbl = tibble_zchr, 
                             g = 10, 
                             c = '4+30*2+4+6+10', 
                             colour_palette = cp) + 
  facet_wrap(. ~ zchr) +
  theme(strip.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(face = 'italic'))

png(filename = 'plots/clean_PSMC/z-chromosome-effect.png', width = 10, height = 8, units = "in", res = 500)
print(plot_zchr_effect)
dev.off()
