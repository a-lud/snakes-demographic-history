###############################################################################
## PSMC - plot ckean PSMC data
library(tidyverse)
library(magrittr)

source('scripts/Rscripts/demographic_history.R')

dir <- 'plots/clean_PSMC'
dir.create(dir)

###############################################################################
## parameters
c_lvls <- factor(c("4+5*3+4", "4+30*2+4+6+10"))
s <- 100       ## Window size

###############################################################################
## Import data + pre-processing
tibble_psmc <- read_rds(file = 'results/psmc_z-clean.Rds')

###############################################################################
## Colour palettes
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
###############################################################################
## Data preparation

## Original
df <- pmap(tibble_psmc, function(gen, mu, reference, data){
  
  filter(data, sample == reference) %>%
    mutate(gen = gen,
           mu = mu,
           reference = reference) %>%
    select(reference, sample, gen, mu, clock, bootstrap, everything())
  
}) %>%
  bind_rows()

###############################################################################
## Plotting

## Single plot - samples to own references - single clock
psmc_singlePlot <- plotPSMC(tbl = df, 
                            g = 10, 
                            c = '4+30*2+4+6+10', 
                            colour_palette = cp)

## Effect of mutation rate/generation time on curve position
psmc_scalingEffect <- plotPSMC_scalingEffect(tbl = df, 
                       c = '4+30*2+4+6+10', 
                       colour_palette = pc)
psmc_scalingEffect <- shift_legend(p = psmc_scalingEffect)

## Effect of reference selection on PSMC curves
plot_refOnNe <- plotXSMC_refEffect(tbl = tibble_psmc,
                                   colour_palette = cp,
                                   filterSample = 'Hydrophis curtus', 
                                   filterRef = 'Hydrophis curtus')

###############################################################################
## Saving plots

png(filename = paste(dir, 'psmc_gen10_4_30x2_4_6_10.png', sep = '/'), 
    width = 10, 
    height = 10, 
    units = 'in', 
    res = 500)
print(psmc_singlePlot)
dev.off()

ggsave(filename = 'psmc_scalingEffect_gen10_4_30x2_4_6_10.png', 
       plot = psmc_scalingEffect, device = 'png', 
       path = dir, 
       width = 10, 
       height = 10, 
       units = 'in', 
       dpi = 500)

png(filename = paste(dir, 'psmc_refEffect_gen10_4_30x2_4_6_10.png', sep = '/'), 
    width = 12, 
    height = 10, 
    units = 'in', 
    res = 300)
print(plot_refOnNe[[1]])
dev.off()
