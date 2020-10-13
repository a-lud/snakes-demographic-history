###############################################################################
## PSMC - plot psmc
library(tidyverse)
library(magrittr)

source('scripts/Rscripts/demographic_history.R')

dir <- 'plots'

###############################################################################
## parameters
c_lvls <- factor(c("4+5*3+4", "4+10*3+6+8", "4+25*2+4+6", "4+30*2+4+6+10"))
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
# zc <- c('#74ADD1', '#FED976', '#BD0026')

###############################################################################
## Import data
tibble_psmc <- read_rds(file = 'results/psmc.Rds') %>%
  filter(! reference %in% c('Laticauda colubrina', 'Hydrophis curtus')) %>%
  mutate(data = map(data, ~{
    .x %>%
      filter(sample != 'Laticauda colubrina')
  }))

###############################################################################
## Data preparation
tibble_bestRef <- pmap(tibble_psmc, function(gen, mu, reference, data){
  
  filter(data, sample == reference) %>%
    mutate(gen = gen,
           mu = mu,
           reference = reference) %>%
    select(reference, sample, gen, mu, clock, bootstrap, everything())
  
}) %>%
  bind_rows()

###############################################################################
## Plotting
plot_byRef <- plotPSMC_byRef(tbl = tibble_psmc, 
                             colour_palette = cp)

## Single plot - each sample to its own reference with laticaude removed
plot_bestRef_noLati <- plotPSMC_bestRef(tbl = tibble_bestRef,
                                        colour_palette = cp)

## Single plot - samples to own references - single clock
plot_bestRef_noLati_single_clock <- plotPSMC(tbl = tibble_bestRef, 
                                             g = 10, 
                                             c = '4+30*2+4+6+10',
                                             colour_palette = cp)

## Effect of reference on Ne estimate
plot_refOnNe <- plotXSMC_refEffect(tbl = tibble_psmc,
                                   colour_palette = cp,
                                   filterClock = c('4+10*3+6+8', '4+25*2+4+6'))

plot_scalingEffect <- plotPSMC_scalingEffect(tbl = tibble_bestRef, 
                                             c = '4+10*3+6+8', 
                                             colour_palette = pc, 
                                             filterSample = 'Laticauda colubrina')

plot_clockEffect <- plotPSMC_clockEffect(tbl = tibble_bestRef, 
                                         g = 10, 
                                         colour_palette = c('#ff0000', '#00ff00', '#0000ff', '#87cefa'), 
                                         filterSample = 'Laticauda colubrina')

###############################################################################
## Saving plots
map(names(plot_byRef), ~{
  
  gen_mu <- .x
  pth <- paste(dir, 'byReference', sep = '/')
  
  dir.create(path = pth,
             showWarnings = FALSE,
             recursive = TRUE)
  
  map(names(plot_byRef[[gen_mu]]), ~{
    n <- file.path(pth, paste0(.x,'_', gen_mu, '.pdf'))
    pdf(file = n, width = 10, height = 10)
    print(plot_byRef[[gen_mu]][[.x]])
    dev.off()
  })
  
})


map(names(plot_bestRef_noLati), ~{
  gen_mu <- .x
  pth <- paste(dir, 'bestRef', sep = '/')
  pth
  
  dir.create(path = pth,
             showWarnings = FALSE,
             recursive = TRUE)
  
  n <- file.path(pth, paste0(gen_mu, '.png'))
  png(filename = n, width = 10, height = 10, units = 'in', res = 300)
  print(plot_bestRef_noLati[[gen_mu]])
  dev.off()
})

map(names(plot_refOnNe), ~{
  gen_mu <- .x
  pth <- paste(dir, 'referenceEffect', sep = '/')
  
  dir.create(path = pth,
             showWarnings = FALSE,
             recursive = TRUE)
  
  n <- file.path(pth, paste0(gen_mu, '.png'))
  png(filename = n, width = 1000, height = 800)
  print(plot_refOnNe[[gen_mu]])
  dev.off()
})

png(filename = paste(dir, 'psmc_gen10_4_30x2_4_6_10.png', sep = '/'), width = 700, height = 700)
print(plot_bestRef_noLati_single_clock)
dev.off()

png(filename = paste(dir, 'psmc_10_4_30x2_4_6_10.scalingEffect.png', sep = '/'), width = 700, height = 700)
print(plot_scalingEffect)
dev.off()
