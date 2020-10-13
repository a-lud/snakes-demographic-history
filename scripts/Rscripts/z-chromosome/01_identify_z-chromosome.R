library(tidyverse)

## Functions
shift_legend <- function(p) {
  # check if p is a valid object
  if(!(inherits(p, "gtable"))){
    if(inherits(p, "ggplot")){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]), 
                               USE.NAMES = F)
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  
  # return repositioned legend
  lemon::reposition_legend(p, 'center', panel=names)
  
}

## Palette
palette <- c(RColorBrewer::brewer.pal(n = 11, 'RdYlBu')[c(9,10,11)],
             RColorBrewer::brewer.pal(n = 9, 'YlGnBu')[c(4)],
             RColorBrewer::brewer.pal(n = 9, 'YlOrRd')[5],
             RColorBrewer::brewer.pal(n = 9, 'YlOrRd')[3],
             RColorBrewer::brewer.pal(n = 9, 'YlOrRd')[8])

lvls <- c('Aipysurus laevis', 
          'Hydrophis melanocephalus',
          'Hydrophis curtus',
          'Laticauda laticaudata', 
          'Naja naja', 
          'Notechis scutatus', 
          'Pseudonaja textilis')

names(palette) <- lvls

## List MUMmer4 Coords files
coords <- list.files('z_chr_alignments/mummer4_alignments', pattern = '.coords', full.names = TRUE)
names(coords) <- sub(pattern = '.coords', replacement = '', x = basename(coords))

## Read in data
coords <- map(coords, ~{
  .x %>%
    read_tsv(
      col_names = c('start_ref', 'end_ref', 'start_query', 'end_query', 
                    'aln_ref', 'aln_query', 'pid', 'length_ref', 
                    'length_query', 'cov_ref', 'cov_query', 'ref', 'query'), 
      skip = 4
    ) %>%
    filter(pid >= 80)
})

## Tibble of scaffold lengths
seq_lengths <- map(coords, ~{
  .x %>%
    distinct(query, length_query)
})

## Collapse overlapping alignments: tibble > GRanges > tibble
unique_coords <- map(coords, ~{
  .x %>%
    select(contains('query')) %>%
    mutate(start = ifelse(start_query < end_query, start_query, end_query),
           end = ifelse(end_query > start_query, end_query, start_query )) %>%
    select(-c(start_query, end_query)) %>%
    rename(width = aln_query,
           length = length_query,
           coverage = cov_query, 
           chr = query) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    GenomicRanges::reduce() %>%
    GenomicRanges::as.data.frame() %>%
    as_tibble() %>%
    mutate(seqnames = as.character(seqnames)) %>%
    select(-strand)
})

## Total length of aligned regions
aln_lengths <- map(unique_coords, ~{
  .x %>%
    group_by(seqnames) %>%
    summarise(total_aln_len = sum(width))
})

## Join total aligned lengths with scaffold lengths
data <- map2(.x = aln_lengths, .y = seq_lengths, ~{
  left_join(.x, .y, by = c('seqnames' = 'query')) %>%
    mutate(prop_aln = total_aln_len/length_query * 100) %>%
    arrange(desc(prop_aln))
}) %>%
  bind_rows(.id = 'genome') %>%
  mutate(genome = case_when(
    genome == 'GCA_009733165.1_Nana_v5_genomic-subset' ~ 'Naja naja',
    genome == 'GCF_900518725.1_TS10Xv2-PRI_genomic' ~ 'Notechis scutatus',
    genome == 'GCF_900518735.1_EBS10Xv2-PRI_genomic' ~ 'Pseudonaja textilis',
    genome == 'Hcur1.v1.1' ~ 'Hydrophis curtus',
    genome == 'GCA_004320005.1_hydMel_1.0_genomic' ~ 'Hydrophis melanocephalus',
    genome == 'aipysurusLaevis' ~ 'Aipysurus laevis',
    genome == 'GCA_004320025.1_latLat_1.0_genomic' ~ 'Laticauda laticaudata'
  ),
  genome = factor(x = genome, levels = lvls))

data_50 <- filter(data, prop_aln >= 50)

## Proportion of total sequence length that has aligned: Histogram + density
histDensPid <- data %>%
  ggplot(aes(prop_aln)) +
  stat_bin(aes(y=..density..), 
           breaks = seq(min(data$prop_aln),
                        max(data$prop_aln), 
                        by = .1), 
           colour="black") +
  geom_line(stat="density", 
            size = 1,
            colour = 'red') +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  geom_vline(xintercept = 50, colour = 'red') +
  theme_bw() + 
  labs(x = 'Percentage aligned (%)',
       y = 'Density') +
  facet_wrap(.~genome, scales = 'free') +
  theme(strip.text = element_text(size = 14, face = 'italic'),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

png(filename = '~/Documents/kate/snakes-demographic-history/plots/z-chromosome/z-chromosome-densityHistogram.png', width = 10, height = 10, units = "in", res = 300)
print(histDensPid)
dev.off()

## Hist/Density of length distribution
histDensLen <- data %>%
  group_by(genome) %>%
  mutate(med = median(length_query)) %>%
  ungroup() %>%
  ggplot(aes(length_query)) +
  geom_histogram(bins = 100, 
                 fill = 'white', 
                 colour = 'black') +
  geom_histogram(data = data_50,
                 colour = 'red',
                 bins = 100) +
  theme_bw() +
  labs(x = bquote('Query length (bp - '*log[10]*')'),
       y = 'Frequency') +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_vline(aes(xintercept = med), 
             colour = 'red') +
  facet_wrap(.~genome, scales = 'free') +
  theme(strip.text = element_text(size = 14, face = 'italic'),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

png(filename =  '~/Documents/kate/snakes-demographic-history/plots/z-chromosome/z-chromosome-length-densityHistogram.png', width = 10, height = 10, units = 'in', res = 300)
print(histDensLen)
dev.off()

## Scatter plot: Query length x Query aligned bases
scatter <- data %>%
  group_by(genome) %>%
  mutate(xm = max(length_query),
         ym = max(total_aln_len)) %>%
  ungroup() %>%
  ggplot(aes(x = length_query,
             y = total_aln_len,
             colour = prop_aln)) +
  geom_point(size = 5) +
  geom_point(colour = 'black', size = 5, shape = 1) +
  geom_point(data = data_50, aes(length_query,
                                 total_aln_len,
                                 colour = prop_aln),
             size = 5) + 
  geom_point(data = data_50, 
             colour = 'red', 
             size = 5, 
             shape = 1) + 
  theme_bw() +
  labs(x = 'Query length (bp)',
       y = 'Query aligned bases (bp)',
       colour = 'Percentage aligned (%)') +
  scale_x_continuous(labels = scales::label_number_si()) +
  scale_y_continuous(labels = scales::label_number_si()) +
  # scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  viridis::scale_colour_viridis() +
  guides(colour = guide_colorbar(title.position = 'top',
                                 label.position = 'bottom',
                                 nrow = 1,
                                 title.hjust = 0.5)) +
  facet_wrap(.~genome,
             scales = 'free',
             ncol = 2) +
  theme(legend.direction = 'horizontal',
        strip.text = element_text(size = 14, face = 'italic'),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.width=unit(1.75, "cm"))

scatter <- shift_legend(scatter)

ggsave(filename =  'z-chromosome-scatter.png', 
       plot = scatter, 
       device = 'png', 
       path = '~/Documents/kate/snakes-demographic-history/plots/z-chromosome/', 
       width = 10, 
       height = 10, 
       units = 'in', 
       dpi = 300)

## Proportion of genome remaining using different threshold cutoffs
pidEffect <- map(seq(30, 100, 10), ~{
  data %>%
    group_by(genome) %>%
    mutate(genome_length = sum(length_query)) %>%
    filter(prop_aln < .x ) %>%
    mutate(remaining = sum(length_query),
           percentage_remaining = remaining/genome_length * 100,
           threshold_value = as.character(.x)) %>%
    ungroup() %>%
    distinct(genome, percentage_remaining,
             threshold_value)
}) %>%
  bind_rows() %>%
  mutate(threshold_value = factor(x = threshold_value,
                                  levels = seq(30, 100, 10))) %>%
  ggplot(aes(x = genome, 
             y = percentage_remaining,
             fill = genome)) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           colour = 'black') +
  scale_fill_manual(name = 'Genome',
                    values = palette) +
  labs(x = 'Genome',
       y = 'Genome retained (%)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1), 
        legend.title.align = 0.5) +
  facet_wrap(. ~ threshold_value,
             nrow = 2) +
  guides(fill = guide_legend(ncol = 2, 
                             title.position = 'top',
                             title.hjust = 0.5)) +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16, face = 'italic'),
        legend.title = element_text(size = 16),
        legend.position = 'bottom')

png(filename = 'plots/z-chromosome/z-chromosome-percentageGenomeRemain.png', width = 8, height = 10, units = 'in', res = 300)
print(pidEffect)
dev.off()

## Plotting raw length vs PID
# fragment_alignments <- map(names(coords), ~{
#   
#   df <- high_prop <- coords[[.x]] %>%
#     mutate(prop_aln = ((aln_query/length_query) * 100))
#   
#   high_prop <- filter(df, prop_aln >= 50)
#   
#   ggplot(data = df, 
#          aes(x = pid, 
#              y = aln_query, 
#              colour = prop_aln)) +
#     geom_point() +
#     geom_point(data = high_prop, 
#                aes(pid,
#                    aln_query,
#                    colour = prop_aln)) +
#     geom_point(data = high_prop, 
#                colour = 'red',
#                shape = 1) +
#     theme_bw() +
#     viridis::scale_colour_viridis() +
#     ggtitle(as.character(.x)) +
#     labs(x = 'Percentage identity (%)',
#          y = 'Length of alignment',
#          colour = '% of aligned query')
# })
# 
# pdf(file = 'plots/z-chromosome-scatter-fragmentAlignments.pdf', width = 10, height = 10)
# print(fragment_alignments)
# dev.off()

## Get sequence identifiers of Z-aligned scaffolds: 50% and 90%
map(names(data), ~{
  ids_50 <- data[[.x]] %>%
    filter(prop_aln < 50) %>%
    pull('seqnames')
  
  ids_90 <- data[[.x]] %>%
    filter(prop_aln < 90) %>%
    pull('seqnames')
  
  write_lines(x = ids_50, path = paste0('z_chr_alignments/psmc_filter_ids/', .x, '_50.txt'))
  write_lines(x = ids_90, path = paste0('z_chr_alignments/psmc_filter_ids/', .x, '_90.txt'))
})
