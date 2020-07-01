###############################################################################
## PSMC - Functions for demographic history analysis
## 
## Author: Alastair Ludington
## Institution: University of Adelaide

###############################################################################
## Functions

## Parse PSMC results
psmc.result <- function(file, mu, s = 100, g){
  
  ## Rounds start at 0 - 31 rounds in total
  RD30 <- 31
  
  ## Read PSMC file as character vector
  X <- scan(file=file,
            what="",
            sep="\n",
            quiet=TRUE)
  
  ## Start and end of information chunk - flanked by RD and '\\'
  START <- grep(pattern = '^RD', x = X)
  END <- grep(pattern = '//', x = X)
  
  ## Need to be the same length
  if(length(START) != length(END)){stop('Number of blocks do not match')}
  
  ## List of chunks
  BLOCKS <- map2(.x = START, .y = END, ~{
    X[.x:.y]
  })
  
  ## Split list up into rounds - Each file has 1 main + 100 bootstrap files
  BLOCKS <- split(BLOCKS, ceiling(seq_along(BLOCKS)/RD30)) %>%
    set_names(nm = 0:100)
  
  ## Subset main and bootstrap files for round 30 data
  map(BLOCKS, ~{
    
    ## Get 30th iteration
    chunk <- .x[[RD30]]
    
    TR <- grep("^TR", chunk, value=TRUE)
    RS <- grep("^RS", chunk, value=TRUE)
    
    theta0 <- as.numeric(read.table(text = TR, header = FALSE)[1,2])
    N0 <- theta0/4/mu/s
    
    a <- read.table(text = RS, header = FALSE)
    Generation <- as.numeric(2*N0*a[,3])
    Ne <- as.numeric(N0*a[,4])
    
    n.points <- length(Ne)
    YearsAgo <- c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
                  Generation[n.points])*g
    Ne <- c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
            Ne[n.points])
    
    tibble(YearsAgo = YearsAgo, Ne = Ne)
    
  }) %>%
    bind_rows(.id = 'bootstrap')
  
}

## Parse bcftools stats files
readBcftoolsStats <- function(fileList){
  
  ## Files to vectors
  files <- map(fileList, scan, what = 'character', sep = '\n', quiet = TRUE)
  
  ## Data to subset
  s <- c('SN', 'TSTV', 'SiS', 'AF', 'QUAL', 'IDD', 'ST', 'DP')
  
  c <- cross2(.x = files, .y = s)
  names(c) <- cross2(names(files), s) %>% 
    map(~{
      unlist(.x) %>%
        paste0(collapse = '::')
    })
  
  list <- map(c, ~{
    f <- .x[[1]]
    v <- paste0('^', .x[[2]])
    
    ## Find matches
    r <- grep(v, f)
    
    if(length(r) != 0){
      
      h <- r[[1]] - 1
      h <- f[h]
      
      ## Cleaning up header
      h <- str_remove(string = h, pattern = '# ')
      h <- gsub("\\[|\\]|[[:digit:]]", '', h)
      
      ## Data
      r <- f[r]
      read_tsv(file = c(h,r), col_names = TRUE, col_types = cols())
      
    } else {
      return(NULL)
    }
    
  })
  
  ## Bind data
  map(s, ~{
    list[str_detect(string = names(list), pattern = paste0(.x, '$'))]
  }) %>%
    set_names(s) %>%
    map(~{
      d <- .x %>%
        bind_rows(.id = 'sample')
      if(nrow(d) != 0){
        d %>%
          select(-2) %>%
          separate(col = sample, into = c('reference', 'sample', 'condition'), sep = '::')
      }
    })
}

## Plot effect of selected reference on Ne
plotXSMC_refEffect <- function(tbl, yearFilter = NULL, colour_palette, filterSample=NULL, filterRef=NULL, filterClock = NULL){
  
  ne <- tbl %>%
    unnest(cols = data) %>%
    ungroup() %>%
    group_by(gen, mu) %>%
    nest()
  
  nms <- ne %>%
    pmap(function(gen, mu, data){
      paste(gen, mu, sep = '_')
    }) %>%
    unlist()
  
  pmap(ne, function(gen, mu, data){
    
    if(!is.null(filterSample)){
      data <- filter(data, ! sample %in% filterSample)
    }
    
    if(!is.null(filterRef)){
      data <- filter(data, ! reference %in% filterRef)
    }
    
    if(!is.null(filterClock)){
      data <- filter(data, ! clock %in% filterClock)
    }
    
    # if(inputSource == 'psmc') {
    #   mt <- 'PSMC'
    # } else {
    #   mt <- 'MSMC'
    # }
    # st <- paste('Generation:', gen, 'Mutation:', mu, sep = ' ')
    
    if(!is.null(yearFilter)){
      data <- filter(data, yearsAgo >= yearFilter)
    }
    
    p <- ggplot(data = data,
                aes(x = YearsAgo, y = Ne,
                    colour = reference,
                    alpha = lt == 'Standard',
                    group = interaction(reference, bootstrap))) +
      geom_step() +
      scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
    
    if(gen == 10){
      p <- p + labs(x = expression(paste("Years (g = 10, mu = ",1.25, 'x', 10^-8, ")")),
                    y = expression(paste('Effective population size (N'[e], ')')))
    } else {
      p <- p + labs(x = expression(paste("Years (g = 3, mu = ",7.2, 'x', 10^-9, ")")),
                    y = expression(paste('Effective population size (N'[e], ')')))
    }
    
    p + 
      scale_colour_manual(values = colour_palette) +
      scale_alpha_manual(values = c(0.1, 1), guide = FALSE) +
      scale_size_manual(values = c(0.4, 1)) +
      theme_bw() +
      theme(strip.text.x = element_text(size = 13),
            strip.text.y = element_text(size = 13),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16),
            legend.position = 'bottom',
            legend.title = element_blank(),
            legend.text = element_text(face = 'italic', size = 16)) +
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      # annotation_logticks() +
      facet_grid(clock ~ sample,
                 scales = 'free_y')
  }) %>%
    set_names(value = nms)
}

## Plot PSMC by reference
plotPSMC_byRef <- function(tbl, max_Ne = NULL, filterSample = NULL, filterRef = NULL, colour_palette){
  
  tbl <- tbl %>%
    ungroup() %>%
    group_by(mu, gen) %>%
    group_split()
  
  nms <- tbl %>%
    map(~{
      .x %>%
        select(mu, gen) %>%
        distinct() %>%
        unite(col = 'nm', mu:gen, sep = '_') %>%
        extract2('nm')
    }) %>%
    unlist()
  
  names(tbl) <- nms
  
  map(tbl, ~{
    
    df <- .x
    
    if(!is.null(filterRef)){
      df <- filter(df, ! reference %in% filterRef)
    }
    
    pmap(df, function(gen, mu, reference, data){
      
      ## Apply max Ne filter
      if(!is.null(max_Ne)){
        data <- filter(data, Ne <= max_Ne)
      }
      
      if(!is.null(filterSample)){
        data <- filter(data, ! sample %in% filterSample)
      }
      
      ## Plotting PSMC - facet by atomic clock
      p <- ggplot(data = data,
                  aes(x = YearsAgo, y = Ne,
                      colour = sample,
                      alpha = lt == 'Standard',
                      group = interaction(sample, bootstrap))) +
        geom_step() +
        scale_x_log10(
          breaks = scales::trans_breaks("log10", function(x) 10^x),
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        scale_y_log10(
          breaks = scales::trans_breaks("log10", function(x) 10^x),
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
      
      if(gen == 10){
        p <- p + labs(x = expression(paste("Years (g = 10, mu = ",1.25, 'x', 10^-8, ")")),
                      y = expression(paste('Effective population size (N'[e], ')')))
      } else {
        p <- p + labs(x = expression(paste("Years (g = 3, mu = ",7.2, 'x', 10^-9, ")")),
                      y = expression(paste('Effective population size (N'[e], ')')))
      }
      
      p +
        scale_colour_manual(values = colour_palette) +
        scale_alpha_manual(values = c(0.1, 1), guide = FALSE) +
        scale_size_manual(values = c(0.4, 1)) +
        theme_bw() +
        theme(legend.position = 'bottom',
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 12)) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        # annotation_logticks() +
        facet_grid(clock ~ .)
      
    }) %>%
      set_names(value = df[['reference']])
  })
  
}

## Plot single PSMC - each sample to its own reference genome
plotPSMC_bestRef <- function(tbl, maxNe = NULL, filterSample = NULL, colour_palette){
  
  if(!is.null(maxNe)){
    tbl <- filter(tbl, Ne <= maxNe)
  }
  
  if(!is.null(filterSample)){
    tbl <- filter(tbl, ! sample %in% filterSample)
  }
  
  ## Split by mutation + generation combo
  tbl <- tbl %>%
    ungroup() %>%
    group_by(mu, gen) %>%
    group_split()
  
  nms <- tbl %>%
    map(~{
      .x %>%
        select(mu, gen) %>%
        distinct() %>%
        unite(col = 'nm', mu:gen, sep = '_') %>%
        extract2('nm')
    }) %>%
    unlist()
  
  names(tbl) <- nms
  
  map(tbl, ~{
    
    g <- unique(.x[['gen']])
    
    p <- ggplot(data = .x,
                aes(x = YearsAgo, y = Ne,
                    colour = sample,
                    alpha = lt == 'Standard',
                    group = interaction(sample, bootstrap))) +
      geom_step() +
      scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
    
    if(g == 10){
      p <- p + labs(x = expression(paste("Years (g = 10, mu = ",1.25, 'x', 10^-8, ")")),
                    y = expression(paste('Effective population size (N'[e], ')')))
    } else {
      p <- p + labs(x = expression(paste("Years (g = 3, mu = ",7.2, 'x', 10^-9, ")")),
                    y = expression(paste('Effective population size (N'[e], ')')))
    }
    
    p +
      scale_colour_manual(values = colour_palette) +
      scale_alpha_manual(values = c(0.1, 1), guide = FALSE) +
      scale_size_manual(values = c(0.4, 1)) +
      theme_bw() +
      theme(strip.text.x = element_text(size = 13),
            strip.text.y = element_text(size = 13),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16),
            legend.position = 'bottom', 
            legend.title = element_blank(),
            legend.text = element_text(face = 'italic', size = 16)) +
      guides(colour = guide_legend(override.aes = list(size = 5))) +
      # annotation_logticks() +
      facet_grid(clock ~ .)
  })
  
}

## Single PSMC plot - one time interval
plotPSMC <- function(tbl, g, c, maxNe = NULL, filterSample = NULL, colour_palette){
  
  tbl <- tbl %>%
    filter(gen == g) %>%
    filter(clock == c)
  
  if(!is.null(maxNe)){
    tbl <- filter(tbl, Ne <= maxNe)
  }
  
  if(!is.null(filterSample)){
    tbl <- filter(tbl, ! sample %in% filterSample)
  }
  
  p <- ggplot(data = tbl,
              aes(x = YearsAgo, y = Ne,
                  colour = sample,
                  alpha = lt == 'Standard',
                  group = interaction(sample, bootstrap))) +
    geom_step() +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  
  if(g == 10){
    p <- p + labs(x = expression(paste("Years (g = 10, mu = ",1.25, 'x', 10^-8, ")")),
                  y = expression(paste('Effective population size (N'[e], ')')))
  } else {
    p <- p + labs(x = expression(paste("Years (g = 3, mu = ",7.2, 'x', 10^-9, ")")),
                  y = expression(paste('Effective population size (N'[e], ')')))
  }
  
  p +
    scale_colour_manual(values = colour_palette) +
    scale_alpha_manual(values = c(0.1, 1), guide = FALSE) +
    scale_size_manual(values = c(0.4, 1)) +
    theme_bw() +
    theme(legend.justification = c(1, 0), 
          legend.position = c(1, 0),
          legend.box.margin = margin(c(25,25,25,25)),
          legend.title = element_blank(),
          legend.text = element_text(face = 'italic', size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16)) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    annotation_logticks()
}

## Effect of scaling values - facet by samples
plotPSMC_scalingEffect <- function(tbl, c, colour_palette, filterSample = NULL){
  
  tbl <- filter(tbl, clock == c)
  
  if(!is.null(filterSample)){
    tbl <- filter(tbl, ! sample %in% filterSample)
  }
  
  tbl <- mutate(tbl, gen = as.character(gen))
  
  p <- ggplot(data = tbl,
              aes(x = YearsAgo, y = Ne,
                  colour = gen,
                  alpha = lt == 'Standard',
                  group = interaction(sample, bootstrap, gen))) +
    geom_step() +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + 
    labs(x = 'Years',
         y = expression(paste('Effective population size (N'[e], ')'))) +
    scale_colour_manual(values = colour_palette) +
    scale_alpha_manual(values = c(0.1, 1), guide = FALSE) +
    scale_size_manual(values = c(0.4, 1)) +
    theme_bw() +
    theme(
      legend.position = 'none',
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16)) +
    # guides(colour = guide_legend(override.aes = list(size = 3))) +
    annotation_logticks() +
    facet_wrap(sample ~ .)
}

## Time segment patterning effect
plotPSMC_clockEffect <- function(tbl, g, colour_palette, filterSample = NULL){
  
  tbl <- filter(tbl, gen != g)
  
  if(!is.null(filterSample)){
    tbl <- filter(tbl, ! sample %in% filterSample)
  }
  
  # tbl <- mutate(tbl, gen = as.character(gen))
  
  p <- ggplot(data = tbl,
              aes(x = YearsAgo, y = Ne,
                  colour = clock,
                  alpha = lt == 'Standard',
                  group = interaction(sample, bootstrap, clock))) +
    geom_step() +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) + 
    labs(x = 'Years',
         y = expression(paste('Effective population size (N'[e], ')'))) +
    scale_colour_manual(values = colour_palette) +
    scale_alpha_manual(values = c(0.1, 1), guide = FALSE) +
    scale_size_manual(values = c(0.4, 1)) +
    theme_bw() +
    theme(
      legend.position = 'bottom',
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    annotation_logticks() +
    facet_wrap(sample ~ .)
  
}
