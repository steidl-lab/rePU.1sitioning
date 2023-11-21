suppressPackageStartupMessages({
  library(tidyverse)
  library(universalmotif)
  library(BSgenome.Hsapiens.UCSC.hg38)
})


GetMotifSequences <- function(
  gr,                          #GRanges of peaks
  motif,                       #motif query
  motifcenter,                 #central base in motif for finding distance to peak center and visualization purposes 
  flank = 25,                  #motif flank length
  return = "GRupdate",         #"GRupdate" adds central match info to input gr metadata and returns updated gr; "full" returns full results from scan_sequences; "seqs" returns central match motif seq + flanks
  genome = BSgenome.Hsapiens.UCSC.hg38, #genome object
  threshold = 0.45,            # pass to scan_sequences
  threshold.type = "logodds"   # pass to scan_sequences
  ) {
  
  #setup variables
  motifname <- motif@name
  motiflen <- ncol(motif)
  peakcenter <- 101 + flank + 10 
  
  #get seqs
  seqs <- getSeq(genome, gr + flank + 10)
  names(seqs) <- gr$PeakID #names needed to pull sequences
  
  #find motif match
  motifsfull <- scan_sequences(motif, seqs, RC = TRUE, threshold = threshold,  threshold.type = threshold.type) %>% 
    as_tibble() %>%
    mutate(center.dist = case_when(strand == "+" ~ abs(start + motifcenter - peakcenter),
                                   strand == "-" ~ abs(start - motifcenter - peakcenter))) %>%
    filter(center.dist < 100) 
  
  #select match closest to center and convert to GRanges
  motifscentral <- motifsfull %>%
    arrange(center.dist) %>%
    distinct(sequence, .keep_all = TRUE) %>%
    mutate(end = pmax(start , stop), #hacky setup to fix RC hits and get makeGRanges working
           start = pmin(start , stop)) %>%
    select(sequence, start, end, strand, score, center.dist) %>%
    makeGRangesFromDataFrame(seqnames.field = "sequence", keep.extra.columns = TRUE)
  
  #Add flanking sequences to motif match 
  motifseq <- motifscentral %>%
    promoters(upstream = flank - motifcenter + 2, downstream = flank + motiflen - motifcenter + 4) %>%
    getSeq(seqs, .)
  
  motifscentral$motifseq <- motifseq
  
  #Add results to input gr
  gr1 <- gr %>% 
    as_tibble() %>% 
    left_join(
      motifscentral %>%
        as_tibble() %>%
        select(seqnames, score, center.dist, motifseq) %>%
        rename_with( ~ paste(motifname, .x, sep = "_"), !seqnames)
      , by = c("PeakID" = "seqnames")) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  #return
  if (return == "GRupdate") {
    return(gr1)
  } else if(return == "full") {
    return(motifsfull)
  } else if(return == "seqs") {
    return(motifseq)
  }
  
}

#Plotting functions
#Note: The functions assume a 4bp central motif (GGAA). Only affects visualization
BaseFreqbyBase <- function(
  x,                    #DNAStringSet of equal lenths (output of GetMotifSequences(return = "seqs"))
  bases = DNA_BASES,    #Vector of letters or letter conbinations 
  plot = FALSE,         #return plot
  len = 54,             #x axis setting
  flank = 25            #x axis setting
  ) {
  if(var(width(x)) != 0){
    stop("unequal seq lenths")
  }
    
  #find letter frequency by base
  l <- list()
  for(i in 1:width(x)[1]) {
    l[[i]] <- x %>% 
      subseq(i, i) %>% 
      paste0(collapse = "") %>% 
      DNAString() %>% 
      letterFrequency(bases, as.prob = TRUE)
  }
  
  df <- bind_rows(l) %>%
    mutate(position = 1:width(x)[1]) %>%
    pivot_longer(cols = !position, names_to = "base", values_to = "frequency")
  
  if(plot == TRUE){
    df %>% 
      ggplot(aes(x=position, y=frequency, color=base, fill=base)) +
      geom_point() +
      geom_smooth(span = .15, alpha = 0.2) +
      theme_minimal() +
      scale_x_continuous(breaks=c(seq(0, flank, by=5), seq(flank+4, len, by=5)),
                         labels=c(seq(-flank,0, by=5), seq(0, flank, by=5)))
  } else {
    return(df)
  }
}

BaseFreqbyBasebyGroup <- function(
    gr,                         #GRanges (or df) containing sequences to plot
    group_by,                   #metadata column to group
    bases = "GC",               #letter or letter group to plot
    seqs_col = "motifseq",      #metadata column containing sequences
    plot = FALSE,               #return plot
    flank = 25,                 #x axis setting
    len = 54                    #x axis setting
    ) {
  if(length(bases) > 1){
    stop("bases must be vector of length 1")
  }
  df <- gr %>%
    as_tibble() %>%
    group_by(!!sym(group_by)) %>%
    group_modify(~ .x %>%
                select(contains(seqs_col)) %>%
                na.omit() %>%
                unlist() %>%
                DNAStringSet() %>%
                BaseFreqbyBase(bases = bases)) 
  
  if(plot == TRUE) {
    df %>% 
      ggplot(aes(x=position, y=frequency, color=!!sym(group_by), fill=!!sym(group_by))) +
      geom_point() +
      geom_smooth(span = .15, alpha = 0.2) +
      theme_minimal() +
      scale_x_continuous(breaks=c(seq(0, flank, by=5), seq(flank+4, len, by=5)),
                         labels=c(seq(-flank,0, by=5), seq(0, flank, by=5)))
  } else {
    return(df)
  }
  
}
