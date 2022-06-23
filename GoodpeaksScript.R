# install.packages(readr, tidyverse, stringr, reshape2)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("rtracklayer","BSgenome.Hsapiens.UCSC.hg38", "ChIPpeakAnno", "TxDb.Hsapiens.UCSC.hg38.knownGene"))

library(readr)
library(tidyverse)
library(rtracklayer)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(reshape2)

# this is where output goes
#setwd("~/Desktop/")

#####
#functions
BedInputToGRanges <- function(file) {
  df <- read_delim(file, delim = "\t", escape_double = FALSE, col_names = FALSE, skip = 1,
                   trim_ws = TRUE)
  df <- df[,1:6]
  colnames(df) <- c("chr", "start", "end", "name", "score", "strand")
  gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
  gr <- keepStandardChromosomes(gr, pruning.mode="tidy")
  seqlevels(gr, pruning.mode="coarse") <- seqlevels(gr)[seqlevels(gr) != "chrM"]
  return(gr)
}
BedGraphInputToRle <- function(file) {
  gr <- import.bedGraph(file, genome = "hg38")
  gr <- keepStandardChromosomes(gr, pruning.mode="coarse")
  seqlevels(gr, pruning.mode="coarse") <- seqlevels(gr)[seqlevels(gr) != "chrM"]
  coverage <- coverage(gr, weight = "score")
  return(coverage)
}
binnedSum <- function(bins, numvar, mcolname){
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewSums(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}
binnedMax <- function(bins, numvar, mcolname){
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewMaxs(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}
recenterPeaks <- function(consensusPeaks, DBbedgraph, VEHbedgraph, pad = 200) {
  pileup <- DBbedgraph+VEHbedgraph
  bins_per_chrom <- split(ranges(consensusPeaks), seqnames(consensusPeaks))
  summit_list <- lapply(names(pileup),
                        function(seqname) {
                          views <- Views(pileup[[seqname]],
                                         bins_per_chrom[[seqname]])
                          viewWhichMaxs(views)
                        })
  summits <- unsplit(summit_list, as.factor(seqnames(consensusPeaks)))
  GrangesDF <- data.frame(seqnames = seqnames(consensusPeaks),
                          start = summits - pad,
                          end = summits + pad)
  centeredConsensusPeaks <- makeGRangesFromDataFrame(GrangesDF, ignore.strand = TRUE)
  centeredConsensusPeaks
}
GCanalysis <- function(peaks, pad = 50) {
  regions <- resize(peaks, 2*pad+1, fix = "center")
  seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, regions)
  GC <- letterFrequency(seq, "GC", as.prob = TRUE)[,1]
  GC
}

#####
#test parameter
# DBpeakcalls <- BedFiles$MOLM13_DB
# VEHpeakcalls <- BedFiles$MOLM13_VEH
# DBbedgraph <- BedGraphFiles$MOLM13_DB
# VEHbedgraph <- BedGraphFiles$MOLM13_VEH
# Prefix <- "MOLM13"
# minPeakScore <- 50
# FCcutoff <- 3
# minSummit <- 5

#main function 
GoodpeaksFromBedgraph <- function(DBpeakcalls, #GRanges object from bed file
                                  VEHpeakcalls, #GRanges object from bed file
                                  DBbedgraph, #Rle object from bedgraph
                                  VEHbedgraph, #Rle object from bedgraph
                                  FilePaths = FALSE, #If TRUE will import files from path provided instead
                                  Prefix = "CellLine", #prefix for output files
                                  minPeakScore = 50,
                                  FCcutoff = 3,
                                  minSummit = 5,
                                  pad = 200, # pad around summit for re-centered peaks; 0 = skip recenter
                                  GCpad = 100 #pad around summit for GC analysis; 0 = skip
){
  if(FilePaths == TRUE){
    DBpeakcalls <- BedInputToGRanges(DBpeakcalls)
    VEHpeakcalls <- BedInputToGRanges(DBpeakcalls)
    DBbedgraph <- BedGraphInputToRle(DBbedgraph)
    VEHbedgraph <- BedGraphInputToRle(VEHbedgraph)
  }
  stopifnot(is(DBpeakcalls, "GRanges"))
  stopifnot(is(VEHpeakcalls, "GRanges"))
  stopifnot(is(DBbedgraph, "RleList"))
  stopifnot(is(VEHbedgraph, "RleList"))
  print(paste0("loaded ", length(DBpeakcalls), " DB peak calls and ", length(VEHpeakcalls), " VEH peak calls"))
  #combine
  MergedPeaks <- c(DBpeakcalls, VEHpeakcalls)
  N1 <- length(MergedPeaks)
  #drop low score peaks
  MergedPeaks <- MergedPeaks[MergedPeaks$score > minPeakScore]
  print(paste0("minPeakScore Filtered ", N1 - length(MergedPeaks), " peak calls, ", length(MergedPeaks), " peak calls remaining"))
  #create consensus peaks and recenter
  MergedPeaks <- GenomicRanges::reduce(MergedPeaks)
  if(pad > 0){MergedPeaks <- recenterPeaks(MergedPeaks, DBbedgraph, VEHbedgraph, pad = pad)}
  print(paste0(length(MergedPeaks), " consensus peaks found"))
  N1 <- length(MergedPeaks)
  #gather data
  MergedPeaks <- MergedPeaks %>% 
    binnedAverage(VEHbedgraph, "VEH_mean") %>%
    binnedAverage(DBbedgraph, "DB_mean") %>% 
    binnedSum(VEHbedgraph, "VEH_sum") %>% 
    binnedSum(DBbedgraph, "DB_sum") %>% 
    binnedMax(VEHbedgraph, "VEH_max") %>%
    binnedMax(DBbedgraph, "DB_max")
  #drop low summit peaks
  MergedPeaks <- MergedPeaks[pmax(MergedPeaks$VEH_max,MergedPeaks$DB_max) > minSummit]
  print(paste0("minSummit Filtered ", N1 - length(MergedPeaks), " consensus peaks, ", length(MergedPeaks), " peaks remaining"))
  #calculate change
  MergedPeaks$Log2FC <- log2(pmax(MergedPeaks$DB_sum,1) / pmax(MergedPeaks$VEH_sum,1))
  MergedPeaks$Change <- "Unchanged"
  MergedPeaks$Change[MergedPeaks$Log2FC > log2(FCcutoff)] <- "Gained"
  MergedPeaks$Change[MergedPeaks$Log2FC < -log2(FCcutoff)] <- "Lost"
  print(table(MergedPeaks$Change))
  if(GCpad > 0){
    MergedPeaks$GC <- GCanalysis(MergedPeaks, GCpad)
    barplot(c(mean(MergedPeaks$GC[MergedPeaks$Change == "Gained"]),
              mean(MergedPeaks$GC[MergedPeaks$Change == "Unchanged"]),
              mean(MergedPeaks$GC[MergedPeaks$Change == "Lost"])),
            main = Prefix,
            xlab = paste0("minPeakScore=", minPeakScore,
                          " FCcutoff=", FCcutoff,
                          " minSummit=", minSummit,
                          " pad=", pad,
                          " GCpad=", GCpad),
            ylab = "GC content",
            names.arg=c(paste("Gained", table(MergedPeaks$Change)[1]),
                        paste("Unchanged", table(MergedPeaks$Change)[3]),
                        paste("Lost", table(MergedPeaks$Change)[2])))
  }
  #write bed files
  export.bed(MergedPeaks[MergedPeaks$Change == "Gained"], paste0(Prefix,"_gained.bed"))
  export.bed(MergedPeaks[MergedPeaks$Change == "Lost"], paste0(Prefix,"_lost.bed"))
  export.bed(MergedPeaks[MergedPeaks$Change == "Unchanged"], paste0(Prefix,"_unchanged.bed"))
  write.csv(MergedPeaks, paste0(Prefix,"_Goodpeaks.csv"))
  return(MergedPeaks)
}

#read files -- skip if reading data in next code chunk
TEST1_peaks <- BedInputToGRanges("~/Desktop/CUTnTAG/SUMMARY FILES for FIGS/MOLM13/1st/MACS2/sig2/IgG-VEH-Control_J1158.bedgraph_peakssig2.bed")
TEST1_bedgraph <- BedGraphInputToRle("~/Desktop/CUTnTAG/SUMMARY FILES for FIGS/MOLM13/1st/BedGraph/IgG-VEH-Control_J1158.bedgraph")
TEST2_peaks <- BedInputToGRanges("~/Desktop/CUTnTAG/SUMMARY FILES for FIGS/MOLM13/1st/MACS2/sig2/IgG-DB-Control_J1158.bedgraph_peakssig2.bed")
TEST2_bedgraph <- BedGraphInputToRle("~/Desktop/CUTnTAG/SUMMARY FILES for FIGS/MOLM13/1st/BedGraph/IgG-DB-Control_J1158.bedgraph")
TEST3_peaks <- BedInputToGRanges("~/Desktop/CUTnTAG/SUMMARY FILES for FIGS/MOLM13/1st/MACS2/p0_05/IgG-DB-Control_J1158_p0_05.bed")
TEST3_bedgraph <- BedGraphInputToRle("~/Desktop/CUTnTAG/SUMMARY FILES for FIGS/MOLM13/1st/Bedgraph/IgG-DB-Control_J1158.bedgraph")

#run goodpeaks function

#LOW-stringe: minPeak=5, FC=3,minSum=2
#HIGH-stringe: minPeak=7.5, FC=4,minSum=3
#ForChIP= minPeak=50, FC=3,minSum=15 (03/08/22)
Goodpeaks1 <- GoodpeaksFromBedgraph(DBpeakcalls = TEST2_peaks, 
                                       VEHpeakcalls = TEST1_peaks, 
                                       DBbedgraph = TEST2_bedgraph, 
                                       VEHbedgraph = TEST1_bedgraph, 
                                       Prefix = "IgGVEHvsIgGDBpeaksCHIP_GP",
                                       minPeakScore = 7.5, 
                                       FCcutoff = 4, 
                                       minSummit = 3,
                                       pad = 200,
                                       GCpad = 100)
Goodpeaks2 <- GoodpeaksFromBedgraph(DBpeakcalls = TEST2_peaks, 
                                       VEHpeakcalls = TEST1_peaks, 
                                       DBbedgraph = TEST2_bedgraph, 
                                       VEHbedgraph = TEST1_bedgraph, 
                                       Prefix = "VEHvsDB2313XDaunpeaks_HIGH_STRINGE",
                                       minPeakScore = 7.5, 
                                       FCcutoff = 4, 
                                       minSummit = 2.5,
                                       pad = 200,
                                       GCpad = 50)


