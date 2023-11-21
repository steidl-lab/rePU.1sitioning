#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("universalmotif")

library(tidyverse)
library(universalmotif)
library(BSgenome.Hsapiens.UCSC.hg38)

#load plotting functions
source("letterplotfunctions.R")

#Load or create GRanges object of centered peaks of 201bp. Peaks of varying lengths will not work.
DB2115_all_granges <- readRDS("MOLM13_all_022023_Granges.rds")
#DB2115_all_granges <- makeGRangesFromDataFrame(DB2115_all_granges, keep.extra.columns = TRUE)

#read query motif
pham_matrix <- as.numeric(read.delim("ShortPU1PWM.txt", sep = " ", header = F)[1,])
pham_matrix <- t(matrix(pham_matrix, ncol = 4, byrow = T))
pham_motif <- motif_rc(create_motif(pham_matrix, alphabet = "DNA", name = "PHAMfull", pseudocount = 0))
view_motifs(pham_motif)


#add motif search results to peak GRanges
gr <- GetMotifSequences(gr = DB2115_all_granges, motif = pham_motif, threshold = 0.9, threshold.type = "logodds", motifcenter = 4, flank = 25, return = "GRupdate", genome = BSgenome.Hsapiens.UCSC.hg38)

#get motif sequences
motifseqs <- GetMotifSequences(gr = DB2115_all_granges, motif = pham_motif, threshold = 0.9, threshold.type = "logodds", motifcenter = 4, flank = 25, return = "seqs", genome = BSgenome.Hsapiens.UCSC.hg38)

#plot frequency distributions

#All bases
BaseFreqbyBase(x = motifseqs, bases = DNA_BASES, plot = TRUE)
BaseFreqbyBase(x = motifseqs, bases = "GC", plot = TRUE)

BaseFreqbyBasebyGroup(gr, group_by = "Change", bases = "A", plot = TRUE)
BaseFreqbyBasebyGroup(gr, group_by = "Change", bases = "GC", plot = TRUE)

#From Granges to make 4 nucleotide plots
flank <- 25
len <- 54

gr %>%
  as_tibble() %>%
  filter(Change == "Gained") %>%
  dplyr::select(PHAMfull_motifseq) %>%
  na.omit() %>%
  unlist() %>%
  DNAStringSet() %>%
  BaseFreqbyBase() %>%
  ggplot(aes(x=position, y=frequency, color=base, fill=base)) +
  geom_point(alpha = 0.5, size = 1.0) +
  geom_smooth(size =1.5, span = .15, alpha = 0.2) +
  theme_classic() +
  labs(x="Distance from PU.1 motif", y="Base Frequency ", title="Base Frequency around PU.1 motifs") +
  theme(axis.text = element_text(size = 13),
        axis.title = element_blank(),
        legend.position = "none") +
  scale_x_continuous(breaks=c(seq(0, flank, by=5), seq(flank+4, len, by=5)),
                     labels=c(seq(-flank,0, by=5), seq(0, flank, by=5)))

#Create new motif
motif_lost <- DB2115_all_granges %>%
  as_tibble() %>%
  filter(Change == "Lost") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  GetMotifSequences(motif = pham_motif, 
                    motifcenter = 2, 
                    flank = 12, #due to inflexibility in the code, flank cannot be too short
                    return = "seqs") %>%
  subseq(start = 7, end = 26) %>% # this will resize the sequences to compensate
  create_motif(alphabet = "DNA", name = "Lost") 


motif_gained <- DB2115_all_granges %>%
  as_tibble() %>%
  filter(Change == "Gained") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  GetMotifSequences(motif = pham_motif, 
                    motifcenter = 2, 
                    flank = 12, #due to inflexibility in the code, flank cannot be too short
                    return = "seqs") %>%
  subseq(start = 7, end = 26) %>% # this will resize the sequences to compensate
  create_motif(alphabet = "DNA", name = "Gained") 

#motif figure
view_motifs(motifs = c(motif_lost, motif_gained, pham_motif),use.type= "PPM", sort.positions = 1)

#save new motif
write_matrix(motif_lost, "motif_lost_PWM.txt", type = "PWM")

