#180413_findCompartmentBoundaries.R

#Using PC1 from matrix2compartments.pl to determine the location of boundaries between A and B compartments

####Load Libraries####
library(tidyverse)
library(corrplot)
library(psych)
library(Hmisc)
library(GGally)
library(GenomicRanges)
library(rtracklayer)


#compartmentBoundaries function

compartmentBoundaries = function(compTsvFile, outDir) {
  
  #extract chromosome number from compartment file name
  pattern = "chr[0123456789XY]{1,2}"
  match = regexpr(pattern, basename(compTsvFile))
  chr = regmatches(basename(compTsvFile), match)
  
  compTsv = read_tsv(file = compTsvFile)
  
  #Calculate which bins flank sign changes
  compBounds = tibble(Bin1Start = compTsv$start[1:(length(compTsv$start)-1)],
                      Bin1End = compTsv$end[1:(length(compTsv$end)-1)],
                      Bin2Start = compTsv$start[2:length(compTsv$start)],
                      Bin2End = compTsv$end[2:length(compTsv$end)],
                      Boundary = (compTsv$eigen1[1:(length(compTsv$eigen1)-1)] * compTsv$eigen1[2:length(compTsv$eigen1)]) <= 0)
  
  #Taking entire range (Bin1Start to Bin2End) of all pairs which are TRUE in boundary column. Making genomic ranges of this. So each boundary in this GR will be 2 bins in size.
  
  compBoundBins = na.omit(compBounds[compBounds$Boundary == TRUE,])
  compBoundGR = GRanges(seqnames = c(chr), ranges = IRanges(start = compBoundBins$Bin1Start, end = compBoundBins$Bin2End))
  
  #Now narrowing to 100kb region around center of each boundary range (can adjust this as needed by changing the value of the width variable)
  compBoundGRNarrow = resize(compBoundGR, width = 100000, fix = "center")
  
  #Saving as .bed
  export(compBoundGR, con = file.path(outDir, paste(chr, "_CompartmentBoundary2bins.bed", sep = "")), format = "BED")
  export(compBoundGRNarrow, con = file.path(outDir, paste(chr, "_CompartmentBoundary100kb.bed", sep = "")), format = "BED")
  #Saving as .rda
  save(compBoundGR, file = file.path(outDir, paste(chr, "_CompartmentBoundaryGR2bins.rda", sep = "")))
  save(compBoundGRNarrow, file = file.path(outDir, paste(chr, "_CompartmentBoundaryGR100kb.rda", sep = "")))
}
