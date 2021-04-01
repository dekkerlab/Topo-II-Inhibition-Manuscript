# Script to implement some permutation/randomisation of interactions as a proof of principle
rm(list=ls())
gc()

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(rtracklayer)
library(igraph)

# Load the interactions
intDataDir = "/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/walks"
load(file.path(intDataDir, "201221_exp4.30and4.39_AllSamples_frames_for_plotting_stricter_R1R2combined.RData"))

seed <- as.integer(commandArgs(trailingOnly = TRUE)[[1]])
set.seed(seed)

# I first thought of two approaches to do the permutation:
# 1) Keep the fragments, randomise interactions between fragments from the same walk, keeping walk length
# 2) Keep the fragments, but  create new interactions that can span fragments from different walks. Keep the distribution of walk lengths
# Let's try to get them done.


# 1) Keep the fragments, randomise interactions between fragments from the same walk, keeping walk length
# Split direct_ints by the pacbio_frag_ID
direct_ints_split <- split(direct_ints_new_R1R2, direct_ints_new_R1R2$pacbio_frag_ID)

# Get the fragment information, avoiding repeating frags that are only visited once
fragment_rangeslist <- mclapply(direct_ints_split,
                                function(x) {
                                    n_rows <- dim(x)[[1]]
                                    ranges <- append(with(x,
                                                          GRanges(seqnames = V1,
                                                                  ranges = IRanges(start = as.numeric(V2),
                                                                                   end = as.numeric(V3)),
                                                                  strand = "*",
                                                                  pacbio_frag_ID = pacbio_frag_ID,
                                                                  walk_length = steps,
                                                                  walk_class = class,
                                                                  full_dataset_ID = full_dataset_ID)),
                                                     with(x[n_rows,],
                                                          GRanges(seqnames = V4,
                                                                  ranges = IRanges(start = as.numeric(V5),
                                                                                   end = as.numeric(V6)),
                                                                  strand = "*",
                                                                  pacbio_frag_ID = pacbio_frag_ID,
                                                                  walk_length = steps,
                                                                  walk_class = class,
                                                                  full_dataset_ID = full_dataset_ID)))
                                    unique(ranges)
                                },
                                mc.cores = getOption("mc.cores", 64L))

# Change the order of the ranges and output as an interaction data.frame
permutated_direct_ints <- mclapply(fragment_rangeslist,
                                   function(x) {
                                       reshufled_ranges <- sample(x, replace = FALSE)
                                       to_return <- cbind(as.data.frame(reshufled_ranges[1:(length(reshufled_ranges)-1)])[, 1:3],
                                                          as.data.frame(reshufled_ranges[2:(length(reshufled_ranges))])[, c(1:3, 6:9)])
                                       return(to_return)

},
				   mc.cores = getOption("mc.cores", 64L))

permutated_direct_ints <- do.call(rbind, permutated_direct_ints)

permutated_direct_ints$seed <- seed

# Save the permutated interactions in a file for later retrieval and analysis
outDataDir = "/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/permutated_direct_ints"
save(permutated_direct_ints, file = file.path(outDataDir, paste0("permutated_direct_ints_4.30and4.39_", seed, ".RData")))

