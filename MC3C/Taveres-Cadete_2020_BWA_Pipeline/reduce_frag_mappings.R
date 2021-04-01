
library(GenomicRanges)

args = commandArgs(trailingOnly = TRUE)

pacbio_subread_frags <- read.table(args[1],
				   header = FALSE, stringsAsFactors = FALSE)

pacbio_subread_frags_split <- split(pacbio_subread_frags, pacbio_subread_frags$V1)

all_results_frames <- mclapply(pacbio_subread_frags_split, function(teste) {

	result_frame <- data.frame(seqnames = teste[1,3],
				   start = teste[1,4],
				   end = teste[1,5],
				   strand = teste[1, 6],
				   stringsAsFactors = FALSE)

	if (dim(teste)[1] == 1) {

		result_frame <- data.frame(fragment = teste[1,1], result_frame)

		return(result_frame)

	}

	for ( i in 2:dim(teste)[1]) {

		rm("evaluation_row",
		   "add_to_row",
		   "evaluation_range",
		   "add_to_range")

		evaluation_row <- teste[i,]

		if (evaluation_row$V4 == "*") {
			result_frame <- rbind(result_frame, data.frame(seqnames="*", start="*", end="*", strand="*"))
			next
		}

		add_to_row <- result_frame[dim(result_frame)[1],]

		if (add_to_row$seqnames == "*") {
			result_frame <- rbind(result_frame, data.frame(seqnames=evaluation_row[,3],
								       start=evaluation_row[,4],
								       end=evaluation_row[,5],
								       strand=evaluation_row[,6]))
			next
		}

		if (evaluation_row$V3 != add_to_row$seqnames) {
			result_frame <- rbind(result_frame, data.frame(seqnames=evaluation_row[,3],
								       start=evaluation_row[,4],
								       end=evaluation_row[,5],
								       strand=evaluation_row[,6]))
			next
		}

		evaluation_range <- with(evaluation_row,
					 GRanges(seqnames = V3,
						 ranges = IRanges(start = as.integer(V4),
								  end = as.integer(V5)),
						 strand = V6))

		add_to_range <- with(add_to_row,
				     GRanges(seqnames = seqnames,
					     ranges = IRanges(start = as.integer(start),
							      end = as.integer(end)),
					     strand = strand))

		reduced_ranges <-  reduce(c(evaluation_range, add_to_range))

		if (length(reduced_ranges) == 1) {

			result_frame <- result_frame[-dim(result_frame)[1], ]

			result_frame <- rbind(result_frame,
					      as.data.frame(reduced_ranges)[, c("seqnames", "start", "end", "strand")])

		} else {

			result_frame <- rbind(result_frame,
					      as.data.frame(evaluation_range)[, c("seqnames", "start", "end", "strand")])

		}

	}

	result_frame <- data.frame(fragment = teste[1,1], result_frame)

	return(result_frame)


}, mc.cores = getOption("mc.cores", 8)

)

#save(all_results_frames, file="all_results_frames_test.RData")


all_direct_ints <- mclapply(all_results_frames, function(fragment) {

			    if (dim(fragment)[1] > 1) {

				    return_frame <- data.frame(stringsAsFactors=FALSE)

				    for (i in 1:(dim(fragment)[1]-1)) {

					    if ((fragment[i, "seqnames"] == "*") & (fragment[i+1, "seqnames"] == "*")) next

					    return_frame <- rbind(return_frame,
								  data.frame(seqnames1 = fragment[i, "seqnames"],
									     start1 = fragment[i, "start"],
									     end1 = fragment[i, "end"],
									     seqnames2 = fragment[i+1, "seqnames"],
									     start2 = fragment[i+1, "start"],
									     end2 = fragment[i+1, "end"],
									     fragmentID = fragment[i, "fragment"],
									     stringsAsFactors=FALSE))

				    }

				    return(return_frame)
			    }
}, mc.cores = getOption("mc.cores", 8)
)

all_direct_ints <- do.call(rbind, all_direct_ints)

write.table(all_direct_ints,
	    row.names = FALSE, col.names = FALSE,
	    quote = FALSE, file = gsub("fragments", "interactions", args[1]),
	    sep = "\t")
