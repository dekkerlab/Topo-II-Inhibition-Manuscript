
import HTSeq
from collections import Counter
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import itertools
from matplotlib.backends.backend_pdf import PdfPages
import sys
from collections import OrderedDict
import numpy
import pickle

fragment_file_path = "/nl/umw_job_dekker/archive/cshare/genome/restrictionFragments/hg19/hg19__GATC.txt"

bam_reader = HTSeq.BAM_Reader( sys.argv[1] )
#
#alignment_read_names = []
#
#alignment_dict = {}
#
#alignment_quals = []
#
#for a in bam_reader:
#
#    if a.aligned:
#        alignment_read_names.append( a.read.name )
#
#        seq_name_split = a.read.name.split("/")
#
#        if alignment_dict.has_key( seq_name_split[1] ):
#            alignment_dict[ seq_name_split[1] ].append(a)
#        else:
#            alignment_dict[ seq_name_split[1] ] = [a]
#
#        alignment_quals.append(a.aQual)
#
#number_of_alignments = Counter(alignment_read_names)
#
#pp = PdfPages(sys.argv[1].replace(".bam", ".pdf"))
#
#plt.figure()
#plt.hist(number_of_alignments.values(), max(number_of_alignments.values())-1)
#pp.savefig()
#
#plt.figure()
#plt.hist(alignment_quals, max(alignment_quals)-1)
#pp.savefig()
#




fragment_ranges = {}
for i in range(1, 23):
    fragment_ranges[ "chr"+str(i)] = []

fragment_ranges["chrX"] = []
fragment_ranges["chrM"] = []
fragment_ranges["chrY"] = []

for fragment_range in HTSeq.BED_Reader( fragment_file_path ):

    fragment_range.iv.start -= 1
    fragment_range.iv.end -= 1

    fragment_ranges[fragment_range.iv.chrom].append( fragment_range )

 
frag_pacbio_dict = {}

for a in bam_reader:
    read_id = a.read.name.split("/")[1]
    if len(a.read.name.split("/")) == 4:
        digest_id = a.read.name.split("/")[3]
    else:
        digest_id = "0"
    if (read_id in frag_pacbio_dict) == False:
        frag_pacbio_dict[read_id] = OrderedDict()
    if (digest_id in frag_pacbio_dict[read_id]) == False:
        frag_pacbio_dict[read_id][digest_id] = []
    if (a.aligned) & (a.aQual == 60):
        for frag in fragment_ranges[a.iv.chrom]:
            if a.iv.overlaps(frag.iv):
                overlap_length = float((min([a.iv.end, frag.iv.end]) - max([a.iv.start, frag.iv.start]))) / (frag.iv.end - frag.iv.start)
                if overlap_length >= 0.8:
                    frag_pacbio_dict[read_id][digest_id].append((frag, a.iv.strand))

pickle.dump( frag_pacbio_dict, open(sys.argv[1].replace(".bam", "_stricter.pickle"), "wb")) 


#proportion_digests = []
#
#step_sizes = {}
#step_sizes_list = numpy.array([])
#
#for pacbio_read in frag_pacbio_dict.keys():
#    found_frags = 0
#    for digest_read in frag_pacbio_dict[pacbio_read]:
#        if len(frag_pacbio_dict[pacbio_read][digest_read]) > 0:
#            found_frags += 1
#    fraction_found_frags = float(found_frags) / len(frag_pacbio_dict[pacbio_read])
#    proportion_digests.append( fraction_found_frags )
#    if (len(frag_pacbio_dict[pacbio_read]) > 1) & (fraction_found_frags == 1):
#        midpoints = []
#        for digest_read in frag_pacbio_dict[pacbio_read]:
#           min_pos = 1e100
#           max_pos = 0
#           for bed_range in frag_pacbio_dict[pacbio_read][digest_read]:
#               if bed_range.iv.start < min_pos:
#                   min_pos = bed_range.iv.start
#               if bed_range.iv.end > max_pos:
#                   max_pos = bed_range.iv.end
#           midpoints.append( ((max_pos - min_pos) / 2.0) + min_pos )
#        step_sizes[pacbio_read] = abs(numpy.diff(midpoints))
#        step_sizes_list = numpy.concatenate( (step_sizes_list, abs(numpy.diff(midpoints))) )
#
#
#
#plt.figure()
#plt.hist( proportion_digests )
#pp.savefig()
#
#pp.close()
#

