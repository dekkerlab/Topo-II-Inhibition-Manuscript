
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

fragment_file_path = "/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/reference/hg19__GATC.txt" #changed to new path

bam_reader = HTSeq.BAM_Reader( sys.argv[1] ) #uses HTSeq.BAM_Reader to read in alignments - digested MC-3C to undigested genome
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



#making dict of all of the DpnII digested fragments in genome, will compare the alignments to these fragments. 

fragment_ranges = {}
for i in range(1, 23): #first adding chromosomes - this will depend on the genome that is used
    fragment_ranges[ "chr"+str(i)] = []

fragment_ranges["chrX"] = []
fragment_ranges["chrM"] = []
fragment_ranges["chrY"] = []

for fragment_range in HTSeq.BED_Reader( fragment_file_path ): #next adding start and end of each fragment iv stands for interval - part of bed reader output (https://htseq.readthedocs.io/en/master/otherparsers.html#bed-reader)

    fragment_range.iv.start -= 1 #subtract 1 from the start and the end of each interval - probably due to 1 vs 0 based intervals? Python vs R?
    fragment_range.iv.end -= 1

    fragment_ranges[fragment_range.iv.chrom].append( fragment_range ) #Add fragment range to correct chromosome 

 
frag_pacbio_dict = {} #initialize dict for fragments that match MC-3C alignments

for a in bam_reader: #from earlier HTSeq_BAMReader call - this is the alignments
    read_id = a.read.name.split("/")[1] #BAM read.name contains moviename/holenumber/ccs/subreadnumber (https://pacbiofileformats.readthedocs.io/en/5.0/BAM.html)
    #holenumber should be unique for each ccs in a given movie - so this is used as the read_id
    if len(a.read.name.split("/")) == 4: #determine if there are subreads - if not, will only be length 3
        digest_id = a.read.name.split("/")[3] #digest_id number from bam read.name of this alignment - where within css it came from (order)
    else:
        digest_id = "0"
    if (read_id in frag_pacbio_dict) == False: #check if this read_id has already been added - if not, add it. 
        frag_pacbio_dict[read_id] = OrderedDict()
    if (digest_id in frag_pacbio_dict[read_id]) == False: #check if this digest id has already been added, if not, add it
        frag_pacbio_dict[read_id][digest_id] = []
    if (a.aligned) & (a.aQual == 60): #only add an aligned fragment if alignment is high quality
        for frag in fragment_ranges[a.iv.chrom]: #Go through all genome fragments in the matching chromosome to see if they overlap with alignment
            if a.iv.overlaps(frag.iv): #if they overlap, calculate how much of fragment overlaps with alignment
                overlap_length = float((min([a.iv.end, frag.iv.end]) - max([a.iv.start, frag.iv.start]))) / (frag.iv.end - frag.iv.start)
                if overlap_length >= 0.8: #if it is over 80% overlap, add this fragment to the frag_pacbio_dict as the alignment
                    frag_pacbio_dict[read_id][digest_id].append((frag, a.iv.strand))

pickle.dump( frag_pacbio_dict, open(sys.argv[1].replace(".bam", "_stricter.pickle"), "wb")) #save dict


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

