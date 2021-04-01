
import sys
import HTSeq
import pickle

frag_pacbio_dict = pickle.load( open(sys.argv[1], "rb") )

output_file = open(sys.argv[1].replace(".pickle", ".fragments.txt"), "w")

for pacbio_read in frag_pacbio_dict.keys():
    
    for subread in frag_pacbio_dict[pacbio_read].keys():

        if len( frag_pacbio_dict[pacbio_read][subread] ) == 0:

            output_file.write( "\t".join([pacbio_read, subread, "*", "*", "*", "*"]) + "\n" )

        else:

            for interval in frag_pacbio_dict[pacbio_read][subread]:
                
                output_file.write( "\t".join([pacbio_read, subread, interval[0].iv.chrom, str(interval[0].iv.start), str(interval[0].iv.end), interval[1]]) + "\n")


output_file.close()



