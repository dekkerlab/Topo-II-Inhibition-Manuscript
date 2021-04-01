
import sys
from Bio import Restriction
from Bio.Seq import Seq
from Bio import SeqIO
import gzip
from mimetypes import guess_type
from functools import partial

print( sys.argv )

file_paths = sys.argv[1]

output_file_path = sys.argv[2]

encoding = guess_type(file_paths)[1] #uses file extention
_open = partial(gzip.open, mode = 'rt') if encoding == 'gzip' else open

output_handle = open(output_file_path, "w")

handle = _open(file_paths)

records = list(SeqIO.parse(handle, "fastq"))

handle.close()

for record in records:
        
    restriction_site_len = len( Restriction.DpnII.site )
        
    cut_sites = Restriction.DpnII.search(record.seq)
        
    if len(cut_sites) == 0:

        output_handle.write(record.format("fastq"))

    else :

        subread_counter = 0

        for i in range(0, len(cut_sites)):
            
            if i == 0:
            
                subread = record[0 : cut_sites[i]-1]

            else:
            
                subread = record[cut_sites[i-1]-1 : cut_sites[i]-1]

            if len(subread.seq) <= 4:
                continue
           
            subread.id += "/0" + str(subread_counter)
            subread.name += "/0" + str(subread_counter)
            subread.description += "/0" + str(subread_counter)

            subread_counter += 1

            output_handle.write(subread.format("fastq"))
            
        subread = record[cut_sites[i]-1 :]

        subread.id += "/0" + str(subread_counter)
        subread.name += "/0" + str(subread_counter)
        subread.description += "/0" + str(subread_counter)

        if len(subread.seq) <= 4:
            continue

        output_handle.write(subread.format("fastq"))

        
output_handle.close()
