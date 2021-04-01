module load fastqc/0.11.5

FASTQ_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/fastq/exp4-39

DIGESTED_READS_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/digested/exp4-39

#probably don't need as many cores/as much memory next time
bsub -q short -n 8 -R "rusage[mem=32000]" -R "span[hosts=1]" -W 240 \
		"fastqc -o $FASTQ_DIR -f fastq $FASTQ_DIR/*.fastq.gz"

bsub -q short -n 8 -R "rusage[mem=32000]" -R "span[hosts=1]" -W 240 \
		"fastqc -o $DIGESTED_READS_DIR -f fastq $DIGESTED_READS_DIR/*.fastq"