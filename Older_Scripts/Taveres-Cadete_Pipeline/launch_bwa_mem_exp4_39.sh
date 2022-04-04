module load bwa/0.7.12
module load samtools/1.3

OUTPUT_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/alignments/exp4-39

DIGESTED_READS_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/digested/exp4-39

REFERENCE_DIR=/nl/umw_job_dekker/cshare/reference/bwa/hg19

for f in $DIGESTED_READS_DIR/*fastq; do
	echo "$f"
	dataset="$(basename -- $f .digested.fastq)"
	echo "$dataset"
	echo "$dataset.bam"
	bsub -q short -n 8 -R "rusage[mem=32000]" -R "span[hosts=1]" -W 240 \
		"bwa mem -t 8 $REFERENCE_DIR/hg19.fa $f | samtools view -bS - > $OUTPUT_DIR/$dataset.bam 
		samtools sort -n $OUTPUT_DIR/$dataset.bam > $OUTPUT_DIR/$dataset.sorted.bam 
		samtools flagstat $OUTPUT_DIR/$dataset.sorted.bam > $OUTPUT_DIR/$dataset.flagstat_out"
done