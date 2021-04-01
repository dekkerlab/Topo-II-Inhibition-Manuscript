OUTPUT_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/digested/exp4-39

FASTQ_READS_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/fastq/exp4-39

for f in $FASTQ_READS_DIR/*fastq.gz; do
	echo "$f"
	dataset="$(basename -- $f .fastq.gz)"
	echo "$dataset"
	echo "$dataset.digested.fastq"
	bsub -q short -n 4 -R "rusage[mem=8000]" -R "span[hosts=1]" -W 240 \
		"python digest_roi_python3.py $f $OUTPUT_DIR/$dataset.digested.fastq"
done