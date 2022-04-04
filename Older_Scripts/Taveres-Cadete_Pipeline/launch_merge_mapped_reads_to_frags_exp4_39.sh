module load samtools/1.3

INPUT_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/alignments/exp4-39

SCRIPT_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/scripts


for f in $INPUT_DIR/*.sorted.bam; do
	echo "$f"
	bsub -q long -W 32:00 -e /home/eh37w/lsf_jobs/LSB_%J.err -o /home/eh37w/lsf_jobs/LSB_%J.log -n 2 -R "select[ib]" -R "rusage[mem=32000]" -R "span[hosts=1]" -N -u erica.hildebrand@umassmed.edu "python $SCRIPT_DIR/merge_mapped_reads_to_frags_stricter.py $f"
done