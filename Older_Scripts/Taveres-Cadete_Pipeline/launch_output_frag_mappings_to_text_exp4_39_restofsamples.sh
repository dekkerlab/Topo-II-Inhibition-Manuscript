INPUT_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/alignments/exp4-39

SCRIPT_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/scripts


for f in $INPUT_DIR/*.sorted_stricter.pickle; do
	echo "$f"
	bsub -q short -W 4:00 -e /home/eh37w/lsf_jobs/LSB_%J.err -o /home/eh37w/lsf_jobs/LSB_%J.log -n 2 -R "select[ib]" -R "rusage[mem=2000]" -R "span[hosts=1]" -N -u erica.hildebrand@umassmed.edu "python $SCRIPT_DIR/output_frag_mappings_to_text.py $f"
done