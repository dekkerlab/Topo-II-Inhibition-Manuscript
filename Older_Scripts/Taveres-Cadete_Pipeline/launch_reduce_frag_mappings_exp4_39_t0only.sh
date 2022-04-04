INPUT_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/alignments/exp4-39/done

SCRIPT_DIR=/nl/umw_job_dekker/users/eh37w/Topo-Inhib/MC3C/scripts


for f in $INPUT_DIR/*sorted_stricter.fragments.txt; do
	echo "$f"
	bsub -q short -W 4:00 -e /home/eh37w/lsf_jobs/LSB_%J.err -o /home/eh37w/lsf_jobs/LSB_%J.log -n 8 -R "select[ib]" -R "rusage[mem=2000]" -R "span[hosts=1]" -N -u erica.hildebrand@umassmed.edu "Rscript $SCRIPT_DIR/reduce_frag_mappings.R $f"
done
