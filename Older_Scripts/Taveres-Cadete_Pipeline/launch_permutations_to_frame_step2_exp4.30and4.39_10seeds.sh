for i in $(seq 20201221 20201231); do bsub -R rusage[mem=8000] -n 2 -R "span[hosts=1]" -W 240 -q short -e /home/eh37w/lsf_jobs/LSB_%J.err -o /home/eh37w/lsf_jobs/LSB_%J.log -N -u erica.hildebrand@umassmed.edu "Rscript interactions_to_usable_frame_for_permutations_EMH_R1R2_10seeds.R $i"; done

