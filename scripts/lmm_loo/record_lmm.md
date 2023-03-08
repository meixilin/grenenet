# run lmm models

```bash
# Tue Mar  7 16:56:17 PST 2023
sbatch /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/scripts/lmm_loo/step1_run_lmm_wrapper.sh # JOB: 177292
# check job completion (All done)
cd /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/logs/lmm_loo/step1_run_lmm
grep 'JOB ID' *.out.txt | grep 'Done' | cut -d '.' -f 1 | cut -d '_' -f 5 | sort -V > done_jobs.txt
```

