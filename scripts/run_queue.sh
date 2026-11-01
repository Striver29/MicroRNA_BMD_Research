#!/bin/bash
# Waits for step 8 to finish, then runs the remaining steps sequentially at full core count.
export PATH="/opt/homebrew/bin:$PATH"
cd /Users/abdulazizalshujaa/Documents/MicroRNA_BMD_Research
until grep -q "SHUFFLED-LABEL BASELINE" results/tables/08_shuffle_log.txt 2>/dev/null; do sleep 20; done
echo "[$(date +%H:%M)] step 8 finished"
for job in "11_utr_matched" "10b_detection_limit" "09_bmd_vs_dmd" "13_leave_one_out" "14_subsample_dmd"; do
  echo "[$(date +%H:%M)] starting $job"
  N_CORE=7 Rscript "scripts/${job}.R" > "results/tables/${job%%_*}_log.txt" 2>&1
  echo "[$(date +%H:%M)] $job exit=$?"
done
echo "[$(date +%H:%M)] QUEUE COMPLETE"
