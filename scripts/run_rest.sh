#!/bin/bash
# Runs 10b then 08 at 4 workers (M1 Air has 4 performance cores; the 4 efficiency
# cores add heat for little speed on a fanless chassis).
# Waits on the OUTPUT FILE, not pgrep -- pgrep -f matches this script's own text.
export PATH="/opt/homebrew/bin:$PATH"
cd /Users/abdulazizalshujaa/Documents/MicroRNA_BMD_Research
until [ -f results/rds/bmd_vs_dmd.rds ]; do sleep 30; done
echo "[$(date +%H:%M)] step 9 output present"
echo "[$(date +%H:%M)] START 10b (4 workers)"
N_REP=100 NPERM=10000 N_CORE=4 Rscript scripts/10b_detection_limit.R > results/tables/10b_log.txt 2>&1
echo "[$(date +%H:%M)] END 10b exit=$?"
echo "[$(date +%H:%M)] START 08 (4 workers)"
N_SHUFFLE=500 NPERM_NULL=10000 N_CORE=4 Rscript scripts/08_shuffle_baseline.R > results/tables/08_shuffle_log.txt 2>&1
echo "[$(date +%H:%M)] END 08 exit=$?"
echo "[$(date +%H:%M)] ALL DONE"
