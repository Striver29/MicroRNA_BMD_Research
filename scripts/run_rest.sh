#!/bin/bash
# Step 8 first: its by-product is the per-family PATIENT-label-shuffle p for the
# vs-control contrasts. Step 9 showed fgsea was overconfident by ~8 orders of magnitude
# for BMD-vs-DMD, so we must find out whether the vs-control miR-29 result survives the
# same test. That now decides what the paper says.
export PATH="/opt/homebrew/bin:$PATH"
cd /Users/abdulazizalshujaa/Documents/MicroRNA_BMD_Research
echo "[$(date +%H:%M)] START 08 (500 shuffles, 4 workers)"
N_SHUFFLE=500 NPERM_NULL=10000 N_CORE=4 Rscript scripts/08_shuffle_baseline.R > results/tables/08_shuffle_log.txt 2>&1
echo "[$(date +%H:%M)] END 08 exit=$?"
echo "[$(date +%H:%M)] START 10b (4 workers)"
N_REP=100 NPERM=10000 N_CORE=4 Rscript scripts/10b_detection_limit.R > results/tables/10b_log.txt 2>&1
echo "[$(date +%H:%M)] END 10b exit=$?"
echo "[$(date +%H:%M)] START 09 classification (cheap, needs step 8 output)"
N_SHUFFLE=1 N_CORE=2 Rscript scripts/09_bmd_vs_dmd.R > results/tables/09_classify_log.txt 2>&1
echo "[$(date +%H:%M)] END 09 classification exit=$?"
echo "[$(date +%H:%M)] ALL DONE"
