#!/bin/bash
# Thermal watchdog for the overnight run on a fanless MacBook Air M1.
#
# What it can and cannot see: CPU die temperature needs `sudo powermetrics`, so we use the
# battery/chassis sensors instead (no sudo). Those lag the die by a few degrees but they track
# how hot the machine FEELS, which is the actual concern on a fanless chassis.
#
#   Temperature         battery pack, in 0.01 C
#   VirtualTemperature  Apple's internal hotspot estimate, in 0.01 C
#   pmset -g therm      reports CPU_Speed_Limit < 100 when macOS is throttling
#
# Always appends a row to results/tables/thermal_log.tsv.
# Prints to stdout ONLY on a state change, so it is quiet unless something matters.
# Kills the R jobs if the battery stays at/above CRIT for two consecutive reads.

cd /Users/abdulazizalshujaa/Documents/MicroRNA_BMD_Research
LOG=results/tables/thermal_log.tsv
WARN=40      # C - getting hot for a fanless chassis
CRIT=45      # C - Apple's battery operating ceiling; stop the job
INTERVAL=${INTERVAL:-300}

[ -f "$LOG" ] || printf "timestamp\tbatt_C\tvirt_C\tcpu_speed_limit\tr_procs\tstate\n" > "$LOG"

read_temps () {
  ioreg -r -n AppleSmartBattery -a 2>/dev/null > /tmp/_batt.plist
  python3 -I - <<'PY'
import plistlib
try:
    d = plistlib.load(open('/tmp/_batt.plist','rb'))[0]
    print(f"{d.get('Temperature',0)/100:.1f} {d.get('VirtualTemperature',0)/100:.1f}")
except Exception:
    print("0.0 0.0")
PY
}

prev_state=""; hot_streak=0
while true; do
  set -- $(read_temps); BATT=$1; VIRT=$2
  LIMIT=$(pmset -g therm 2>/dev/null | awk -F'= *' '/CPU_Speed_Limit/{print $2}')
  [ -z "$LIMIT" ] && LIMIT=100
  NPROC=$(pgrep -f "file=scripts/" | wc -l | tr -d ' ')

  state=OK
  awk -v b="$BATT" -v w="$WARN" 'BEGIN{exit !(b>=w)}' && state=WARN
  awk -v b="$BATT" -v c="$CRIT" 'BEGIN{exit !(b>=c)}' && state=CRIT
  [ "$LIMIT" -lt 100 ] 2>/dev/null && [ "$state" = "OK" ] && state=THROTTLED

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$(date +%H:%M:%S)" "$BATT" "$VIRT" "$LIMIT" "$NPROC" "$state" >> "$LOG"

  if [ "$state" = "CRIT" ]; then
    hot_streak=$((hot_streak+1))
    if [ "$hot_streak" -ge 2 ]; then
      echo "CRITICAL ${BATT}C for 2 reads - stopping the R jobs now"
      pkill -f "file=scripts/"
      printf "%s\tKILLED by watchdog at %sC\n" "$(date +%H:%M:%S)" "$BATT" >> "$LOG"
      exit 1
    fi
  else
    hot_streak=0
  fi

  # only speak when something changes
  if [ "$state" != "$prev_state" ]; then
    echo "$(date +%H:%M) thermal: ${BATT}C batt / ${VIRT}C virt / cpu_limit ${LIMIT}% / ${NPROC} R procs -> $state"
    prev_state=$state
  fi

  # stop watching once the work is done
  if [ "$NPROC" -eq 0 ] && grep -q "ALL DONE" results/tables/queue_log.txt 2>/dev/null; then
    echo "$(date +%H:%M) queue finished, final temp ${BATT}C - watchdog standing down"
    exit 0
  fi
  sleep "$INTERVAL"
done
