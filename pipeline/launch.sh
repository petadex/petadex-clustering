#!/bin/bash
# Run 2 launcher, fourth version (2026-10-27 relaunch 4, after the 10-26 23:58-00:17 OOM storm + disk-guard kill).
# Earlier versions are in archive_20261026/. Still running from earlier launchers and NOT restarted here:
#   abhydrolase chain (pid 3508285) and the big pool xargs (pid 3555352: cytochrome-p450, peptidase, trypsin).
# This script requeues what died (finished steps are kept via .done markers) and writes ALL CHAINS FINISHED only
# after its own pools AND those pids have exited.
#   big   - betalactamase (OOM-killed in L90 step 3), amidase-signature (disk guard, L90 step 3); queue on the lock
#   small - the 8 unfinished small domains; DIAMOND_OOM_ADJ=1000 on their DIAMOND processes only, 16G caps, P6=2;
#           SMALL_MAX=700000 so their L60 step 8 (258k-654k seqs) skips the lock instead of starving the pool
cd /watson/rnalab/dzhu/petadex-clustering/run2
export HEAVY_MEM=580G HEAVY_THREADS=200 SMALL_MAX=250000
WAIT_PIDS="3508285 3555352"
fail() { echo "[$(date '+%F %T')] CHAIN $1 FAILED" >> logs/chain.log; }
export -f fail
pgrep -u $USER -f 'run2/disk_guard.sh|\./disk_guard.sh' > /dev/null || nohup setsid ./disk_guard.sh > logs/disk_guard.out 2>&1 < /dev/null &
echo "[$(date '+%F %T')] RELAUNCH 4b: requeue betalactamase, amidase-signature + 8 small domains after OOM storm" >> logs/chain.log
pool() { xargs -P $1 -I{} bash -c './run_chain.sh {} > logs/chain.{}.out 2>&1 && echo "[$(date "+%F %T")] CHAIN {} done" >> logs/chain.log || fail {}'; }
export -f pool
printf '%s\n' betalactamase amidase-signature \
  | DIAMOND_THREADS=64 DIAMOND_MEM=40G P6=4 T6=16 pool 2 &
printf '%s\n' acid-protease d-ala-esterase oxygenase-like cytochromec-like peroxidase sgnh-hydrolase cupredoxins \
  metallo-phosphatase \
  | SMALL_MAX=700000 DIAMOND_OOM_ADJ=1000 DIAMOND_THREADS=64 DIAMOND_MEM=16G P6=2 T6=24 pool 3 &
wait
for p in $WAIT_PIDS; do while kill -0 $p 2>/dev/null; do sleep 60; done; done
echo "[$(date '+%F %T')] ALL CHAINS FINISHED" >> logs/chain.log
