#!/bin/bash
# 2026-10-27: when peptidase-s8-s53 L90 step 6 finishes, give its 4 extra chunk slots to amidase-signature L90 step 6
# (xargs pid 1390616, P 4 -> up to 8 via SIGUSR1), one slot per 3 min and only while MemAvailable > 150 GB.
cd /watson/rnalab/dzhu/petadex-clustering/run2
X=1390616
memavail() { awk '/^MemAvailable/{print int($2/1048576)}' /proc/meminfo; }
until [[ -e L90/.done/peptidase-s8-s53.6 ]]; do kill -0 $X 2>/dev/null || exit 0; sleep 120; done
n=0
while (( n < 4 )) && kill -0 $X 2>/dev/null; do
  if (( $(memavail) > 150 )); then kill -USR1 $X; n=$((n+1)); fi
  sleep 180
done
(( n > 0 )) && echo "[$(date '+%F %T')] amidase-signature L90 step 6 parallelism 4->$((4+n)) (xargs SIGUSR1, peptidase step 6 finished)" >> logs/chain.log
