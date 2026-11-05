#!/bin/bash
# Run 2 launcher, fifth version (resume after stop_run2.sh, 2026-10-27; resume_run2.sh installs this as launch.sh).
# Earlier versions are in archive_20261026/. Requeues every domain without a "CHAIN <d> done" line in chain.log, in one
# pool. Finished steps (.done) and finished Logan chunks (run_unit.sh chunk_done) are skipped. Heavy deepclust is still
# serialised by logs/.heavy.lock; P6=5 because 23 Logan chunks at once filled /dev/shm (378 GB) on 10-27 19:52.
cd /watson/rnalab/dzhu/petadex-clustering/run2
export HEAVY_MEM=580G HEAVY_THREADS=200 SMALL_MAX=250000 FORCE_LINEAR_NR="L90:amidase-signature"
fail() { echo "[$(date '+%F %T')] CHAIN $1 FAILED" >> logs/chain.log; }
export -f fail
pgrep -u $USER -f 'run2/disk_guard.sh|\./disk_guard.sh' > /dev/null || nohup setsid ./disk_guard.sh > logs/disk_guard.out 2>&1 < /dev/null &
pending=$(for d in $(cat domains); do grep -q "CHAIN $d done" logs/chain.log || echo $d; done)
echo "[$(date '+%F %T')] RELAUNCH 5 (resume): $(echo $pending)" >> logs/chain.log
printf '%s\n' $pending | DIAMOND_THREADS=64 DIAMOND_MEM=40G P6=5 T6=20 \
  xargs -P $(wc -w <<< "$pending") -I{} bash -c './run_chain.sh {} > logs/chain.{}.out 2>&1 && echo "[$(date "+%F %T")] CHAIN {} done" >> logs/chain.log || fail {}'
echo "[$(date '+%F %T')] ALL CHAINS FINISHED" >> logs/chain.log
