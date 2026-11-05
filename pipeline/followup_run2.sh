#!/bin/bash
# Run 2 follow-up, v3 (2026-11-02): amidase-signature is the only unfinished domain (14 of 15 done). Runs at ABOUT HALF the
# server so the other half stays free for k2-circ and Thomas: Logan back-map chunks P6=8 x 30 threads = 240 threads, all
# at nice 10 with temp files on /watson (SHM_START_MIN=9999 keeps the wrapper off /dev/shm, which Thomas may want), and
# the heavy deepclust steps capped at HEAVY_MEM=350G / HEAVY_THREADS=200 (was 580G). v2 is in archive_20261026/.
# Writes "ALL CHAINS FINISHED (final)" -- the only finish line notify_watch.sh acts on.
cd /watson/rnalab/dzhu/petadex-clustering/run2
export HEAVY_MEM=350G HEAVY_THREADS=200 SMALL_MAX=250000 FORCE_LINEAR_NR="L90:amidase-signature"
export SHM_START_MIN=9999 SHM_GAP=120
echo 30 > logs/.chunk_threads
chain() { if ./run_chain.sh $1 > logs/chain.$1.out 2>&1; then echo "[$(date '+%F %T')] CHAIN $1 done" >> logs/chain.log
          else echo "[$(date '+%F %T')] CHAIN $1 FAILED" >> logs/chain.log; fi; }
pgrep -u $USER -f 'run2/disk_guard.sh|\./disk_guard.sh' > /dev/null || nohup setsid ./disk_guard.sh > logs/disk_guard.out 2>&1 < /dev/null &
touch logs/.followup_started   # no retry is queued: a FAILED from here is final
echo "[$(date '+%F %T')] FOLLOW-UP v3: resume amidase-signature at ~half the server (P6=8 x 30 threads, nice 10, temp on /watson, heavy steps 350G/200 threads)" >> logs/chain.log
DIAMOND_THREADS=64 DIAMOND_MEM=40G P6=8 T6=30 chain amidase-signature
echo "[$(date '+%F %T')] ALL CHAINS FINISHED (final)" >> logs/chain.log
