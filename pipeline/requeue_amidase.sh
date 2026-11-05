#!/bin/bash
# Requeue amidase-signature (2026-10-27): its L90 step 3 on upstream settings was killed by the disk guard twice
# (00:17, 12:26; >690 GB unlinked spill on 818k seqs), so L90 step 3 runs in linear mode (deviations.tsv).
# Same chain settings as launch.sh's big pool. CHAIN done/FAILED lines go to chain.log for notify_watch.sh.
cd /watson/rnalab/dzhu/petadex-clustering/run2
export HEAVY_MEM=580G HEAVY_THREADS=200 SMALL_MAX=250000 FORCE_LINEAR_NR="L90:amidase-signature"
echo "[$(date '+%F %T')] REQUEUE amidase-signature: L90 step 3 forced to faster_lin fast_lin" >> logs/chain.log
if DIAMOND_THREADS=64 DIAMOND_MEM=40G P6=4 T6=16 ./run_chain.sh amidase-signature > logs/chain.amidase-signature.out 2>&1
then echo "[$(date '+%F %T')] CHAIN amidase-signature done" >> logs/chain.log
else echo "[$(date '+%F %T')] CHAIN amidase-signature FAILED" >> logs/chain.log; fi
