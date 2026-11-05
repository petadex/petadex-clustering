#!/bin/bash
set -e  # Run 2 local: upstream lacked it, so a failed/killed DIAMOND was masked by the trailing command
source cfg.sh
f="nr_vs_pazy_$1.unmapped.faa"
/watson/rnalab/dzhu/petadex-clustering/run2/unid/cluster-full.sh $1 nr $f
cat pazy_$1.reps.faa nr_$1.reps.faa > pazy_nr_$1.reps.faa
