#!/bin/bash
set -e
domain=$1
code=$2
file=$3
id=${code}_${domain}
source cfg.sh
# Run 2 local: run_unit.sh's abhydrolase fallback (option B) overrides the step list; unset = upstream behaviour
steps_nr=${STEPS_NR_OVERRIDE:-$steps_nr}
diamond deepclust -d $file -o $id.clust --cluster-steps $steps_nr -M 2000G \
        --ext full --comp-based-stats 0 --masking 0 --log --id $id_cutoff --member-cover 90

seqtk subseq $file <(cut -f1 $id.clust | uniq) > $id.reps.faa
