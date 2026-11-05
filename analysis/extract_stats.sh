#!/bin/bash
# Per level x domain cluster statistics from run2/L*/<domain>.all.tsv (rep <tab> member, grouped by rep).
# Output data/cluster_stats.tsv: level domain clusters singletons members PL NR LG other   (centroid source by ID prefix)
R=/watson/rnalab/dzhu/petadex-clustering/run2; D=$(dirname $(dirname $(readlink -f $0)))/data; mkdir -p $D/parts
one() { L=${1%%:*}; d=${1#*:}; f=$R/$L/$d.all.tsv
  awk -F'\t' -v L=$L -v d=$d 'function flush(){ if(n==0) return; c++; if(n==1) s++; src=substr(p,1,2); if(src=="PL")pl++; else if(src=="NR")nr++; else if(src=="LG")lg++; else ot++ }
    $1!=p { flush(); p=$1; n=0 } { n++; m++ }
    END{ flush(); printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", L, d, c, s, m, pl, nr, lg, ot }' $f > $D/parts/$L.$d.tsv; }
export -f one; export R D
for L in L90 L60 L30; do for d in $(cut -f1 $R/$L/counts.tsv); do [[ -s $R/$L/$d.all.tsv ]] && echo $L:$d; done; done | xargs -P 6 -I{} nice -n 19 bash -c 'one {}'
{ printf "level\tdomain\tclusters\tsingletons\tmembers\tPL\tNR\tLG\tother\n"; cat $D/parts/*.tsv | sort -k1,1r -k2,2; } > $D/cluster_stats.tsv
# per-level input/cluster counts written by run_unit.sh (domain input rows allrows members clusters reps), kept with the figure data
for L in L90 L60 L30; do cp $R/$L/counts.tsv $D/counts_$L.tsv; done
