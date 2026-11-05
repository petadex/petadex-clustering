#!/bin/bash
# Exact cluster-size histogram per level (size = direct members: cores at 90%, 90% centroids at 60%, 60% centroids at 30%).
# Output ../../data/diagnostics/cluster_size_hist.tsv: level size clusters
R=/watson/rnalab/dzhu/petadex-clustering/run2; D=$(dirname $(dirname $(dirname $(readlink -f $0))))/data/diagnostics; mkdir -p $D/size_parts
one() { L=${1%%:*}; d=${1#*:}; LC_ALL=C awk -F'\t' -v L=$L '$1!=p{ if(n) h[n]++; p=$1; n=0 } {n++} END{ if(n) h[n]++; for(s in h) print L"\t"s"\t"h[s] }' $R/$L/$d.all.tsv > $D/size_parts/$L.$d.tsv; }
export -f one; export R D
for L in L90 L60 L30; do for d in $(cut -f1 $R/$L/counts.tsv); do echo $L:$d; done; done | xargs -P 6 -I{} nice -n 19 bash -c 'one {}'
{ printf "level\tsize\tclusters\n"; cat $D/size_parts/*.tsv | awk -F'\t' '{h[$1"\t"$2]+=$3} END{for(k in h) print k"\t"h[k]}' | sort -k1,1r -k2,2n; } > $D/cluster_size_hist.tsv
rm -rf $D/size_parts; echo DONE > $D/size_hist.done
