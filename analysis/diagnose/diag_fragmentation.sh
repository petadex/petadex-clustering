#!/bin/bash
# Is Logan more fragmented than nr, and does that make its sequences harder to cluster? 90% level (inputs are raw cores).
# Per domain, from run2/L90/<domain>.all.tsv (rep<TAB>member, grouped by rep; IDs carry hmm_start|hmm_end):
#   lengths.<d>.tsv : domain source core_length n_cores                       (every input core)
#   clusters.<d>.tsv: domain centroid_source centroid_length singleton own_source_singleton n_clusters
#     own_source_singleton = the centroid has no other member from its own source (fair for nr, whose clusters also absorb Logan)
R=/watson/rnalab/dzhu/petadex-clustering/run2; D=$(dirname $(dirname $(dirname $(readlink -f $0))))/data/fragmentation
one() { d=$1; LC_ALL=C awk -F'\t' -v d=$d -v o1=$D/parts/lengths.$d.tsv -v o2=$D/parts/clusters.$d.tsv '
  function len(id,  a){ split(id,a,"|"); return a[5]-a[4]+1 }
  function flush(){ if(n==0) return; k=rs SUBSEP rl SUBSEP (n==1) SUBSEP (no==1); c[k]++ }
  $1!=p { flush(); p=$1; rs=substr($1,1,2); rl=len($1); n=0; no=0 }
  { n++; ms=substr($2,1,2); if(ms==rs) no++; h[ms SUBSEP len($2)]++ }
  END{ flush(); for(k in h){ split(k,q,SUBSEP); print d"\t"q[1]"\t"q[2]"\t"h[k] > o1 }
       for(k in c){ split(k,q,SUBSEP); print d"\t"q[1]"\t"q[2]"\t"q[3]"\t"q[4]"\t"c[k] > o2 } }' $R/L90/$d.all.tsv; }
export -f one; export R D
cut -f1 $R/L90/counts.tsv | xargs -P 6 -I{} nice -n 19 bash -c 'one {}'
{ printf "domain\tsource\tcore_length\tn_cores\n"; cat $D/parts/lengths.*.tsv; } > $D/core_lengths_L90.tsv
{ printf "domain\tsource\tcentroid_length\tsingleton\town_source_singleton\tclusters\n"; cat $D/parts/clusters.*.tsv; } > $D/clusters_by_length_L90.tsv
rm -rf $D/parts; echo DONE > $D/frag.done
