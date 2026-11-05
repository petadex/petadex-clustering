#!/bin/bash
# Rows assigned by each mapping step, per level, summed over the 15 domains (wc -l of run2/L*/<files>).
# Output data/method_flow_counts.tsv: level  nr_to_pazy  logan_step4  logan_step6
R=/watson/rnalab/dzhu/petadex-clustering/run2; D=$(dirname $(dirname $(readlink -f $0)))/data; mkdir -p $D/flow_parts
one() { L=$1; d=$2
  printf "%s\t%s\t%d\t%d\t%d\n" $L $d $(wc -l < $R/$L/nr_vs_pazy_$d.aln) $(wc -l < $R/$L/logan_vs_nr_$d.aln) $(wc -l < $R/$L/logan_vs_nr_$d.aln2) > $D/flow_parts/$L.$d.tsv; }
export -f one; export R D
for L in L90 L60 L30; do for d in $(cut -f1 $R/$L/counts.tsv); do echo "$L $d"; done; done | xargs -P 6 -L1 nice -n 19 bash -c 'one $0 $1'
{ printf "level\tdomain\tnr_to_pazy\tlogan_step4\tlogan_step6\n"; cat $D/flow_parts/*.tsv | sort -k1,1r -k2,2; } > $D/method_flow_by_domain.tsv
awk -F'\t' 'NR>1{a[$1]+=$3;b[$1]+=$4;c[$1]+=$5} END{print "level\tnr_to_pazy\tlogan_step4\tlogan_step6"; for(l in a) print l"\t"a[l]"\t"b[l]"\t"c[l]}' $D/method_flow_by_domain.tsv | sort -r > $D/method_flow_counts.tsv
rm -rf $D/flow_parts; echo DONE > $D/method_flow.done
