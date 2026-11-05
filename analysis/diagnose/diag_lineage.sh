#!/bin/bash
# Pass B: how many original cores sit under each cluster, following centroids down the hierarchy.
#   weight(90% cluster) = its size; weight(60% cluster) = sum of weights of its 90% members; weight(30%) likewise.
# Output data/diagnostics/lineage_by_class.tsv: level domain source singleton clusters cores  w1 w2 w3_5 w6_10 w11_100 w101_1000 wgt1000
R=/watson/rnalab/dzhu/petadex-clustering/run2; ROOT=$(dirname $(dirname $(dirname $(readlink -f $0)))); D=$ROOT/data/diagnostics/parts
W=${WORK:-/tmp/lineage_work}; mkdir -p $W
one() { d=$1; export LC_ALL=C
  awk -F'\t' '$1!=p{ if(n) print p"\t"n; p=$1; n=0 } {n++} END{ if(n) print p"\t"n }' $R/L90/$d.all.tsv > $W/$d.w90
  awk -F'\t' 'NR==FNR{w[$1]=$2; next} $1!=p{ if(n) print p"\t"n"\t"s; p=$1; n=0; s=0 } {n++; s+=w[$2]} END{ if(n) print p"\t"n"\t"s }' $W/$d.w90 $R/L60/$d.all.tsv > $W/$d.w60
  awk -F'\t' 'NR==FNR{w[$1]=$3; next} $1!=p{ if(n) print p"\t"n"\t"s; p=$1; n=0; s=0 } {n++; s+=w[$2]} END{ if(n) print p"\t"n"\t"s }' $W/$d.w60 $R/L30/$d.all.tsv > $W/$d.w30
  for L in 90 60 30; do
    if [[ $L == 90 ]]; then cols='$1,$2,$2'; else cols='$1,$2,$3'; fi
    awk -F'\t' -v L=L$L -v d=$d 'function wb(x){return x==1?0:x==2?1:x<=5?2:x<=10?3:x<=100?4:x<=1000?5:6}
      { s=substr($1,1,2); g=($2==1); k=s SUBSEP g; c[k]++; cores[k]+=$3; b[k,wb($3)]++ }
      END{ for(k in c){ split(k,q,SUBSEP); printf "%s\t%s\t%s\t%d\t%d\t%d", L,d,q[1],q[2],c[k],cores[k]; for(i=0;i<7;i++) printf "\t%d", b[k,i]; printf "\n" } }' \
      <(awk -F'\t' "{print $cols}" OFS='\t' $W/$d.w$L) > $D/B.L$L.$d.tsv
  done
  rm -f $W/$d.w90 $W/$d.w60 $W/$d.w30; }
export -f one; export R W D
cut -f1 $R/L30/counts.tsv | xargs -P 5 -I{} nice -n 19 bash -c 'one {}'
{ printf "level\tdomain\tsource\tsingleton\tclusters\tcores\tw1\tw2\tw3_5\tw6_10\tw11_100\tw101_1000\tw_gt1000\n"; cat $D/B.*.tsv; } > $D/../lineage_by_class.tsv
echo DONE > $D/../diag_B.done
