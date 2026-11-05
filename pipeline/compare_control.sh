#!/bin/bash
# Grouping effect at 30 %: per-component (control30) vs per-CATH-domain (unid30_run), same input + recipe.
# Singleton = 30 % cluster holding exactly one input domain (n_domains summed from Run 1 cluster_stats_60).
set -euo pipefail
P=/watson/rnalab/dzhu/petadex-clustering
OUT=$P/run2/control30/grouping_effect.tsv
awk -F'\t' '
  FNR==1 { file++ }
  file==1 { if (FNR>1) nd[$1]=$5; next }                      # cluster_stats_60: centroid -> n_domains
  file==2 { dom[$2]=$3; d_n[$2]+=nd[$1]; next }   # per-domain run: c60 c30 domain
  file==3 { c_n[$1]+=nd[$2]; next }        # control: c30 c60
  END {
    for (k in d_n) { split(k,a,"|"); D=dom[k]; dc[D]++; if (d_n[k]==1) ds[D]++ }
    printf "unit\tclusters\tsingletons\n"
    for (D in dc) printf "domain:%s\t%d\t%d\n", D, dc[D], ds[D]
    for (k in c_n) { split(k,a,"|"); C=a[6]; cc[C]++; if (c_n[k]==1) cs[C]++ }
    for (C in cc) printf "component:%s\t%d\t%d\n", C, cc[C], cs[C]
  }' $P/results/cluster_stats_60.tsv $P/unid30_run/results/map_60_to_30.unid.tsv <(cat $P/run2/control30/L30/C*.all.tsv) > $OUT.tmp
# roll components up to their domain
awk -F'\t' 'NR==FNR{m[$1]=$2; next} FNR==1{next} $1~/^domain:/{sub("domain:","",$1); dc[$1]=$2; ds[$1]=$3; next}
  {sub("component:","",$1); D=m[$1]; cc[D]+=$2; cs[D]+=$3}
  END{printf "domain\tper_domain_clusters\tper_component_clusters\tper_domain_singletons\tper_component_singletons\n";
      for(D in dc){printf "%s\t%d\t%d\t%d\t%d\n",D,dc[D],cc[D],ds[D],cs[D]; a+=dc[D];b+=cc[D];c+=ds[D];d+=cs[D]}
      printf "TOTAL\t%d\t%d\t%d\t%d\n",a,b,c,d}' $P/run2/component_domain.tsv $OUT.tmp | (IFS= read -r h; echo "$h"; sort -t$'\t' -k2,2nr) > $OUT
rm $OUT.tmp; column -t $OUT
