#!/bin/bash
# Concatenate Thomas's per-component catalytic-core FASTAs into per-CATH-domain files in unid/ layout, then count.
set -euo pipefail
IN=/watson/rnalab/tquigley/petadex.v260701/dataset_generation/7_petadex_clustering/output
R=/watson/rnalab/dzhu/petadex-clustering/run2
O=$R/input90; mkdir -p $O/pazy $O/nr $O/logan
build() {  # domain
  local d=$1 comps; comps=$(awk -v d=$d '$2==d{print $1}' $R/component_domain.tsv)
  for c in $comps; do cat $IN/plasticome_catalytic_domain_per_component/petadex-plasticome.v260701.catalytic_domain.$c.fastaa; done > $O/pazy/${d}_pazy_catalytic_domains.fa
  for c in $comps; do cat $IN/nr_catalytic_domain_per_component/petadex-nr.v260701.catalytic_domain.$c.fastaa; done > $O/nr/${d}_blastnr_catalytic_domains.fa
  for c in $comps; do cat $IN/logan_catalytic_domain_per_component/petadex-logan.v260701.catalytic_domain.$c.fastaa; done > $O/logan/${d}_logan_catalytic_domains.fa
  printf '%s\t%s\t%s\t%s\n' $d $(grep -c '>' $O/pazy/${d}_pazy_catalytic_domains.fa || true) $(grep -c '>' $O/nr/${d}_blastnr_catalytic_domains.fa || true) $(grep -c '>' $O/logan/${d}_logan_catalytic_domains.fa || true)
}
export -f build; export IN R O
# smallest first so cytochrome-p450 is ready early
cat $R/domains | xargs -P 8 -I{} bash -c 'build {}' | sort > $O/counts.tsv.tmp
awk '{p+=$2;n+=$3;l+=$4} END{print "TOTAL",p,n,l,p+n+l}' $O/counts.tsv.tmp > $O/total.txt
mv $O/counts.tsv.tmp $O/counts.tsv
cat $O/total.txt
