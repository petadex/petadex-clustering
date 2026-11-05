#!/bin/bash
source cfg.sh
domain=$1
f="$data_dir/pazy/${domain}_pazy_catalytic_domains.fa"
echo $f
echo $id_cutoff
/watson/rnalab/dzhu/petadex-clustering/run2/unid/cluster-all-vs-all.sh $domain pazy $f
