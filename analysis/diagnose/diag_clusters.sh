#!/bin/bash
# Pass A: per level x domain cluster properties from run2/L*/<domain>.all.tsv (rep<TAB>member, grouped by rep).
# Writes data/diagnostics/parts/A.<level>.<domain>.tsv (class table) and S.<level>.<domain>.tsv (size table).
R=/watson/rnalab/dzhu/petadex-clustering/run2; D=$(dirname $(dirname $(dirname $(readlink -f $0))))/data/diagnostics/parts
one() { L=${1%%:*}; d=${1#*:}
  LC_ALL=C awk -F'\t' -v L=$L -v d=$d -v out=$D/A.$L.$d.tsv -v out2=$D/S.$L.$d.tsv '
   function lb(x){return x<100?0:x<150?1:x<250?2:x<400?3:4}
   function eb(x){return x<=-50?0:x<=-30?1:x<=-20?2:x<=-10?3:4}
   function sb(n){return n==1?0:n==2?1:n<=5?2:n<=10?3:n<=100?4:n<=1000?5:n<=10000?6:7}
   function flush(   a,s,g,len,ev,b){ if(n==0) return; split(p,a,"|"); s=substr(p,1,2); g=(n==1)
     len=a[5]-a[4]+1; ev=(a[3]+0>0)?log(a[3]+0)/log(10):-300
     k=s SUBSEP g; cnt[k]++; sl[k]+=len; lbn[k,lb(len)]++; ebn[k,eb(ev)]++
     z=s SUBSEP sb(n); sc[z]++; sm[z]+=n }
   $1!=p { flush(); p=$1; n=0 } { n++ }
   END{ flush(); for(k in cnt){ split(k,q,SUBSEP); printf "%s\t%s\t%s\t%d\t%d\t%d", L,d,q[1],q[2],cnt[k],sl[k] > out
          for(i=0;i<5;i++) printf "\t%d", lbn[k,i] > out; for(i=0;i<5;i++) printf "\t%d", ebn[k,i] > out; printf "\n" > out }
          for(z in sc){ split(z,q,SUBSEP); printf "%s\t%s\t%s\t%d\t%d\t%d\n", L,d,q[1],q[2],sc[z],sm[z] > out2 } }' $R/$L/$d.all.tsv; }
export -f one; export R D
for L in L90 L60 L30; do for d in $(cut -f1 $R/$L/counts.tsv); do echo $L:$d; done; done | sort -t: -k1,1r | xargs -P 6 -I{} nice -n 19 bash -c 'one {}'
{ printf "level\tdomain\tsource\tsingleton\tclusters\tsum_len\tlen_lt100\tlen_100_149\tlen_150_249\tlen_250_399\tlen_ge400\tev_le50\tev_50_30\tev_30_20\tev_20_10\tev_gt10\n"; cat $D/A.*.tsv; } > $D/../clusters_by_class.tsv
{ printf "level\tdomain\tsource\tsizebin\tclusters\tmembers\n"; cat $D/S.*.tsv; } > $D/../clusters_by_size.tsv
echo DONE > $D/../diag_A.done
