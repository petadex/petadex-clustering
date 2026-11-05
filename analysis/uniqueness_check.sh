#!/bin/bash
# Run 2 assignment/centroid uniqueness check (2026-10-29). Read-only on run2/L*/<d>.all.tsv (rep \t member).
# IDs: ORF|component|evalue|start|end|component. "core" = ORF|start|end (one catalytic core); "protein" = ORF.
set -euo pipefail
export LC_ALL=C
R=/watson/rnalab/dzhu/petadex-clustering/run2; W=$(cd "$(dirname "$0")/.." && pwd)/data/uniqueness; T=$W/tmp; mkdir -p $T
S="sort -S 3G -T $T --parallel 8"
for L in L90 L60 L30; do
  ds=$(cut -f1 $R/$L/counts.tsv)
  # members: protein dom core ; reps: protein dom core (distinct) ; selfreps
  for d in $ds; do
    awk -F'\t' -v d=$d '{split($2,m,"|"); print m[1]"\t"d"\t"m[1]"|"m[4]"|"m[5] > "/dev/stdout"
                          if($1==$2) print m[1]"\t"d"\t"m[1]"|"m[4]"|"m[5] > "/dev/stderr"}' $R/$L/$d.all.tsv 2>>$W/$L.self.raw
    cut -f1 $R/$L/$d.all.tsv | awk -F'|' -v d=$d '{print $1"\t"d"\t"$1"|"$4"|"$5}' | $S -u >> $W/$L.reps.raw
  done | $S > $W/$L.members.sorted
  $S -u $W/$L.reps.raw > $W/$L.reps.sorted; $S $W/$L.self.raw > $W/$L.self.sorted; rm $W/$L.reps.raw $W/$L.self.raw
  {
  echo "### $L ($(wc -w <<< "$ds") domains)"
  awk -F'\t' '{n++; if($3!=pc){uc++; pc=$3}else dupc++
               if($1!=pp){ if(pp!="") {np++; if(cnt>1){mp++; if(nd>1)xd++; else sd++}}; pp=$1; cnt=0; nd=0; pd=""}
               cnt++; if($2!=pd){nd++; pd=$2}}
       END{np++; if(cnt>1){mp++; if(nd>1)xd++; else sd++}
           print "member rows", n; print "distinct cores", uc; print "core rows beyond the first (same core in >1 row)", dupc+0
           print "distinct proteins", np; print "proteins in >1 row", mp+0; print "  ... within one domain", sd+0; print "  ... across domains", xd+0}' $W/$L.members.sorted
  echo "distinct centroids (rep column) $(wc -l < $W/$L.reps.sorted)"
  echo "self rows (rep==member) $(wc -l < $W/$L.self.sorted)"
  echo "centroids without a self row $(cut -f2,3 $W/$L.reps.sorted | sort -u | comm -23 - <(cut -f2,3 $W/$L.self.sorted | sort -u) | wc -l)"
  awk -F'\t' '{k=$1; if(k==p){c++} else {if(c>1){m++; if(nd>1)x++; else s++}; p=k; c=1; nd=1; pd=$2; next}; if($2!=pd){nd++; pd=$2}}
       END{if(c>1){m++; if(nd>1)x++; else s++}; print "proteins that are the centroid of >1 cluster", m+0; print "  ... within one domain", s+0; print "  ... across domains", x+0}' $W/$L.reps.sorted
  } > $W/$L.report
done
# nesting: each finished level's member set == the level above's centroid set, per domain
for pair in "L90 L60" "L60 L30"; do set -- $pair
  for d in $(cut -f1 $R/$2/counts.tsv); do
    a=$(grep -P "\t$d\t" $W/$1.self.sorted | cut -f3 | sort -u | md5sum); b=$(grep -P "\t$d\t" $W/$2.members.sorted | cut -f3 | sort -u | md5sum)
    [[ $a == "$b" ]] && echo "nest $1->$2 $d OK" || echo "nest $1->$2 $d MISMATCH"
  done
done > $W/nesting.report
echo DONE > $W/DONE
