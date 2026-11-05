#!/bin/bash
# Appends a live "Morning snapshot" to the Run 2 overnight report and pushes it to the vault (2026-10-29).
# Sleeps until $AT (default 11:30 UTC), so it is launched once under watchjob. REPORT=<file> NOPUSH=1 AT=now for a test.
R=/watson/rnalab/dzhu/petadex-clustering/run2; V=/watson/rnalab/dzhu/obsidian-vault
F=${REPORT:-"$V/PETadex/Projects/Catalytic-Core Clustering v260701/03 Run 2 Overnight Report 2026-10-29.md"}
t=$(date -d "${AT:-2026-10-29 11:30 UTC}" +%s); n=$(date +%s); (( t > n )) && sleep $((t - n))
cd $R
since="2026-10-29 04:49"; t0=$(date -d "$since UTC" +%s)
chunks_done() { local d=$1 c=0 i; for i in $(ls L90/split-$d/seqs_*.fasta | sed 's/.*_//; s/\.fasta//'); do
                  tail -n 1 L90/logs/$d.6.$i.log 2>/dev/null | grep -q 'Peak RSS' && c=$((c+1)); done; echo $c; }
ndone=$(grep -o 'CHAIN [a-z0-9-]* done' logs/chain.log | sort -u | wc -l)
{
echo
echo "## Morning snapshot ($(date -u '+%F %H:%M') UTC, written automatically)"
echo
echo "- **Domains done:** $ndone of 15."
for d in betalactamase amidase-signature; do
  if grep -q "CHAIN $d done" logs/chain.log; then echo "- **$d:** done ($(grep "CHAIN $d done" logs/chain.log | tail -1 | cut -c2-17) UTC)."; continue; fi
  last=""; for L in L30 L60 L90; do last=$(grep " $d " $L/logs/progress.log 2>/dev/null | tail -1); [[ -n $last ]] && break; done
  echo "- **$d:** L90 Logan chunks $(chunks_done $d)/$(ls L90/split-$d/seqs_*.fasta | wc -l); last step line: \`${last}\`"
done
a=$(chunks_done amidase-signature); el=$(( ($(date +%s) - t0) / 3600 ))
if (( a > 7 && el > 0 )); then r=$(awk -v d=$((a-7)) -v h=$el 'BEGIN{printf "%.2f", d/h}'); left=$((41 - a))
  echo "- **Amidase rate:** $((a-7)) chunks finished in ${el} h since the 04:49 restart (${r}/h); ${left} left → about $(awk -v l=$left -v r=$r 'BEGIN{printf "%.0f", l/r}') h more for its L90 Logan step, then its L60/L30 levels."
else echo "- **Amidase rate:** no new chunk finished yet since the 04:49 restart (${el} h; each takes ~16–24 h, so the first ones land around midday or later)."; fi
echo "- **OOM kills since the 10-27 resume:** $(( $(awk '/^oom_kill /{print $2}' /proc/vmstat) - 390 ))."
echo "- **Machine:** load $(cut -d' ' -f1 /proc/loadavg) / 480; $(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo) GB memory available; /dev/shm $(df -BG --output=avail /dev/shm | tail -1 | tr -dc 0-9) GB free; /watson $(df -BG --output=avail /watson | tail -1 | tr -dc 0-9) GB free."
echo "- **Thomas (tquigley):** $(ps -u tquigley -o pcpu= 2>/dev/null | awk '{s+=$1} END{printf "%d", s/100}') cores in use."
echo
echo "**Run 2 log since the 04:49 restart:**"
echo
echo '```'
awk -v s="[$since" '$0 >= s' logs/chain.log | tail -n 25
echo '```'
echo
echo "**Pushes sent since then:**"
echo
awk -F'\t' -v s="$since" '$1 >= s && ($4 ~ /Run 2|watson/) {print "- " $1 " — " $4}' ~/.config/run-notify/sent.log
} >> "$F"
[[ -n $NOPUSH ]] && exit 0
cd $V && git pull -q --rebase --autostash && git add "$F" && git commit -qm "Run 2 overnight report: morning snapshot" && git push -q
