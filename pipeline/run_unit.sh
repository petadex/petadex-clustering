#!/bin/bash
# Run Benjamin Buchfink's unid/ pipeline (upstream commit in unid/UPSTREAM_COMMIT) for one unit (CATH domain, or
# component for the control) at one identity level.
#   usage: run_unit.sh <level dir> <unit>      e.g. run_unit.sh /…/run2/L90 cytochrome-p450
# The level dir holds cfg.sh; its parent holds `domains` (unit list, read by 6.sh) and unid/.
# Deviations from upstream, none of which change what is aligned or with which settings:
#  - Step 6 (Logan back-map) runs its chunks in parallel ($P6 processes x $T6 threads), split into ~$NCHUNK chunks by
#    residue count when the input exceeds $SPLIT_MIN residues (upstream: array job, per-domain split sizes).
#  - Steps 9-13 (self-alignment + greedy vertex cover of the Logan reps) are skipped: upstream step 14 reads
#    unmapped2.clust, not unmapped2.clust_merged (changed in upstream commit 71a6e9e), so their output is unused.
#  - Resumable: each finished step leaves .done/<unit>.<step>.
#  - Non-linear deepclust (step 3 at every level; step 8 when its cluster steps are not all *_lin) runs one at a time
#    across Run 2 under logs/.heavy.lock with $HEAVY_MEM / $HEAVY_THREADS. With 4 in parallel at a 350 GB cap
#    they spilled >1 TB of unlinked temp files to /watson in 10 min (2026-10-26). disk_guard.sh kills a heavy
#    process (it is the only one with -M $HEAVY_MEM) when /watson free space drops below 400 GB.
#    A heavy step whose input has < $SMALL_MAX sequences skips the lock and runs with the chain's own
#    DIAMOND_MEM/THREADS (spill is a few GB), so small domains use the cores the lock would leave idle.
#  - Step 3 for $FALLBACK_UNITS (abhydrolase, amidase-signature): always under the heavy lock (so the disk guard
#    covers their spill), and if killed (guard SIGTERM or OOM) retried with --cluster-steps faster_lin fast_lin
#    (Run 2 option B), recorded in deviations.tsv. $FORCE_LINEAR_NR ("<level>:<unit> ...") skips straight to linear
#    mode, without the lock (linear spill is small): set for L90 amidase-signature after the guard killed its default
#    run twice (2026-10-27 00:17, 12:26).
set -euo pipefail
L=$(readlink -f $1); d=$2
R=$(dirname $L)
export PATH=$R/bin:$R/unid:/watson/rnalab/dzhu/.conda/envs/unid-clustering/bin:$PATH
export LC_ALL=C
: ${DIAMOND_THREADS:=64} ${P6:=8} ${T6:=24} ${NCHUNK:=40} ${SPLIT_MIN:=200000000}
: ${HEAVY_MEM:=580G} ${HEAVY_THREADS:=200} ${SMALL_MAX:=250000}
: ${FALLBACK_UNITS:="abhydrolase amidase-signature"} ${FORCE_LINEAR_NR:=}
export DIAMOND_THREADS
cd $L; mkdir -p .done logs reps/pazy reps/nr reps/logan
source cfg.sh
line=$(grep -nx "$d" ../domains | cut -d: -f1)
lvl=$(basename $L)
log() { echo "[$(date '+%F %T')] $lvl $d $*" | tee -a $L/logs/progress.log; }
done_() { [[ -e .done/$d.$1 ]]; }
step() { local n=$1; shift; if done_ $n; then return; fi; local t0=$SECONDS
         "$@" > logs/$d.$n.log 2>&1 || return $?; touch .done/$d.$n; log "step $n done $((SECONDS-t0))s"; }
heavy() { local n=$1 in=$2; shift 2; if done_ $n; then return; fi
          local k=$(nseq $in)
          if (( k < SMALL_MAX )); then log "step $n input $k seqs < $SMALL_MAX: no heavy lock"; step $n "$@"; return; fi
          local DIAMOND_MEM=$HEAVY_MEM DIAMOND_THREADS=$HEAVY_THREADS; export DIAMOND_MEM DIAMOND_THREADS
          exec 9>>$R/logs/.heavy.lock; log "step $n waiting for heavy lock"; flock 9; log "step $n holding heavy lock"
          local rc=0; step $n "$@" || rc=$?
          flock -u 9; exec 9>&-; return $rc; }
disk_free() { df -BG --output=avail /watson | tail -1 | tr -dc 0-9; }
nseq() { if [[ -s $1 ]]; then grep -c '^>' "$1" || true; else echo 0; fi; }
nres() { if [[ -s $1 ]]; then grep -v '^>' "$1" | tr -d '\n' | wc -c; else echo 0; fi; }

trap 'log "FAILED at line $LINENO (exit $?)"' ERR
log "begin (line $line) data_dir=$data_dir"
step 1 ../unid/1-cluster-pazy.sh $d
step 2 ../unid/2-map-nr-vs-pazy.sh $d
# Small units can map every nr core onto PAZy; DIAMOND rejects the empty file, so do step 3's no-op by hand.
if ! done_ 3 && [[ $(nseq nr_vs_pazy_$d.unmapped.faa) -eq 0 ]]; then
  : > nr_$d.clust; : > nr_$d.reps.faa; cat pazy_$d.reps.faa > pazy_nr_$d.reps.faa
  touch .done/$d.3; log "step 3 skipped: every nr core mapped to PAZy"
fi
# step 3 in linear mode: $1 = reason (deviations.tsv), $2 = lock|nolock
lin3() { log "step 3: --cluster-steps faster_lin fast_lin ($1)"
         printf '%s\t%s\tstep 3\t--cluster-steps faster_lin fast_lin (upstream: %s)\t%s\t%s\n' \
           $lvl $d "$steps_nr" "$1" "$(date '+%F %T')" >> $R/deviations.tsv
         [[ -e logs/$d.3.log ]] && mv logs/$d.3.log logs/$d.3.killed.log
         export STEPS_NR_OVERRIDE="faster_lin fast_lin"
         if [[ $2 == lock ]]; then while (( $(disk_free) < 600 )); do sleep 60; done; heavy 3 nr_vs_pazy_$d.unmapped.faa ../unid/3-cluster-nr.sh $d
         else step 3 ../unid/3-cluster-nr.sh $d; fi
         unset STEPS_NR_OVERRIDE; }
if [[ " $FALLBACK_UNITS " == *" $d "* ]]; then
  SMALL_MAX=0   # always take the heavy lock (step 3 and step 8) so the disk guard covers the spill
  if ! done_ 3; then
    if [[ " $FORCE_LINEAR_NR " == *" $lvl:$d "* ]]; then
      lin3 "default mode killed by the disk guard twice (exit 143)" nolock
    else
      rc=0; heavy 3 nr_vs_pazy_$d.unmapped.faa ../unid/3-cluster-nr.sh $d || rc=$?
      if (( rc != 0 )); then
        if (( rc != 137 && rc != 143 )); then log "step 3 FAILED (exit $rc), not a kill: no fallback"; exit $rc; fi
        log "step 3 killed (exit $rc): falling back to --cluster-steps faster_lin fast_lin"
        lin3 "default mode killed, exit $rc" lock
      fi
    fi
  fi
fi
heavy 3 nr_vs_pazy_$d.unmapped.faa ../unid/3-cluster-nr.sh $d
step 4 ../unid/4-mapback-logan-lin.sh $d

# A chunk is finished when its log ends with DIAMOND's closing "Peak RSS" line (6.sh runs under set -e; out_<i> is
# created when the chunk starts, so its existence alone proves nothing). A re-run keeps the split (split-$d/.split_ok)
# and runs only unfinished chunks, so a kill loses only the chunks in flight (2026-10-27, before Thomas's priority run).
chunk_done() { [[ -e split-$d/out_$1 ]] && tail -n 1 logs/$d.6.$1.log 2>/dev/null | grep -q 'Peak RSS'; }
if ! done_ 6; then
  u=logan_vs_nr_$d.unmapped.faa
  if [[ $(nseq $u) -gt 0 ]]; then
    if [[ ! -e split-$d/.split_ok ]]; then
      rm -rf split-$d; mkdir -p split-$d
      res=$(nres $u)
      if (( res > SPLIT_MIN )); then size=$(( res / NCHUNK + 1 )); else size=$(( res + 1 )); fi
      # a chunk must hold the longest record
      split-fasta $u $size split-$d seqs > logs/$d.5.log 2>&1
      touch split-$d/.split_ok
    fi
    nch=$(ls split-$d/seqs_*.fasta | wc -l)
    todo=$(ls split-$d/seqs_*.fasta | sed 's/.*_//; s/\.fasta//' | while read -r i; do chunk_done $i || echo $i; done)
    ntodo=$(wc -w <<< "$todo")
    (( ntodo < nch )) && log "step 6 resuming: $(( nch - ntodo )) of $nch chunks already finished"
    t0=$SECONDS
    if (( nch > 1 )); then
      [[ -n $todo ]] && printf '%s\n' $todo \
        | DIAMOND_THREADS=$T6 xargs -P $P6 -I{} sh -c "bash ../unid/6.sh $line {} > logs/$d.6.{}.log 2>&1"
    elif [[ -n $todo ]]; then
      bash ../unid/6.sh $line 1 > logs/$d.6.1.log 2>&1
    fi
    # every chunk must have finished
    for f in split-$d/seqs_*.fasta; do i=${f##*_}; i=${i%.fasta}; chunk_done $i || { log "FAILED step 6 chunk $i"; exit 1; }; done
    log "step 6 done $((SECONDS-t0))s ($nch chunks)"
  else
    rm -rf split-$d; mkdir -p split-$d
    : > split-$d/out_0; log "step 5/6 skipped: no Logan left after step 4"
  fi
  touch .done/$d.6
fi
step 7 ../unid/7-get-unmapped2.sh $d
if [[ $(nseq logan_vs_nr_$d.unmapped2.faa) -gt 0 ]]; then
  if tr ' ' '\n' <<< "$steps_logan_linclust" | grep -qv '_lin$'; then heavy 8 logan_vs_nr_$d.unmapped2.faa ../unid/8-linclust.sh $d
  else step 8 ../unid/8-linclust.sh $d; fi
else
  : > logan_vs_nr_$d.unmapped2.clust; touch .done/$d.8; log "step 8 skipped: no Logan left after step 6"
fi
step 14 ../unid/14-merge-all.sh $d
step 15 ../unid/15-get-reps.sh $d
if ! done_ check; then
  n_in=$(( $(nseq $data_dir/pazy/${d}_pazy_catalytic_domains.fa) + $(nseq $data_dir/nr/${d}_blastnr_catalytic_domains.fa) + $(nseq $data_dir/logan/${d}_logan_catalytic_domains.fa) ))
  n_mem=$(cut -f2 $d.all.tsv | sort -u | wc -l); n_rows=$(wc -l < $d.all.tsv)
  n_cl=$(cut -f1 $d.all.tsv | sort -u | wc -l)
  n_rep=$(( $(nseq reps/pazy/${d}_pazy_catalytic_domains.fa) + $(nseq reps/nr/${d}_blastnr_catalytic_domains.fa) + $(nseq reps/logan/${d}_logan_catalytic_domains.fa) ))
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' $d $n_in $n_rows $n_mem $n_cl $n_rep >> counts.tsv
  if [[ $n_in -ne $n_rows || $n_in -ne $n_mem || $n_cl -ne $n_rep ]]; then log "CHECK FAILED: input $n_in rows $n_rows members $n_mem clusters $n_cl reps $n_rep"; exit 1; fi
  touch .done/$d.check
  log "done: input $n_in, clusters $n_cl"
fi
