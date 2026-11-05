#!/bin/bash
# Standalone Run 2 watcher (start with nohup setsid). Sends ntfy pushes via ~/bin/notify.sh for:
#   a domain chain done / FAILED, a step FAILED (tagged OOM if the kernel OOM-killed something since the last check),
#   a linear-mode fallback (abhydrolase, amidase-signature), any kernel OOM kill on watson while Run 2 runs, /watson < 800 GB free,
#   launcher exit, all done.
# Drain mode (2026-10-27, $R/DRAIN exists): sends "Run 2 drained" once when no Run 2 DIAMOND is computing (the wrapper
# holds new ones), and does not push FAILED lines then: chains whose Logan step lost chunks fail on purpose at the end
# of their in-flight chunks and are resumed by resume_run2.sh.
# Follow-up (2026-10-28): only "ALL CHAINS FINISHED (final)" (written by followup_run2.sh) counts as the finish; a
# FAILED before followup_run2.sh starts its rerun says the domain will be retried; followup_run2.sh counts as a launcher.
# Only "CHAIN <d> done|FAILED" and "ALL CHAINS FINISHED" lines in chain.log count; other lines (notes, guard) are
# ignored. Messages carry only domain, level, step, state and exit code, never raw log lines.
# (2026-10-26 rewrite: the old version sent every chain.log line, so manual notes went out as "domain done".)
R=${RUN2_DIR:-/watson/rnalab/dzhu/petadex-clustering/run2}
N=~/bin/notify.sh
ndom=$(wc -l < $R/domains)
touch $R/logs/chain.log
c0=$(wc -l < $R/logs/chain.log)
declare -A seen
low_sent=0
drained_sent=0
shm_sent=0
stall_since=0
stall_sent=0
draining() { [[ -e $R/DRAIN ]]; }
retry_queued() { pgrep -u $USER -f 'followup_run2\.sh' > /dev/null && [[ ! -e $R/logs/.followup_started ]]; }
# a Run 2 DIAMOND that is computing: the real binary (not the wrapper waiting in drain mode) with its cwd in Run 2.
# Amidase's L90 Logan chunks are ignored (2026-10-27: ~25% done after 8 h, ~1 day left each, lost tonight regardless).
busy() { local p; for p in $(pgrep -u $USER -x diamond); do
           tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -q 'split-amidase-signature/' && continue
           [[ $(readlink /proc/$p/cwd 2>/dev/null) == $R/* && $(readlink /proc/$p/exe 2>/dev/null) != *bash ]] && return 0
         done; return 1; }
oom() { awk '/^oom_kill /{print $2}' ${VMSTAT:-/proc/vmstat}; }
oom0=$(oom)
# progress.log FAILED lines already present when the watcher starts are not re-sent
while read -r l; do seen[$l]=1; done < <(grep -h "FAILED\|falling back" $R/L90/logs/progress.log $R/L60/logs/progress.log $R/L30/logs/progress.log 2>/dev/null)
while true; do
  o=$(oom); oom_new=0
  if (( o > oom0 )); then
    oom_new=1
    $N "watson OOM kill during Run 2" "The kernel killed $((o - oom0)) process(es) for lack of memory. A Run 2 step may fail next." 5
    oom0=$o
  fi
  c=$(wc -l < $R/logs/chain.log)
  if (( c > c0 )); then
    tail -n $((c - c0)) $R/logs/chain.log | while read -r l; do
      if [[ $l =~ CHAIN\ ([a-z0-9-]+)\ done ]]; then
        n=$(grep -c 'CHAIN [a-z0-9-]* done' $R/logs/chain.log)
        $N "Run 2: ${BASH_REMATCH[1]} done" "All three levels finished ($n of $ndom domains)." 3
      elif [[ $l =~ CHAIN\ ([a-z0-9-]+)\ FAILED ]] && ! draining; then
        if retry_queued; then $N "Run 2: ${BASH_REMATCH[1]} FAILED" "Will be retried automatically once the other chains finish." 4
        else $N "Run 2: ${BASH_REMATCH[1]} FAILED" "Domain chain failed; other domains continue." 5; fi
      elif [[ $l == *"ALL CHAINS FINISHED (final)"* ]]; then
        $N "Run 2 finished" "All domain chains have exited (90/60/30)." 4
      fi
    done
    c0=$c
    grep -qF "ALL CHAINS FINISHED (final)" $R/logs/chain.log && exit 0
  fi
  while read -r l; do
    [[ -z $l || -n ${seen[$l]} ]] && continue
    seen[$l]=1
    draining && continue
    read -r _ _ lvl d _ <<< "$l"
    if [[ $l == *"falling back"* ]]; then
      $N "Run 2: $d fallback" "$lvl $d step 3 was killed; rerunning in linear mode." 4
    else
      ex=$(grep -o 'exit [0-9]*' <<< "$l"); st=$(grep -o 'step [0-9a-z]*' <<< "$l")
      tag=""; (( oom_new )) && tag=" Likely OOM: the kernel killed a process just before."
      [[ $ex == "exit 137" ]] && tag=" Exit 137 = SIGKILL, likely OOM."
      retry_queued && [[ $l == *"FAILED step 6 chunk"* ]] && continue   # lost chunks; the chain FAILED push says it will be retried
      $N "Run 2: $d $lvl FAILED" "${st:-a step} failed${ex:+ ($ex)}.$tag" 5
    fi
  done < <(grep -h "FAILED\|falling back" $R/L90/logs/progress.log $R/L60/logs/progress.log $R/L30/logs/progress.log 2>/dev/null)
  if draining && (( drained_sent == 0 )) && ! busy; then
    $N "Run 2 drained" "Only amidase's long chunks still run (about a day from done); safe to stop now." 4; drained_sent=1
  fi
  # /dev/shm nearly full: Logan chunks' temp files fail with "No space left on device" (10-27 19:52: 34 chunks lost)
  sf=$(df -BG --output=avail /dev/shm | tail -1 | tr -dc 0-9)
  if (( sf < 15 && shm_sent == 0 )); then $N "watson /dev/shm nearly full" "${sf} GB left; Run 2 Logan chunks may fail." 5; shm_sent=1; fi
  (( sf > 60 )) && shm_sent=0
  # Logan chunks stopped by shm_guard.sh for over 2 h: the guard is meant to clear this itself (10-29: 22 h lost silently)
  st=0; for p in $(pgrep -u $USER -x diamond); do
    [[ $(awk '{print $3}' /proc/$p/stat 2>/dev/null) == T && $(readlink /proc/$p/cwd 2>/dev/null) == $R/* ]] && st=$((st + 1)); done
  [[ -e $R/PAUSED ]] && st=0   # deliberate pause (pause_run2.sh): no push
  if (( st == 0 )); then stall_since=0; stall_sent=0
  elif (( stall_since == 0 )); then stall_since=$(date +%s)
  elif (( stall_sent == 0 && $(date +%s) - stall_since > ${STALL_LIMIT:-7200} )); then
    $N "Run 2 chunks paused" "$st Logan chunk(s) paused for over 2 hours." 5; stall_sent=1; fi
  a=$(df -BG --output=avail /watson | tail -1 | tr -dc 0-9)
  if (( a < 800 && low_sent == 0 )); then $N "watson disk low" "/watson has ${a} GB free." 5; low_sent=1; fi
  (( a > 1000 )) && low_sent=0
  if ! pgrep -f "run2/launch.sh|\./launch.sh|followup_run2\.sh" > /dev/null; then
    $N "Run 2 launcher stopped" "The launcher is no longer running and not all domains reported done." 5; exit 1
  fi
  sleep ${WATCH_SLEEP:-60}
done
