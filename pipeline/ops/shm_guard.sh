#!/bin/bash
# Run 2 /dev/shm guard (2026-10-29). Logan back-map chunks (6.sh, DIAMOND -t /dev/shm) keep temp files in the 378 GB
# tmpfs; when it fills they fail ("No space left on device", 34 chunks lost 10-27 19:52). Every 30 s: below $LOW GB free,
# and still falling, SIGSTOP the newest running Run 2 chunk (it keeps its files but stops growing); while free space is
# not falling, SIGCONT the oldest paused one. Lossless: a paused chunk just finishes later. Exits when followup_run2.sh
# is gone and nothing is paused. notify_watch.sh pushes separately if /dev/shm still drops below 15 GB.
# 2026-10-30: was "pause whenever below LOW, resume only above HIGH", which deadlocked: a chunk's temp files shrink as
# it aligns, so pausing every chunk at a flat 24 GB froze the only thing that could free space (10-29 16:27, 22 h lost).
# The guard now adopts chunks already stopped when it starts.
R=/watson/rnalab/dzhu/petadex-clustering/run2
: ${LOW:=25}
log() { echo "[$(date '+%F %T')] SHM GUARD $*" >> $R/logs/chain.log; }
free_gb() { df -BG --output=avail /dev/shm | tail -1 | tr -dc 0-9; }
paused=()
chunks() { for p in $(pgrep -u $USER -x diamond); do   # running Run 2 chunk DIAMONDs, oldest first
             [[ $(readlink /proc/$p/cwd 2>/dev/null) == $R/* ]] || continue
             tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null | grep -q -- '-t /dev/shm' || continue
             [[ $(awk '{print $3}' /proc/$p/stat) == T ]] && continue
             [[ $(readlink /proc/$p/exe 2>/dev/null) == *bash ]] && continue   # still waiting at the bin/diamond gate
             # amidase chunks sort last, so they are paused first: the other domains are near done (10-29 14:30)
             k=0; tr '\0' ' ' < /proc/$p/cmdline | grep -q 'split-amidase-signature/' && k=1
             echo "$k $(awk '{print $22}' /proc/$p/stat) $p"; done | sort -k1,1n -k2,2n | awk '{print $3}'; }
for p in $(pgrep -u $USER -x diamond); do   # adopt chunks a previous guard left stopped
  [[ $(readlink /proc/$p/cwd 2>/dev/null) == $R/* && $(awk '{print $3}' /proc/$p/stat) == T ]] && paused+=($p); done
log "started (pause below ${LOW} GB free while falling, resume while not falling; ${#paused[@]} adopted)"
prev=$(free_gb)
while pgrep -u $USER -f 'followup_run2\.sh' > /dev/null || (( ${#paused[@]} > 0 )); do
  live=(); for p in "${paused[@]}"; do kill -0 $p 2>/dev/null && live+=($p); done; paused=("${live[@]}")
  f=$(free_gb)
  if (( f < LOW && f < prev )); then
    p=$(chunks | tail -1)
    if [[ -n $p ]]; then kill -STOP $p; paused+=($p); log "/dev/shm ${f} GB free and falling: paused chunk pid $p (${#paused[@]} paused)"; fi
  elif (( f >= prev && ${#paused[@]} > 0 )); then
    p=${paused[0]}; kill -CONT $p; paused=("${paused[@]:1}"); log "/dev/shm ${f} GB free: resumed chunk pid $p (${#paused[@]} paused)"
  fi
  prev=$f
  sleep 30
done
log "exiting"
