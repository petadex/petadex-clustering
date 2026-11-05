#!/bin/bash
# Undo pause_run2.sh: SIGCONT the frozen Run 2 processes, clear PAUSED, restart shm_guard.
R=/watson/rnalab/dzhu/petadex-clustering/run2; cd $R
declare -A anc; p=$$; while (( p > 1 )); do anc[$p]=1; p=$(awk '{print $4}' /proc/$p/stat 2>/dev/null || echo 1); done
sel() { local p c; for p in $(pgrep -u $USER -x diamond) $(pgrep -u $USER -x xargs) $(pgrep -u $USER -f 'run2/bin/diamond'); do
          [[ -n ${anc[$p]} ]] && continue
          c=$(readlink /proc/$p/cwd 2>/dev/null) || continue
          [[ $c == $R || $c == $R/* ]] && echo $p; done | sort -u; }
rm -f PAUSED
pgrep -u $USER -x shm_guard.sh > /dev/null || { nohup setsid nice -n 10 ./shm_guard.sh > logs/shm_guard.out 2>&1 < /dev/null & }
for p in $(sel); do kill -CONT $p; sleep 0.5; done   # staggered: the guard re-pauses if /dev/shm dips
echo "[$(date '+%F %T')] UNPAUSE: unpause_run2.sh resumed Run 2" >> logs/chain.log
echo resumed
