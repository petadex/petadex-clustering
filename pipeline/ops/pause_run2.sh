#!/bin/bash
# Pause Run 2 WITHOUT losing work (2026-10-30): SIGSTOP every Run 2 DIAMOND, its gate wrappers and the xargs that starts
# chunks, so no CPU is used and nothing new starts. Memory (RSS) and /dev/shm stay held by the frozen chunks. Stops
# shm_guard (it would resume them) and tells the watcher not to push "chunks paused". Undo with ./unpause_run2.sh.
R=/watson/rnalab/dzhu/petadex-clustering/run2; cd $R
declare -A anc; p=$$; while (( p > 1 )); do anc[$p]=1; p=$(awk '{print $4}' /proc/$p/stat 2>/dev/null || echo 1); done
touch PAUSED
for p in $(pgrep -u $USER -x shm_guard.sh); do kill $p; done
sleep 1
sel() { local p c; for p in $(pgrep -u $USER -x diamond) $(pgrep -u $USER -x xargs) $(pgrep -u $USER -f 'run2/bin/diamond'); do
          [[ -n ${anc[$p]} ]] && continue
          c=$(readlink /proc/$p/cwd 2>/dev/null) || continue
          [[ $c == $R || $c == $R/* ]] && echo $p; done | sort -u; }
pids=$(sel); kill -STOP $pids 2>/dev/null; sleep 2
echo "[$(date '+%F %T')] PAUSE: pause_run2.sh froze $(wc -w <<< "$pids") Run 2 processes (SIGSTOP, nothing killed, memory and /dev/shm still held); shm_guard stopped; resume with unpause_run2.sh" >> logs/chain.log
n=0; for p in $pids; do [[ $(awk '{print $3}' /proc/$p/stat 2>/dev/null) == T ]] && n=$((n+1)); done
echo "frozen $n of $(wc -w <<< "$pids")"; echo "still running: $(for p in $(sel); do [[ $(awk '{print $3}' /proc/$p/stat) != T ]] && echo -n "$p "; done)"
