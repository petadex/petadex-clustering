#!/bin/bash
# Stop Run 2 completely (2026-10-27, to free watson for Thomas's priority DIAMOND work). Safest after the "Run 2
# drained" push (touch DRAIN first): then nothing is in flight. Finished steps (.done markers) and finished Logan chunks
# are kept either way; ./resume_run2.sh picks up from there. Any chunk or step still running is lost.
# Selects Run 2 processes by script/binary name AND cwd inside run2, never this script's own ancestors (so the shell
# you run it from survives); freezes them all (SIGSTOP) so nothing reacts, then SIGKILLs them.
R=/watson/rnalab/dzhu/petadex-clustering/run2
declare -A anc; p=$$; while (( p > 1 )); do anc[$p]=1; p=$(awk '{print $4}' /proc/$p/stat 2>/dev/null || echo 1); done
pat='launch\.sh|run_chain\.sh|run_unit\.sh|requeue_amidase\.sh|rebalance_amidase\.sh|disk_guard\.sh|notify_watch\.sh|unid/|bin/diamond|xargs -P|^diamond '
run2_pids() { local p c; for p in $(pgrep -u $USER -f "$pat") $(pgrep -u $USER -x diamond); do
                [[ -n ${anc[$p]} ]] && continue
                c=$(readlink /proc/$p/cwd 2>/dev/null) || continue
                [[ $c == $R || $c == $R/* ]] && echo $p; done | sort -u; }
pids=$(run2_pids); n=$(wc -w <<< "$pids")
(( n == 0 )) && { echo "no Run 2 processes running"; exit 0; }
echo "[$(date '+%F %T')] STOP: stop_run2.sh killing $n Run 2 processes (priority DIAMOND work); resume with resume_run2.sh" >> $R/logs/chain.log
kill -STOP $pids 2>/dev/null; kill -KILL $pids 2>/dev/null
sleep 5; left=$(run2_pids)
[[ -n $left ]] && { kill -KILL $left 2>/dev/null; sleep 5; left=$(run2_pids); }
echo "killed $n processes; still running: ${left:-none}"
free -g | sed -n 2p; df -h /dev/shm | tail -1
