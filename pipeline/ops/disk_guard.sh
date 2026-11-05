#!/bin/bash
# Run 2 disk guard (start with nohup setsid; launch.sh does this). Every 20 s checks /watson free space:
#   < 400 GB : SIGTERM the heavy deepclust (the only DIAMOND with -M $HEAVY_MEM); its unlinked temp files are freed
#              at exit. run_unit.sh then fails that step (abhydrolase L90 step 3 falls back to linear mode instead).
#   still < 400 GB 60 s later, or < 250 GB at any time : stop the whole run (launcher, chains, every Run 2 DIAMOND).
# Exits when the launcher is gone and no Run 2 DIAMOND is left.
R=/watson/rnalab/dzhu/petadex-clustering/run2
N=~/bin/notify.sh
: ${HEAVY_MEM:=580G} ${LOW:=400} ${FLOOR:=250}
log() { echo "[$(date '+%F %T')] DISK GUARD $*" | tee -a $R/logs/chain.log; }
free_gb() { df -BG --output=avail /watson | tail -1 | tr -dc 0-9; }
run2_diamonds() { for p in $(pgrep -u $USER -f '^/watson/rnalab/dzhu/bin/diamond '); do
                    [[ $(readlink /proc/$p/cwd 2>/dev/null) == $R/* ]] && echo $p; done; }
stop_all() {
  log "$1: stopping the whole run"; $N "Run 2 STOPPED by disk guard" "$1" 5
  pkill -TERM -u $USER -f 'run2/launch.sh|\./launch.sh'; pkill -TERM -u $USER -f 'run_chain.sh|run_unit.sh'
  pkill -TERM -u $USER -f 'xargs -P [0-9]+ -I\{\} bash -c \./run_chain.sh'
  kill -TERM $(run2_diamonds) 2>/dev/null; sleep 30; kill -KILL $(run2_diamonds) 2>/dev/null
  exit 2
}
log "started (low $LOW GB, floor $FLOOR GB, heavy -M $HEAVY_MEM), free $(free_gb) GB"
killed_at=0
while true; do
  a=$(free_gb)
  (( a < FLOOR )) && stop_all "/watson free ${a} GB < ${FLOOR} GB"
  if (( a < LOW )); then
    if (( killed_at == 0 )); then
      h=$(pgrep -u $USER -f "^/watson/rnalab/dzhu/bin/diamond deepclust .*-M $HEAVY_MEM")
      log "/watson free ${a} GB < ${LOW} GB: SIGTERM heavy deepclust ${h:-(none running)}"
      $N "Run 2 disk guard fired" "/watson free ${a} GB; killed the heavy clustering step." 5
      [[ -n $h ]] && kill -TERM $h
      killed_at=$SECONDS
    elif (( SECONDS - killed_at >= 60 )); then
      stop_all "/watson free still ${a} GB 60 s after killing the heavy step"
    fi
  else
    killed_at=0
  fi
  # 2026-11-02: also stay while a chain/unit script runs (between steps no DIAMOND runs for a few seconds)
  if ! pgrep -u $USER -f 'run2/launch.sh|\./launch.sh|run_chain\.sh|run_unit\.sh' > /dev/null && [[ -z $(run2_diamonds) ]]; then
    log "launcher gone and no Run 2 DIAMOND left: exiting"; exit 0
  fi
  sleep 20
done
