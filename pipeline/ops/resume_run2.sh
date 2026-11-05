#!/bin/bash
# Resume Run 2 after stop_run2.sh: clears DRAIN, installs launch_resume.sh as launch.sh (the old one is archived),
# starts it (it starts the disk guard) and the watcher. Refuses if any Run 2 process is still running.
R=/watson/rnalab/dzhu/petadex-clustering/run2; cd $R
if pgrep -u $USER -f 'launch\.sh|run_chain\.sh|run_unit\.sh|notify_watch\.sh|disk_guard\.sh' > /dev/null; then
  echo "Run 2 processes still running; run ./stop_run2.sh first"; exit 1; fi
rm -f DRAIN
cp -p launch.sh archive_20261026/launch.$(date +%Y%m%d-%H%M%S).sh
cp launch_resume.sh launch.new && chmod +x launch.new && mv launch.new launch.sh
nohup setsid ./launch.sh > logs/launch.out 2>&1 < /dev/null &
sleep 5
nohup setsid ./notify_watch.sh > logs/notify_watch.out 2>&1 < /dev/null &
sleep 2; tail -1 logs/chain.log
