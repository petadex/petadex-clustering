#!/bin/bash
# One CATH domain through all three levels: 90 -> 60 -> 30. Each level resumes from its .done markers.
set -euo pipefail
cd /watson/rnalab/dzhu/petadex-clustering/run2
for L in L90 L60 L30; do ./run_unit.sh $L $1; done
