# Benjamin Buchfink's `unid/` clustering pipeline, as used by Run 2

Distilled 2026-11-04 from the upstream clone `/watson/rnalab/dzhu/clustering/`
(`https://github.com/bbuchfink/clustering`, commit `ccb98d1`, 2026-10-15; no README upstream).
Run 2 vendors the scripts in `run2/unid/` (commit recorded in `run2/unid/UPSTREAM_COMMIT`).
This note is the reference for what each script does and how Run 2 differs. Run 1 (`docs/run_history.md`,
`archive/run1_scripts/cluster_pipeline.sh`) is a different, self-written design and does not use these scripts.

## What the pipeline does (per CATH domain, per identity level)

One domain at a time, three sources in priority order **PAZy → nr → Logan**. The centroid of a cluster is always
from the highest-priority source present. All DIAMOND calls use `--ext full --comp-based-stats 0 --masking 0`;
coverage is one-directional (**member** covered ≥ 90 %); `$id_cutoff` is the level's identity (90 / 60 / 30).

| Step | Script | What it does | Key DIAMOND settings |
|---|---|---|---|
| 1 | `1-cluster-pazy.sh` → `cluster-all-vs-all.sh` | PAZy cores: all-vs-all alignment, then greedy vertex cover → `pazy_<d>.clust`, `pazy_<d>.reps.faa` | `blastp -k0 --ultra-sensitive --id $id --query-or-subject-cover 90`; `greedy-vertex-cover --member-cover 90` |
| 2 | `2-map-nr-vs-pazy.sh` | Map every nr core onto the PAZy reps; the unmapped remainder goes on → `nr_vs_pazy_<d>.aln`, `.unmapped.faa` | `blastp -k1 -c1 --ultra-sensitive --id $id --query-cover 90` |
| 3 | `3-cluster-nr.sh` → `cluster-full.sh` | Cluster the unmapped nr cores → `nr_<d>.clust`; PAZy reps + nr reps = `pazy_nr_<d>.reps.faa` | `deepclust --cluster-steps $steps_nr --id $id --member-cover 90 -M 2000G` |
| 4 | `4-mapback-logan-lin.sh` | Map all Logan cores onto `pazy_nr` reps with a fast linear search → `logan_vs_nr_<d>.aln`, `.unmapped.faa` | `blastp --fast --linsearch -k1 -c1 --id $id --query-cover 90` |
| 5 | `5-split.sh` / `split-fasta` | Split the Logan leftovers into chunks by residue count (per-domain sizes in `cfg.sh`) | — |
| 6 | `6.sh` (array job, one chunk per task) | Second, more sensitive Logan map-back of each chunk onto `pazy_nr` reps → `split-<d>/out_<i>` | `blastp $sens_backmap_logan -k1 -c1 -b1 --id $id --query-cover 90 -t /dev/shm` |
| 7 | `7-get-unmapped2.sh` | Concatenate chunk outputs → `logan_vs_nr_<d>.aln2`; unmapped again → `.unmapped2.faa` | — |
| 8 | `8-linclust.sh` | Cluster the still-unmapped Logan cores → `logan_vs_nr_<d>.unmapped2.clust` | `deepclust --cluster-steps $steps_logan_linclust --id $id --member-cover 90` |
| 9–13 | `9-split`, `10.sh`, `11-concat`, `12-gvc`, `13-merge` | Self-align the step-8 Logan reps and greedy-vertex-cover them into `…unmapped2.clust2` / `clust_merged` | **Dead code since upstream commit `71a6e9e`**: step 14 reads `unmapped2.clust`, not `clust_merged`. Run 2 skips them |
| 14 | `14-merge-all.sh` | Concatenate all six assignment files, `LC_ALL=C` sort on column 1 → **`<d>.all.tsv`** | — |
| 15 | `15-get-reps.sh` | Pull the centroid sequences per source into `reps/{pazy,nr,logan}/` (next level's input) | `seqtk subseq` |

### Output format

`<d>.all.tsv`: `representative<TAB>member`, one row per core, including each centroid's self row.
It is the concatenation of: `pazy_<d>.clust`, `nr_vs_pazy_<d>.aln`, `nr_<d>.clust`, `logan_vs_nr_<d>.aln`,
`logan_vs_nr_<d>.aln2`, `logan_vs_nr_<d>.unmapped2.clust`. (The `-f 6 sseqid qseqid` alignment files already have
rep, member order.)

Reps come from exactly three places, which is why a centroid is never a lower-priority core than a member:
`pazy_<d>.clust`, `nr_<d>.clust` and `unmapped2.clust`. An nr core mapped onto a PAZy rep, or a Logan core mapped onto
a PAZy/nr rep, is never a centroid.

### Levels chain through `reps/`

| Level dir | Input (`data_dir`) | `id_cutoff` | `steps_nr` (step 3) | `sens_backmap_logan` (step 6) | `steps_logan_linclust` (step 8) |
|---|---|---|---|---|---|
| L90 | all cores (`input90/`) | 90 | `default` | `--fast` | `faster_lin fast_lin shapes-30x10_lin` |
| L60 | `L90/reps` | 60 | `default` | (DIAMOND default) | `faster_lin default` |
| L30 | `L60/reps` | 30 | `ultra-sensitive` | `--ultra-sensitive` | `ultra-sensitive` |

Upstream also has `unid90_30/` (90 → 30 directly: `id 30`, `default ultra-sensitive`, `--ultra-sensitive`,
`faster_lin default ultra-sensitive`), which Run 2 mirrors as the `control30/` comparison.
The L60 and L30 inputs are the previous level's centroids, so each level is a clustering of centroids, and
nesting is exact by construction (verified 30 of 30, see `analysis/uniqueness_check.sh`, reports in `data/uniqueness/`).

## Helpers and dependencies

- `fasta2tsv.sh`: `bioawk -c fastx '{print $name"\t"$seq}'`. `sort.sh`: `LC_ALL=C sort --parallel 16 -S 20% -T . --batch-size=1021`.
  `merge.sh`: `join` of two cluster files to compose them (used only by the dead steps 9–13).
- `split-fasta`: a compiled binary, **no source in the repo**. Usage `split-fasta <fasta> <max residues per chunk> <outdir> <prefix>`
  → `<outdir>/<prefix>_<i>.fasta`; it needs the chunk size ≥ the longest record.
- Needs `diamond` (2.2.4 in Run 2), `seqtk`, `bioawk`. Conda env `unid-clustering`.
- `domains`: one domain name per line; `6.sh` and `10.sh` pick theirs by line number (`sed "$1!d" ../domains`).
  Upstream lists 11 domains; Run 2 uses 15. Shared: abhydrolase, amidase-signature, betalactamase, trypsin-like,
  cupredoxins, cytochromec-like, d-ala-esterase, oxygenase-like. Upstream only: arylesterase, cyclase-like, uridase-like.
  Run 2 only: sgnh-hydrolase, peptidase-s8-s53, gh-deacetylase, peroxidase, cytochrome-p450, acid-protease, metallo-phosphatase.

## Upstream quirks to know about

- **`unid90/cfg.sh` has `id_cutoff=60`**, not 90 (known issue CL-09). Run 2's `L90/cfg.sh` restores 90. Anyone re-running
  upstream's `unid90` as is would produce a 60 % clustering.
- Paths are hard-coded to `/root/tmp/clustering/...` (Benjamin's machine). Run 2 rewrote them to `run2/`.
- `2-map-nr-vs-pazy.sh` and `3-cluster-nr.sh` lack `set -e`, so a killed DIAMOND is masked by the trailing `rm`/`cat`.
  Run 2 added `set -e`.
- `-M 2000G` is Benjamin's RAM; watson has 754 GB, so Run 2's `bin/diamond` caps it at `$DIAMOND_MEM`.
- `10.sh` hard-codes `--id 90 --id2 17` (dead step, irrelevant).
- Step 5/9 split sizes are per-domain in `cfg.sh` (`split_amidase`, `split_abhydrolase`, `split_arylesterase`); they size an
  array job and are overridden in Run 2 (below).

## Run 2 deviations (none change what is aligned or the identity/coverage settings)

Full detail is in the header of `run2/run_unit.sh`, `run2/bin/diamond` and `run2/deviations.tsv`.

| Change | Why |
|---|---|
| `set -e` in steps 2 and 3; `cluster-full.sh` honours `STEPS_NR_OVERRIDE` | see above; allows the fallback below |
| Step 6 runs its chunks in parallel via `xargs -P $P6` with `$T6` threads each, split to ~`$NCHUNK` chunks by residue count | replaces the upstream array job |
| Steps 9–13 skipped | dead output (above) |
| Resumable: `.done/<unit>.<step>` markers; a step-6 chunk is finished when its log ends with DIAMOND's `Peak RSS` line | survives kills |
| Non-linear `deepclust` (step 3 always; step 8 when its steps are not all `*_lin`) serialised under `logs/.heavy.lock` with `$HEAVY_MEM`/`$HEAVY_THREADS` | four in parallel spilled > 1 TB of temp files on 2026-10-26 |
| Step 3 fallback to `--cluster-steps faster_lin fast_lin` if killed (exit 137/143) or forced via `FORCE_LINEAR_NR` | recorded in `deviations.tsv`; **only L90 amidase-signature step 3 (2026-10-27)** used it |
| `bin/diamond` wrapper: thread count from `$DIAMOND_THREADS`, `-M` capped, start gate for `-t /dev/shm` chunks, drain flag | watson resource limits |
| Final per-unit check: input rows = `all.tsv` rows = distinct members; clusters = reps | catches lost cores |

Run-time guards (`disk_guard.sh`, `shm_guard.sh`, `notify_watch.sh`) and the pause/stop/resume scripts are operational
only and do not alter the clustering.

## Where things are

| Item | Path |
|---|---|
| Upstream repo (local clone deleted 2026-11-04) | `https://github.com/bbuchfink/clustering` |
| Vendored scripts used by Run 2 | `run2/unid/` (copy in `pipeline/unid/`) |
| Per-level configs | `run2/L*/cfg.sh` (copies in `pipeline/configs/`) |
| Per-domain driver | `run2/run_unit.sh`, `run2/run_chain.sh` (copies in `pipeline/`) |
| Results | `run2/L*/<domain>.all.tsv`, centroids in `run2/L*/reps/` |
| Uniqueness / nesting check | `analysis/uniqueness_check.sh` |
