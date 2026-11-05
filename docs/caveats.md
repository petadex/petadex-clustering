# PETadex catalytic-core clustering (v2): caveats

Things to know before quoting or comparing the v2 clustering numbers. The analysis itself is in
`petadex_clustering_experiment.ipynb`; this file collects the caveats and open issues separately.

## Input

### PL559 is missing from the input, probably by mistake (open, PET-25)
PL559 (AII11221.1, cupredoxin, component C015) has two PETadex HMM hits. Its best hit (e = 9.0e-27, residues 103–229) is
incomplete; its second hit (e = 8.7e-20, residues 217–344) covers all four catalytic histidines. Thomas's 2026-10-18
rewrite keeps only each protein's single best hit, and only if that hit is complete, so PL559 drops out. A "best
*complete* hit" rule would keep it (487 Plasticome proteins instead of 486). Raised with Thomas; not resolved.

### The Plasticome input changed between Run 1 and v2
Run 1 (2026-10-15) clustered 503 PAZy cores from 487 proteins. Thomas regenerated the PAZy FASTAs on 2026-10-18 (one core
per protein), giving 486 cores from 486 proteins, and v2 used that version. The 17-core difference is the 16 extra cores of
proteins that had more than one complete hit, plus PL559. nr and Logan inputs are identical between the two runs.
The display output in Thomas's `7_petadex_clustering.ipynb` still shows the Oct 11 numbers (503 / 487).

## Comparison with v1.1

### v2 has far more centroids than v1.1 at every level
| Level | v2 | v1.1 | v2 / v1.1 |
|---|---:|---:|---:|
| 90 % | 32,022,716 | 18,172,960 | 1.76× |
| 60 % | 3,230,662 | 1,814,296 | 1.78× |
| 30 % | 65,769 | 22,235 | 2.96× |

The same pipeline design should not produce nearly three times as many superfamilies from a similar number of sequences.
The excess is not explained by the extra domains (the 8 domains both runs share show the same pattern) and is only partly
singletons. It is unexplained.

### The two runs are not like-for-like
The inputs and parts of the method differ, so the v1.1 numbers are context, not a control:
- **What is clustered:** v1.1 clustered **full ORFs**, v2 clusters **catalytic cores** (the HMM hit range only). The v1.1 30 %
  centroids average 368 aa (median 314) against 274 aa for v2's; 21 % of them start with M. Source: `archive/compare_old_c30/compare.py` (copy in `data/v11/`)
  ("old ones (full ORFs, clustered 2026-05-26)") and `archive/compare_old_c30/old.fa`. With the same 90 % member-coverage rule, a core
  and its full-length protein can cluster differently, so this alone can change the cluster counts.
- **Inputs:** v1.1 had 211 PAZy / 4.72 M nr / 302.4 M Logan sequences; v2 has 486 / 6.76 M / 291.9 M cores, defined by
  different HMMs.
- **Domains:** v1.1 covered 11 CATH domains, v2 covers 15 (8 shared, 3 only in v1.1, 7 only in v2).
- **Settings:** v2 runs Benjamin's `unid/` scripts with the stated settings. The config that actually produced v1.1 is not
  known; the upstream `unid90/cfg.sh` ships with `id_cutoff=60` (issue CL-09), so the published configs are not reliable
  evidence of what v1.1 ran.
- **Source of the v1.1 numbers:** `archive/compare_old_c30/old.tsv` (Benjamin's v1.1 output), tabulated in
  `data/v11/singletons_by_domain_vs_v1.1.tsv` and `data/v11/source_breakdown_vs_benjamin.tsv`.

## Method

### Linear clustering for one step (amidase-signature, 90 %, nr)
Step 3 for amidase-signature at 90 % (clustering the 818k nr cores not assigned to PAZy) ran with
`--cluster-steps faster_lin fast_lin` instead of the default, after the default mode was stopped twice for filling the disk
(`pipeline/configs/deviations.tsv`). The identity and coverage thresholds are unchanged. On a 1M-sequence test in Run 1, linear mode gave
2.9 % more clusters than default, so this step probably adds a few thousand extra nr centroids (amidase has 280,236 nr
centroids at 90 %, out of 32.0 M clusters overall). Negligible for the totals.

### Compute figures exclude lost work
The 33,200 thread-hours reported in the notebook are the DIAMOND runs whose output was kept. Runs that were killed and
redone (two default-mode amidase clusterings, Logan chunks lost to /dev/shm exhaustion and to the 10-30 stop) are not
counted, so the true cost was higher.

## Unexplained observations
- **Trypsin-like superfamilies:** trypsin-like is 3.2 % of 90 % clusters but 20.0 % of 30 % clusters (13,186), second only
  to α/β hydrolases, and has the highest 30 % singleton rate of any domain (35.5 %). Not investigated.
