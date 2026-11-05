# Run history: Run 1 (per CATH component, own settings), 2026-10-14 to 10-15

Before v2 (Benjamin's `unid/` pipeline, run as "Run 2" in `run2/`), the cores were clustered with a self-written pipeline that
reproduced the *design* of Benjamin's earlier clustering. This is that run's original README, kept for the record. Its code, outputs
and checks are in `archive/` (local only, not in the repo): paths below that start `scripts/`, `results/`, `preflight/`, `test/` or
`work/` are now under `archive/` (`scripts/` is `archive/run1_scripts/`). Its numbers differ from v2's: different method, per-component
instead of per-domain, and the PAZy input was 503 cores before Thomas's 2026-10-18 rewrite (see `docs/caveats.md`).

Summary of how the runs relate:

| Run | Grouping | Method | 90 % / 60 % / 30 % clusters |
|---|---|---|---|
| v1.1 (Benjamin, 2026-05) | 11 CATH domains | `unid/`, full ORFs, older input | 18,172,960 / 1,814,296 / 22,235 |
| Run 1 (2026-10-15) | 46 CATH components | own pipeline (this file) | 30,453,501 / 2,928,810 / 167,232 |
| unid 30 % step only (2026-10, `archive/unid30_run/`) | 15 CATH domains | Run 1's 60 % centroids → `unid30` | – / – / 77,524 |
| **v2 = Run 2 (2026-11-04)** | 15 CATH domains | `unid/` end to end | **32,022,716 / 3,230,662 / 65,769** |

---

## Original Run 1 README

Priority clustering of every complete catalytic core in PETadex into 90 / 60 / 30 % identity
centroids, run on watson with DIAMOND 2.2.4. Started 2026-10-14.

## Input

Thomas Quigley's per-component catalytic-domain FASTAs from notebook
`petadex.v260701/dataset_generation/7_petadex_clustering.ipynb` (2026-10-11):

`/watson/rnalab/tquigley/petadex.v260701/dataset_generation/7_petadex_clustering/output/`

| Source | Directory | Domains |
|---|---|---|
| PAZy / plasticome (`PL…`) | `plasticome_catalytic_domain_per_component/` | 503 |
| NCBI nr (`NR…`) | `nr_catalytic_domain_per_component/` | 6,764,617 |
| Logan (`LG…`) | `logan_catalytic_domain_per_component/` | 291,922,243 |
| **Total** | 46 CATH components each | **298,687,363** |

Each record is the **catalytic core only**: the PETadex-HMM hit range
`sequence[hmm_start-1:hmm_end]` of a hit with i-evalue < 1e-5 and `complete_catalytic_hit == TRUE`
(a residue in every catalytic match state). Header: `>{petadex_id}|{hmm_name}|{hmm_evalue}|{hmm_start}|{hmm_end}|{component_id}`.
One protein can contribute more than one core, so the full header is the domain ID.

Preflight (`preflight/validate.tsv`, `scripts/validate_file.sh`), all 138 files, all 298,687,363 records:
sequence length == `hmm_end - hmm_start + 1` for every record, header component == file component,
no empty sequences, no duplicate IDs. nr cores were also spot-checked against the full proteins in
`petadex-nr.v260701.proteins.fastaa`: each is the exact substring at the HMM offset.

## Method

Each component is clustered on its own (sequences cannot cluster across CATH components).
At each level the three sources are processed in priority order **PAZy → nr → Logan**:

1. Tier *k* is aligned to the centroids already chosen (`diamond blastp -k1 --max-hsps 1 --approx-id L --query-cover 90`);
   a hit assigns it to that centroid.
2. The unassigned rest is clustered (`diamond cluster --approx-id L --member-cover 90`), and its centroids join the centroid set.

So a centroid always comes from the highest-priority source in its cluster, which keeps metadata-rich PAZy / nr
sequences as representatives, and new sequences can later be added the same way.

| Level | Input | Assignment sensitivity |
|---|---|---|
| 90 % | all catalytic cores | `--fast` |
| 60 % | the 90 % centroids | `--sensitive` |
| 30 % | the 60 % centroids | `--more-sensitive` |
| control | the 90 % centroids straight to 30 % | `--more-sensitive` |

**Exception, C003 at 90 %:** its Logan assignment step uses `--faster -g 512` instead of `--fast`. With `--fast`,
the 17.4M C003 Logan cores against 319k near-identical nr centroids ran at ~60 queries/s (projected ~80 h), because
every query was fully aligned against hundreds of candidate centroids. `-g 512` limits full alignment to the top 512
targets per query. Against the `--fast` output for C003's first 100k Logan queries, `--faster -g 512` recovered
49,388 / 49,510 assignments (99.75 %) with no extra ones, in 36 s rather than ~28 min. The ~0.25 % it misses are
clustered among the Logan leftovers instead, so C003's 90 % count is very slightly inflated. All other components use `--fast`.

**Exception 2, C003 at 90 %, Logan leftovers:** `diamond cluster` on the 9,051,662 unassigned C003 Logan cores was
SIGKILLed (almost certainly out of memory) in its round 2, the non-linear all-vs-all `fast` search over the 3.04M
round-1 representatives. `-g` is not permitted in the cluster workflow, so this one step uses
`--cluster-steps faster_lin fast_lin` (both rounds linear). On a 1M subset of the same leftovers:
default 1,352 s / 50 GB / 476,786 clusters; linear 143 s / 5.6 GB / 490,668 clusters (+2.9 %).
**Parallel split assignment (all levels, any tier > 500k queries; used for C003 60 % and control):** on C003 a single
`blastp --sensitive` process averaged ~8.5 busy cores of 200 for 5 h, because its repetitive sequences funnel the
seed search into a few huge seed buckets that one thread each must process. `assign_tier()` now splits large tiers
round-robin into 40 chunks run as parallel DIAMOND processes with the same settings (~25–35 s per seed-search
chunk instead of ~420 s). On the 3-component test set, forcing the split changed c90/c60/c30 by −1 / −10 / −3
(18,365→18,364, 2,633→2,623, 1,480→1,477): DIAMOND's heuristics depend slightly on which queries share a block.
Components finished before 2026-10-15 10:47 ran unsplit.

So C003's Logan-origin 90 % centroids are inflated by roughly 3 %; it was finished from the saved state by
`scripts/resume_c003_L90.sh` (driver `scripts/finish_c003.sh`). All other cluster steps use DIAMOND's default cascade.

Coverage is uni-directional (member covered ≥ 90 %), as agreed with Benjamin Buchfink for his earlier run.
`--approx-id` is DIAMOND's approximate identity, the measure `diamond cluster` itself uses.

This reproduces the *design* of Benjamin's earlier clustering, not his scripts (those are in his private GitHub repo),
on a newer and different input, so cluster counts are not expected to match his
18,172,960 / 1,814,296 / 22,235.

## Products (`results/`)

| File | Content |
|---|---|
| `centroids_90.faa`, `centroids_60.faa`, `centroids_30.faa` | centroid catalytic cores at each level |
| `hierarchy.tsv` | one row per input domain: `domain_id petadex_id source component c90 c60 c30` |
| `map_90_to_60.tsv`, `map_60_to_30.tsv` | centroid → parent centroid |
| `cluster_stats_{90,60,30}.tsv` | per centroid: `component centroid_source n_children n_domains n_pazy n_nr n_logan` (`n_domains` counts all input domains beneath it) |
| `control_90to30_centroids.faa`, `map_90_to_30_direct.tsv` | control: 90 centroids clustered straight to 30 % (transitivity check) |
| `summary.tsv`, `summary_by_component.tsv` | counts per source, centroids per level, singletons |
| `step_timings.tsv` | wall time and in/out counts for every assign / cluster step |
| `MD5SUMS` | checksums |

Checks run at the end (`scripts/build_products.sh`): each level places every input sequence exactly once,
`hierarchy.tsv` has all 298,687,363 domains, centroid FASTAs match the stats tables.

### Notes on the outputs

- **FASTA line wrapping is mixed.** Sequences that passed through DIAMOND's `--un` output are wrapped at 160
  characters; the rest are single-line. Read the `.faa` files with a multi-line-aware parser. All 33.5M centroid
  sequences were checked: length == `hmm_end − hmm_start + 1` for every record at 90/60/30.
- **Result (2026-10-15):** 298,687,363 domains → 30,453,501 (90 %) → 2,928,810 (60 %) → 167,232 (30 %);
  direct 90→30 control 182,915. Benjamin's earlier run (different, older input): 18,172,960 / 1,814,296 / 22,235.

### Why the 30 % count is high (investigated 2026-10-15)

The 90→60 collapse matches Benjamin's (10.4× vs 10.0×); only 60→30 differs (17.5× vs 81.6×).

- **It is a long tail:** 16,227 of the 167,232 clusters hold 99 % of all domains; 96,037 clusters have ≤ 5 domains,
  95,182 of which contain no PAZy or nr member. 51.5 % of 30 % clusters hold a single 60 % centroid.
- **Not short fragments:** small-cluster centroids are normal length (median 223 aa, 3.9 % < 100 aa).
- **Not the tiering:** on C016, plain `diamond cluster --approx-id 30 --member-cover 90` on the 266,213 60 % centroids
  gives 12,690 clusters vs 10,788 from the tiered pipeline.
- **Member coverage is the biggest lever:** the same C016 run with `--member-cover 80` gives 7,739 (−28 % vs tiered),
  `--member-cover 50` gives 4,731 (−56 %). At 30 % identity, alignments rarely span 90 % of the member.
  Scaling that to all components still leaves ~70–75k, ~3× Benjamin's, so his 30 % step likely differed in
  other parameters too (not known: his pipeline is in his private GitHub repo).

## Re-running

```bash
cd /watson/rnalab/dzhu/petadex-clustering
T=192 MEM=400 ./scripts/cluster_pipeline.sh   # resumable: finished steps leave work/<comp>/<level>/.done
```

`ROOT`, `IN`, `N_COMP`, `EXPECT` are environment-overridable; `test/` is a 3-component dry run (C008, C027, C039).
