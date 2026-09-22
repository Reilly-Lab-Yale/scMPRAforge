# Synthetic-factorial sweep history

The synthetic-factorial design space spans seven axes that jointly determine
the statistical setup of an scMPRA experiment: `n_cells`, `n_cres`,
`bcs_per_cre`, `moi`, `lib_alpha_nb` (NB overdispersion), `minP` (per-CRE
floor activity), and `activity_max_mult` (dynamic range, p95 over minP).

Detection power is characterized from a single 5000-point Latin Hypercube
sample (the **union** sweep) over the full design box. The three earlier
sweeps (full, topup, topup3) were exploratory: they calibrated where the box
boundaries needed to be so the empirical datasets (cohen, shendure, takeshi)
land inside it. Their simulations have since been removed; the union sweep
covers a superset of their ranges on every axis, so nothing is lost in
coverage. This file records that provenance, since the leftover `s`/`t`/`u`
sample-id prefixes and the box-calibration history may otherwise be opaque.

Each union sample is tested under eight arms rather than one, so the sweep
also resolves how the *transfection reporter* interacts with the seven design
axes. See "Test arms" below.

## Sweeps

| Sweep | n_pts | LHS seed | sid prefix | Box | Role |
|---|---|---|---|---|---|
| full   | 1000 | 20260505 | `s` | `AXIS_BOUNDS`        | Exploratory: initial box, calibrated to expected empirical bounds |
| topup  | 100  | 20260506 | `t` | `TOPUP_AXIS_BOUNDS`  | Exploratory: cohen-corner extension (high `bcs_per_cre`, high `moi`) |
| topup3 | 120  | 20260507 | `u` | `TOPUP3_AXIS_BOUNDS` | Exploratory: wider `moi`, `minP`, `activity_max_mult` to reach takeshi |
| **union** | **5000** | 20260511 | `v` | `UNION_AXIS_BOUNDS` | **Canonical: independent LHS over the full union box** |

The union box is the superset of all three exploratory boxes on every axis
(`bcs_per_cre` 3-50k, `moi` 0.5-350, `minP` 0.003-2.0, `activity_max_mult`
1-120; the other three axes were unchanged throughout). So union coverage
strictly contains the earlier sweeps.

## Test arms

Every simulation is tested eight ways. All eight read the same counts, so the
arms are exactly paired within a sample: any between-arm difference is
attributable to the test, never to a difference in simulated data.

| Arm | Test | Reporter |
|---|---|---|
| `mwu` / `mwu_deflated`               | Mann-Whitney U | present / absent |
| `ttest` / `ttest_deflated`           | Welch's t      | present / absent |
| `ks` / `ks_deflated`                 | two-sample KS  | present / absent |
| `pseudobulk` / `pseudobulk_deflated` | Welch's t on per-replicate means | present / absent |

`mwu` is the primary arm and the one the design-space figures use; Welch's t
over-rejects on scMPRA counts by 2-3x. The other six ride along because they
cost only a second pass over counts that already exist, whereas recovering
them later would cost a full re-simulation -- the raw sims are pruned as the
sweep goes and cannot be re-tested after the fact.

### The reporter contrast

A `_deflated` arm drops zero-count observations before testing. That is what a
no-reporter experiment actually sees: the simulator emits one row per
transfection event, so a zero in the raw counts is a transfected-but-silent
event, and only a transfection reporter can tell that apart from a cell that
was never transfected at all. Without the reporter those rows are not zeros in
the data -- they are absent from it.

A CRE that loses every observation to the drop yields a NaN p-value, which BH
adjustment fills to 1.0. Such CREs therefore count as non-detections rather
than disappearing from the denominator: `fc` comes from the ground truth, not
from the test, so the row survives aggregation. This is the conservative
treatment and the correct one -- a CRE you never observe is a CRE you cannot
call. It is mildly *over*-conservative in one respect: the filled rows still
occupy the BH denominator, whereas a real no-reporter experiment would never
have tested them.

Note that this contrast lives entirely at test time. It is a different (and
much cheaper) question than the fit-time reporter evaluation sketched in
`paper_plan.md`, which contrasts an obs-condition NB fit against a
consider-missing + MOIB fit under Wald. Ortho fitting is not feasible at
125,000 simulations, so the sweep answers the reporter question in the
nonparametric frame only.

## Why a single independent union sweep (not a pooled set)

Each individual sweep is a proper Latin Hypercube (`scipy.stats.qmc`), but
concatenating sweeps is **not** itself an LHS. The exploratory sweeps each
occupy a different sub-box, so a pooled dataset acquires structural
cross-axis correlation - e.g. a high-`moi` sample could only have come from
topup/topup3, which also biased its `bcs_per_cre` and `minP`. That confounds
per-axis attribution (a LOESS marginal over a correlated design implicitly
absorbs the co-moving variables).

| Dataset | n | mean &#124;r&#124; | max &#124;r&#124; |
|---|---|---|---|
| full only            | 1000 | 0.033 | 0.104 |
| combined patchwork (full+topup+topup3) | 1220 | 0.106 | 0.444 |
| **union (canonical)** | **5000** | **0.009** | **0.022** |
| pooled (all four)    | 6220 | 0.027 | 0.094 |

Worst-pair entanglement in the patchwork: `bcs_per_cre x moi` r = +0.44,
`moi x minP` r = -0.41. The union sweep, drawn from scratch over the full
box, has cross-axis &#124;r&#124; at the level expected from sampling noise
alone (max 0.022). A pooled set (union + the three exploratory sweeps) would
still be acceptable (max 0.094, since union is 80% of it), but it trades the
cleanest design for ~24% more points whose only added density is in the
corners - density that a `frac=0.4` LOESS smoother barely uses. We therefore
analyze the union sweep alone.

## Figures

All figures derive from the union sweep's aggregated power summary,
`output/samples_power_union.parquet` (5000 rows, MWU power):

- `marginals*_union*.svg` - per-axis marginals, 5 metrics x {no/with takeshi}
  x {log, linear x}.
- `pairwise_heatmaps_union.svg` - top-3-axis pairwise power heatmaps.
- `attribution_{bar,marginals}_<pair>.svg` - per-axis attribution for the
  three anchor pairs (`cohen_vs_shendure`, `takeshi_vs_shendure`,
  `takeshi_vs_cohen`).

The raw per-rep results behind the summary (5000 samples x 5 reps x 5 sims,
64.8M rows) are consolidated into a single parquet under the shared data area
for any future re-aggregation with different metrics:
`/nfs/roberts/project/pi_skr2/shared/tabula_data/simulated/synthetic_factorial_union_2026-05-16/union_cached_results.parquet`
(338M; see the README there for schema + provenance).

## Running the sweep

The sweep is paced to run over weeks at low concurrency rather than saturating
the cluster for a few days. 200 slices of ~25 samples each; every slice driver
stands up one dask cluster and walks its samples serially, keeping worker
submissions far below the YCRC 200/hr sbatch limit.

```bash
sbatch --array=0-199%20 wrap_union_slow.sh      # pace set by %N, nothing else
ls <cache_root>/union/*.parquet | wc -l         # progress, out of 5000
```

Interrupting is free and resuming is automatic. A sample's cache is written
only after all its reps and all eight arms are complete, and it is written to a
temp name and renamed, so a driver killed mid-write cannot leave a truncated
file that a later run would trust. Any sample with a cache is skipped for the
price of one stat, so the array can be cancelled and resubmitted at any time --
to change the pace, after a node failure, or after a 7-day timeout.

Caches live on **project**, not scratch, because scratch purges files untouched
for 30 days and a multi-week sweep would otherwise lose its early samples
before its late ones finished. The raw sims stay on scratch and are deleted per
sample as soon as that sample is cached, so peak disk scales with the number of
samples in flight, not with the size of the sweep.

## Reproducing the figures

```bash
# Marginals + pairwise heatmap from the committed power summary (no re-sim):
sbatch wrap_synthetic_factorial.sh replot union

# Per-axis attribution (3 anchor pairs):
sbatch wrap_attribution.sh
```
