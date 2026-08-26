# Restoring sweep monitoring after a session restart

The sweep itself (array 22641062) is independent of any Claude session and
keeps running across restarts. Only the *monitoring* is session-bound and has
to be recreated. Everything below lives on project storage precisely so a
session that comes back on a different compute node does not lose it -- the
Claude scratchpad is under `/tmp`, which is node-local.

Working dir for all of this:
`<repo>/.claude/worktrees/reporter-arms-sweep/analyses/simulation/design_space/synthetic_factorial`

## 1. Recreate the progress cron (3x/day)

`CronCreate` jobs die with the session. Recreate with cron `23 8,14,20 * * *`,
recurring, and a prompt that: runs `./progress.sh`; counts failures with
`grep -hcE "FAILED v" logs/slow-22641062_*.out`; diffs against
`logs/monitor/sweep_state.txt` (format `cached=<n> failed=<n> epoch=<n>`) and
overwrites it; reports cached/5000, delta, samples/hr, projected finish, and
new-vs-total failures; and sends a PushNotification ONLY on: new failures,
zero progress while tasks run, zero tasks running with work left, scratch
above 1 TB, or completion.

Note the state file path moved from the scratchpad to `logs/monitor/`.

## 2. Restart the pending-reason sampler

    bash logs/monitor/pendwatch.sh &

Samples `squeue -t PD` for `synfac_worker` every 15 min and appends Slurm's
pending Reason to `logs/monitor/pending_reasons.log`. Needed because Slurm
clears `Reason` to `None` once a job starts, so it is only readable live.

## 3. Known state as of 2026-08-22 15:30

- 2744/5000 cached, 73 samples failed, ETA ~Aug 26.
- Failure causes, all diagnosed and fixed:
  - worker walltime 8h < slice duration -> killed slice TAILS (fixed, 3d)
  - worker queue starvation -> kills slice HEADS (mitigated, 6h wait)
- Starvation windows cluster just after midnight (00:13-00:35 twice) with
  waits of 21, 54 and 173 min. Reason is `(Priority)` -- other users' load,
  not our 128G ask, not fairshare (ours is 0.89).
- Maintenance reservation covers 276 nodes from 2026-09-15. After the sweep,
  but relevant to anything scheduled later.

## 4. At the end: the mop-up pass

Failed samples are not retried automatically. Once the array drains, re-run
only the slices that actually contain uncached samples.

Do NOT blanket-resubmit `--array=0-199`. `_make_slurm_client()` runs before
the per-sample skip check, so a slice with no work still stands up a dask
cluster, submits N_WORKERS sbatches, and can sit in `wait_for_workers` for up
to 6h. 200 slices x 5 workers is ~1000 sbatches churning fast enough to trip
the YCRC 200/hr limit -- the same limit that cost 6 of 10 slices in 2026-05.

Compute the slices that need work (slice = row index % n_slices):

    python - <<'EOF'
    import pandas as pd
    from pathlib import Path
    cache = Path('/nfs/roberts/project/pi_skr2/shared/tabula_data/simulated'
                 '/synthetic_factorial_caches/union')
    have = {p.stem for p in cache.glob('*.parquet')}
    s = pd.read_parquet('output/samples_union.parquet').reset_index(drop=True)
    miss = s[~s.sample_id.isin(have)]
    print(','.join(map(str, sorted({i % 200 for i in miss.index}))))
    EOF

then submit just those:

    sbatch --array=<that list>%10 wrap_union_slow.sh

2026-08-26 run: 84 missing samples fell in only 16 slices
(2,3,5,7,10,11,13,19,20,73,74,84,100,101,193,194).
