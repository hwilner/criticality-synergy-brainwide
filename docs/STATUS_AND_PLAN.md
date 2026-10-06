# Status and Plan

## Status (as of repository creation)

| Component | State |
|---|---|
| `critsynergy.avalanches` (binning, avalanches, branching ratio, MLE exponent, benchmarks) | Done, tested |
| `critsynergy.synergy` (discrete/Gaussian entropy, O-information, TC, TE) | Done, tested |
| Unit tests (16, all passing offline) | Done |
| Docs: INTRODUCTION / EXTENDED_INTRODUCTION / METHODS | Done |
| IBL data pipeline | Not started (issue #1) |
| Empirical coupling analysis | Not started (issues #2–#6) |

## Plan

The work is decomposed into atomic issues on the GitHub issue tracker:

1. **IBL downloader** — ONE-API session fetch (spikes, clusters, trials);
   inclusion criteria; local cache.
2. **Avalanche and sigma pipeline** — binning sweep, avalanche
   extraction, branching-ratio estimation with bootstrap CI and
   subsampling-fraction stability.
3. **O-information pipeline** — discrete and Gaussian estimators with
   shuffle-bias correction; region-balanced unit groups.
4. **Coupling analysis** — GAMM of Omega vs distance-to-criticality
   with rate/bin/region controls; interior-optimum bootstrap test.
5. **Behaviour link** — accuracy/RT vs epoch-level sigma; mixed models.
6. **Surrogates and replication** — dithering/trial-shuffle/Galton-Watson
   battery; Allen Visual Coding end-to-end replication.
7. **Figures and write-up** — avalanche distributions, sigma atlas,
   Omega-vs-delta curve with CI, replication panel; results into
   `docs/CURRENT_RESULTS_AND_DISCUSSION.md`.

## Milestones

- **M1 (data-ready):** issues #1–#2 closed; `sessions.csv` + sigma table.
- **M2 (measures):** issue #3 closed; Omega/TE tables per session/epoch.
- **M3 (verdict):** issues #4–#6 closed; falsifier verdict stated.
- **M4 (paper-grade):** issue #7 closed.
