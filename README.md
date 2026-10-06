# Criticality–Synergy Coupling in Brain-Wide Spiking

**Does multivariate information synergy peak when a neural population
sits closest to critical branching?**

Two of the strongest ideas in systems neuroscience have never been tested
against each other at scale. The **criticality hypothesis** says cortical
networks operate near a branching-process critical point (avalanche size
distributions ~ power laws, branching ratio σ ≈ 1), where dynamic range
and information transmission are theoretically maximised. **Multivariate
information theory** says what matters for computation is not raw
activity but *synergy* — information that only emerges from combinations
of neurons (PID, O-information). This project tests, on the largest
standardised spike dataset in the world, whether the two peaks coincide:
is synergy maximised exactly at the edge of critical branching?

- **Core hypothesis.** Behavioural epochs whose estimated branching ratio
  is closest to the data-derived critical point show the highest
  stimulus/choice-related synergistic information; redundancy dominates
  away from criticality; behavioural performance follows an inverted-U in
  distance-from-criticality.
- **Falsifier.** Synergy (O-information and PID-style measures) is flat
  or monotonic in distance-to-criticality — no interior optimum exists —
  once subsampling, firing-rate, and bin-size confounds are controlled.

## Repository layout

```
src/critsynergy/
    avalanches.py  Population binning, avalanche extraction, branching-
                   ratio regression estimator, power-law MLE exponent,
                   exact branching-process benchmarks.
    synergy.py     Discrete & Gaussian Shannon entropies, O-information,
                   total correlation, discrete transfer entropy.
tests/
    test_critsynergy.py  16 unit tests against analytic ground truth
                         (Galton-Watson recovery, XOR = -1 bit, etc.).
docs/
    INTRODUCTION.md              Publication-quality motivation and theory.
    EXTENDED_INTRODUCTION.md     High-school-level walkthrough with links.
    METHODS.md                   Materials and methods (data, pipeline).
    STATUS_AND_PLAN.md           Task board snapshot.
    CURRENT_RESULTS_AND_DISCUSSION.md
```

## Quickstart

```bash
pip install -r requirements.txt
export PYTHONPATH=src        # or: pip install -e .
python -m unittest discover -s tests -v
```

```python
from critsynergy import (simulate_branching_activity, branching_ratio,
                         mean_avalanche_size_branching, o_information_discrete)

sizes = simulate_branching_activity(sigma=0.9, n_ancestors=4000, seed=1)
print(sizes.mean(), mean_avalanche_size_branching(0.9))   # ~10, exactly 10
```

## Empirical target

The [IBL Brain-Wide Map](https://www.internationalbrainlab.com/data):
~620k neurons, 139 mice, 279 brain regions, one standardised decision
task — with [Allen Visual Coding](https://brain-map.org/our-research/circuits-behavior/visual-coding)
as an independent replication. Per session and epoch we estimate the
branching ratio (subsampling-aware), compute O-information over
region-resolved populations, and regress synergy on
distance-from-criticality with behavioural covariates (see
`docs/METHODS.md`).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). All code uses Google-style
docstrings and ships with unit tests.

## License

MIT (see LICENSE).
