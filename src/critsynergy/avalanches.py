"""Avalanche extraction and critical-branching estimators.

A neural *avalanche* is a contiguous run of above-threshold population
activity; under the critical-branching hypothesis the size distribution
approaches a power law with exponent near -3/2 and the branching ratio
sigma approaches 1.  This module bins spike trains into population
activity, extracts avalanches, estimates the branching ratio, and fits
power-law exponents by maximum likelihood.
"""

from __future__ import annotations

import numpy as np


def bin_population(spike_trains: list[np.ndarray], t_min: float, t_max: float,
                   bin_size: float) -> np.ndarray:
    """Bin multiple spike trains into a population activity vector.

    Args:
        spike_trains: List of 1-D arrays of spike times (seconds).
        t_min: Start of the analysis window.
        t_max: End of the analysis window.
        bin_size: Bin width in the same units as spike times.

    Returns:
        Integer array of shape ``(n_bins,)`` with the total spike count
        per bin across all trains.
    """
    if bin_size <= 0:
        raise ValueError("bin_size must be positive")
    n_bins = int(np.floor((t_max - t_min) / bin_size))
    if n_bins < 1:
        raise ValueError("window shorter than one bin")
    counts = np.zeros(n_bins, dtype=int)
    for train in spike_trains:
        idx = np.floor((np.asarray(train) - t_min) / bin_size).astype(int)
        idx = idx[(idx >= 0) & (idx < n_bins)]
        np.add.at(counts, idx, 1)
    return counts


def extract_avalanches(activity: np.ndarray, threshold: float = 0.0) -> list[tuple[int, float]]:
    """Extract avalanches as runs of activity above a threshold.

    Args:
        activity: 1-D population activity time series.
        threshold: Bins with activity strictly greater than this belong
            to an avalanche.

    Returns:
        List of ``(start_bin, size)`` where ``size`` is the summed
        activity over the run.
    """
    activity = np.asarray(activity, dtype=float)
    out, start, size = [], None, 0.0
    for t, a in enumerate(activity):
        if a > threshold:
            if start is None:
                start, size = t, 0.0
            size += a
        elif start is not None:
            out.append((start, size))
            start = None
    if start is not None:
        out.append((start, size))
    return out


def branching_ratio(activity: np.ndarray) -> float:
    """Estimate the branching ratio sigma from population activity.

    The estimator is the least-squares slope of
    ``A_{t+1} = sigma * A_t + const`` over bins with ``A_t > 0``, the
    standard multistep-regression estimator for the expected number of
    descendants per ancestor in a branching process.

    Args:
        activity: 1-D population activity time series.

    Returns:
        Estimated sigma; values near 1 indicate near-critical branching.
    """
    a = np.asarray(activity, dtype=float)
    x, y = a[:-1], a[1:]
    mask = x > 0
    if mask.sum() < 3:
        raise ValueError("too few active bins to estimate sigma")
    x, y = x[mask], y[mask]
    xm, ym = x.mean(), y.mean()
    denom = np.sum((x - xm) ** 2)
    if denom == 0:
        raise ValueError("degenerate activity series")
    return float(np.sum((x - xm) * (y - ym)) / denom)


def avalanche_sizes(activity: np.ndarray, threshold: float = 0.0) -> np.ndarray:
    """Sizes of all avalanches in an activity series.

    Args:
        activity: 1-D population activity time series.
        threshold: Avalanche threshold, see :func:`extract_avalanches`.

    Returns:
        Array of avalanche sizes.
    """
    return np.array([s for _, s in extract_avalanches(activity, threshold)])


def powerlaw_alpha_mle(samples: np.ndarray, x_min: float | None = None) -> float:
    """Maximum-likelihood exponent of a continuous power law.

    For samples drawn from ``p(x) ~ x^{-alpha}`` with ``x >= x_min``, the
    ML estimator is ``alpha = 1 + n / sum(log(x_i / x_min))``.

    Args:
        samples: Positive samples (e.g. avalanche sizes).
        x_min: Lower cut-off; defaults to the sample minimum.

    Returns:
        Estimated exponent alpha (> 1 for normalisable data).
    """
    x = np.asarray(samples, dtype=float)
    x = x[x > 0]
    if x.size < 2:
        raise ValueError("need at least two positive samples")
    xm = float(x.min()) if x_min is None else float(x_min)
    x = x[x >= xm]
    return float(1.0 + x.size / np.sum(np.log(x / xm)))


def mean_avalanche_size_branching(sigma: float) -> float:
    """Mean total progeny of a subcritical branching process.

    For expected offspring sigma < 1 the mean total avalanche size
    (ancestor plus all descendants) is exactly ``1 / (1 - sigma)``.

    Args:
        sigma: Branching ratio.

    Returns:
        Mean avalanche size; ``inf`` for sigma >= 1.
    """
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if sigma >= 1:
        return float("inf")
    return 1.0 / (1.0 - sigma)


def simulate_branching_activity(sigma: float, n_ancestors: int,
                                seed: int = 0) -> np.ndarray:
    """Simulate total sizes of Poisson branching-process avalanches.

    Each individual leaves ``Poisson(sigma)`` offspring; the returned
    array holds the total size (ancestor included) of each of
    ``n_ancestors`` independent avalanches.  For ``sigma < 1`` the sample
    mean converges to ``1 / (1 - sigma)``.

    Args:
        sigma: Expected offspring per individual.
        n_ancestors: Number of avalanches to simulate.
        seed: Random seed.

    Returns:
        Array of total avalanche sizes (integers >= 1).
    """
    rng = np.random.default_rng(seed)
    sizes = np.empty(n_ancestors, dtype=int)
    for i in range(n_ancestors):
        total, current = 0, 1
        while current > 0:
            total += current
            current = int(rng.poisson(sigma, size=current).sum())
        sizes[i] = total
    return sizes
