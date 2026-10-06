"""Criticality-synergy analysis of brain-wide spiking activity.

Public API: avalanche extraction and critical-branching estimators
(:mod:`critsynergy.avalanches`) and multivariate information measures
(:mod:`critsynergy.synergy`).
"""

from critsynergy.avalanches import (
    avalanche_sizes,
    bin_population,
    branching_ratio,
    extract_avalanches,
    mean_avalanche_size_branching,
    powerlaw_alpha_mle,
    simulate_branching_activity,
)
from critsynergy.synergy import (
    entropy_discrete,
    gaussian_entropy,
    joint_entropy_discrete,
    o_information_discrete,
    o_information_gaussian,
    total_correlation_gaussian,
    transfer_entropy_discrete,
)

__all__ = [
    "avalanche_sizes",
    "bin_population",
    "branching_ratio",
    "entropy_discrete",
    "extract_avalanches",
    "gaussian_entropy",
    "joint_entropy_discrete",
    "mean_avalanche_size_branching",
    "o_information_discrete",
    "o_information_gaussian",
    "powerlaw_alpha_mle",
    "simulate_branching_activity",
    "total_correlation_gaussian",
    "transfer_entropy_discrete",
]

__version__ = "0.1.0"
