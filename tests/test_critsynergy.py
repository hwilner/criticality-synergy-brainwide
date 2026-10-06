"""Unit tests for the critsynergy package."""

import unittest

import numpy as np

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


class TestBinningAndAvalanches(unittest.TestCase):
    """Spike binning and avalanche extraction."""

    def test_bin_population_counts(self):
        trains = [np.array([0.05, 0.06, 0.5]), np.array([0.55])]
        counts = bin_population(trains, 0.0, 1.0, 0.1)
        self.assertEqual(counts.sum(), 4)
        self.assertEqual(counts[0], 2)
        self.assertEqual(counts[5], 2)

    def test_extract_avalanches_runs(self):
        activity = np.array([0, 2, 3, 0, 1, 1, 0, 5])
        avs = extract_avalanches(activity)
        self.assertEqual([s for _, s in avs], [5.0, 2.0, 5.0])
        self.assertEqual(avs[0][0], 1)

    def test_empty_activity_no_avalanches(self):
        self.assertEqual(avalanche_sizes(np.zeros(10)).size, 0)

    def test_bad_bin_size_raises(self):
        with self.assertRaises(ValueError):
            bin_population([], 0.0, 1.0, 0.0)


class TestBranchingProcess(unittest.TestCase):
    """Critical-branching estimators against exact theory."""

    def test_branching_ratio_recovers_sigma(self):
        rng = np.random.default_rng(0)
        sigma_true = 0.9
        # Galton-Watson activity series: A_{t+1} = sum Poisson(sigma) over A_t
        a = np.zeros(4000)
        a[0] = 50
        for t in range(len(a) - 1):
            if a[t] > 0:
                a[t + 1] = rng.poisson(sigma_true, size=int(a[t])).sum()
            else:
                a[t + 1] = 5  # small immigration keeps the process alive
        est = branching_ratio(a)
        self.assertAlmostEqual(est, sigma_true, delta=0.05)

    def test_mean_size_matches_exact(self):
        sizes = simulate_branching_activity(0.9, n_ancestors=4000, seed=1)
        self.assertAlmostEqual(sizes.mean(), 10.0, delta=1.5)
        self.assertAlmostEqual(mean_avalanche_size_branching(0.9), 10.0)

    def test_supercritical_mean_is_inf(self):
        self.assertEqual(mean_avalanche_size_branching(1.0), float("inf"))

    def test_powerlaw_mle_recovers_alpha(self):
        rng = np.random.default_rng(2)
        alpha_true = 2.5
        u = rng.uniform(1e-9, 1.0, size=6000)
        x = u ** (-1.0 / (alpha_true - 1.0))  # exact power-law samples
        est = powerlaw_alpha_mle(x, x_min=1.0)
        self.assertAlmostEqual(est, alpha_true, delta=0.15)


class TestSynergyMeasures(unittest.TestCase):
    """O-information and friends on analytically known systems."""

    def test_independent_variables_zero_oinfo(self):
        rng = np.random.default_rng(3)
        data = rng.integers(0, 2, size=(20000, 4))
        self.assertAlmostEqual(o_information_discrete(data), 0.0, delta=0.05)

    def test_identical_copies_redundant(self):
        rng = np.random.default_rng(4)
        x = rng.integers(0, 4, size=20000)
        data = np.stack([x, x, x, x], axis=1)
        self.assertGreater(o_information_discrete(data), 1.0)

    def test_xor_is_synergistic(self):
        rng = np.random.default_rng(5)
        x = rng.integers(0, 2, size=50000)
        y = rng.integers(0, 2, size=50000)
        z = (x + y) % 2
        data = np.stack([x, y, z], axis=1)
        # exact value: (3-2)*H(x,y,z) + 3*(1 - 2) = 2 - 3 = -1 bit
        self.assertAlmostEqual(o_information_discrete(data), -1.0, delta=0.05)

    def test_gaussian_entropy_identity(self):
        d = 5
        expected = 0.5 * d * np.log(2 * np.pi * np.e)
        self.assertAlmostEqual(gaussian_entropy(np.eye(d)), expected, places=10)

    def test_gaussian_oinfo_independent_is_zero(self):
        self.assertAlmostEqual(o_information_gaussian(np.eye(4)), 0.0,
                               places=10)
        self.assertAlmostEqual(total_correlation_gaussian(np.eye(3)), 0.0,
                               places=10)

    def test_gaussian_total_correlation_positive_for_copies(self):
        cov = np.full((3, 3), 0.9) + np.eye(3) * 0.1
        self.assertGreater(total_correlation_gaussian(cov), 0.0)

    def test_non_psd_covariance_raises(self):
        with self.assertRaises(ValueError):
            gaussian_entropy(np.array([[1.0, 2.0], [2.0, 1.0]]))

    def test_transfer_entropy_direction(self):
        rng = np.random.default_rng(6)
        n = 40000
        x = rng.integers(0, 2, size=n)
        y = np.zeros(n, dtype=int)
        y[1:] = x[:-1]  # y copies previous x
        te_xy = transfer_entropy_discrete(x, y, lag=1)
        te_yx = transfer_entropy_discrete(y, x, lag=1)
        self.assertGreater(te_xy, 0.5)
        self.assertAlmostEqual(te_yx, 0.0, delta=0.05)


if __name__ == "__main__":
    unittest.main()
