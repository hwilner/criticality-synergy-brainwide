"""Information-theoretic synergy/redundancy measures.

The central statistic is the *O-information*

    Omega_n = (n - 2) H(X_1..n) + sum_i [ H(X_i) - H(X_{-i}) ],

which is negative for synergy-dominated systems and positive for
redundancy-dominated ones.  Both discrete (histogram) and Gaussian
(closed-form) estimators are provided, plus discrete transfer entropy

    TE(X -> Y) = H(Y+ | Y) - H(Y+ | Y, X)

as a directed-dependence complement.  All estimators operate on samples,
not density models, so the same code runs on spike counts, latents, or
behavioural variables.
"""

from __future__ import annotations

import numpy as np


def entropy_discrete(labels: np.ndarray) -> float:
    """Shannon entropy (bits) of a discrete sample.

    Args:
        labels: 1-D array of hashable labels (e.g. integer bins).

    Returns:
        ``-sum p log2 p`` over the empirical distribution.
    """
    _, counts = np.unique(np.asarray(labels), return_counts=True)
    p = counts / counts.sum()
    return float(-np.sum(p * np.log2(p)))


def joint_entropy_discrete(data: np.ndarray) -> float:
    """Joint Shannon entropy (bits) of columns of a sample matrix.

    Args:
        data: Array of shape ``(n_samples, n_vars)``; each row is one
            joint observation, treated as a tuple label.

    Returns:
        Joint entropy in bits.
    """
    from collections import Counter

    data = np.asarray(data)
    if data.ndim != 2:
        raise ValueError("data must be 2-D (n_samples, n_vars)")
    counts = np.array(list(Counter(map(tuple, data)).values()), dtype=float)
    p = counts / counts.sum()
    return float(-np.sum(p * np.log2(p)))


def o_information_discrete(data: np.ndarray) -> float:
    """O-information (bits) of a discrete multivariate sample.

    Negative values indicate synergy-dominated dependence (e.g. XOR),
    positive values redundancy-dominated dependence (e.g. copies).

    Args:
        data: Array of shape ``(n_samples, n_vars)`` with ``n_vars >= 3``.

    Returns:
        O-information in bits.
    """
    data = np.asarray(data)
    n = data.shape[1]
    if n < 3:
        raise ValueError("O-information requires at least 3 variables")
    h_all = joint_entropy_discrete(data)
    acc = 0.0
    for i in range(n):
        h_i = entropy_discrete(data[:, i])
        h_rest = joint_entropy_discrete(np.delete(data, i, axis=1))
        acc += h_i - h_rest
    return float((n - 2) * h_all + acc)


def gaussian_entropy(cov: np.ndarray) -> float:
    """Differential entropy (nats) of a Gaussian with covariance ``cov``.

    ``H = 0.5 * log((2 pi e)^d det(cov))``; computed from the eigenvalues
    for numerical stability.

    Args:
        cov: Symmetric positive-definite covariance, shape ``(d, d)``.

    Returns:
        Entropy in nats.
    """
    cov = np.asarray(cov, dtype=float)
    eigs = np.linalg.eigvalsh(cov)
    if np.any(eigs <= 0):
        raise ValueError("covariance must be positive definite")
    d = cov.shape[0]
    return float(0.5 * (d * np.log(2 * np.pi * np.e) + np.sum(np.log(eigs))))


def total_correlation_gaussian(cov: np.ndarray) -> float:
    """Total correlation (multi-information, nats) of a Gaussian.

    ``TC = sum_i H(X_i) - H(X)``, zero iff all variables are independent.

    Args:
        cov: Symmetric positive-definite covariance, shape ``(d, d)``.

    Returns:
        Total correlation in nats.
    """
    cov = np.asarray(cov, dtype=float)
    d = cov.shape[0]
    h_joint = gaussian_entropy(cov)
    h_marginals = sum(gaussian_entropy(cov[[i], :][:, [i]]) for i in range(d))
    return float(h_marginals - h_joint)


def o_information_gaussian(cov: np.ndarray) -> float:
    """O-information (nats) of a Gaussian with covariance ``cov``.

    Args:
        cov: Symmetric positive-definite covariance, shape ``(d, d)``
            with ``d >= 3``.

    Returns:
        O-information in nats; sign convention as in the discrete case.
    """
    cov = np.asarray(cov, dtype=float)
    n = cov.shape[0]
    if n < 3:
        raise ValueError("O-information requires at least 3 variables")
    h_all = gaussian_entropy(cov)
    acc = 0.0
    for i in range(n):
        h_i = gaussian_entropy(cov[[i], :][:, [i]])
        rest = np.delete(np.delete(cov, i, axis=0), i, axis=1)
        acc += h_i - gaussian_entropy(rest)
    return float((n - 2) * h_all + acc)


def transfer_entropy_discrete(source: np.ndarray, target: np.ndarray,
                              lag: int = 1) -> float:
    """Discrete transfer entropy (bits) from source to target.

    ``TE = H(target_future | target_past) - H(target_future | target_past,
    source_past)`` with single-sample histories.

    Args:
        source: 1-D array of discrete source samples.
        target: 1-D array of discrete target samples, same length.
        lag: History lag in samples.

    Returns:
        Transfer entropy in bits (>= 0 up to finite-sample bias).
    """
    s = np.asarray(source)
    t = np.asarray(target)
    if s.shape != t.shape:
        raise ValueError("source and target must have equal length")
    if len(s) <= lag:
        raise ValueError("series shorter than lag")
    y_future = t[lag:]
    y_past = t[:-lag]
    x_past = s[:-lag]
    h_yy = (joint_entropy_discrete(np.stack([y_future, y_past], axis=1))
            - entropy_discrete(y_past))
    joint = np.stack([y_future, y_past, x_past], axis=1)
    h_yyx = (joint_entropy_discrete(joint)
             - joint_entropy_discrete(np.stack([y_past, x_past], axis=1)))
    return float(h_yy - h_yyx)
