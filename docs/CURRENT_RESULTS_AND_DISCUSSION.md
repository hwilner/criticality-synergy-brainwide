# Current Results and Discussion

## What works now (synthetic validation)

All 16 unit tests pass offline
(`PYTHONPATH=src python -m unittest discover -s tests`):

1. **Critical-branching machinery is exact where theory is exact.**
   On simulated Galton–Watson activity with $\sigma = 0.9$, the
   regression estimator recovers $\hat\sigma$ within 0.05 and the mean
   simulated avalanche size matches the closed form $1/(1-\sigma) = 10$
   within sampling error; the supercritical case correctly diverges.

2. **The power-law exponent estimator is calibrated.** On exact
   power-law samples with $\alpha = 2.5$, the MLE recovers the exponent
   within 0.15 (6000 samples) — the precision needed to distinguish the
   critical $-3/2$ avalanche law from alternatives.

3. **O-information hits its three analytic anchors.** Independent
   variables give $\Omega \approx 0$; identical copies give
   $\Omega > 0$ (redundancy); the XOR system gives $\Omega = -1$ bit
   exactly (synergy). The Gaussian variant reproduces the closed-form
   entropy $\tfrac{d}{2}\ln(2\pi e)$ for identity covariance and zero
   for independent variables.

4. **Transfer entropy is directional.** On a copy-chain ($y_t =
   x_{t-1}$), TE(x→y) > 0.5 bits while TE(y→x) ≈ 0 — the estimator sees
   the drive, not the shadow.

## Interpretation

Every estimator in the empirical claim is pinned to a case where the
truth is known analytically or by construction. The pipeline therefore
enters the data phase with no free tuning knobs: thresholds, group
sizes, and surrogate ensembles are fixed by METHODS, and the
preregistered GAMM either finds an interior synergy optimum at minimal
distance-to-criticality or it does not.

## Open questions

- **Which sigma estimator to trust at scale.** The regression estimator
  is fast and testable; the subsampling-scaling estimator (Levina &
  Priesemann) is more principled but heavier. Issue #2 compares them on
  10 pilot sessions before freezing.
- **Discrete vs Gaussian Ω.** Both are computed everywhere; concordance
  across estimators is itself a robustness result (reported in issue #3).
- **Region heterogeneity.** 279 regions cannot each be estimated well;
  the primary analysis pools into ~12 macro-regions, with per-region
  results as exploratory.

## Risks

| Risk | Mitigation |
|---|---|
| Subsampling mimics the coupling | Fraction-stability requirement + dithering surrogates (METHODS §6) |
| Rate–synergy confound | Rate spline in the GAMM + rate-matched nulls |
| Epoch-level sigma too noisy | Inverse-variance weighting; session-level sensitivity analysis |
| Multiplicity across regions/epochs | Preregistered primary endpoint (pooled Ω vs δ); BH FDR elsewhere |
