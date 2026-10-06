# Materials and Methods

## 1. Data

**Primary.** [IBL Brain-Wide Map](https://www.internationalbrainlab.com/data)
(ONE API): Neuropixels spike times and cluster annotations, ~620k
neurons / 139 mice / 279 regions, plus aligned task events (stimulus,
choice, feedback). **Replication.** Allen Visual Coding Neuropixels
(AllenSDK). Sessions are included with ≥ 50 good units and ≥ 300
completed trials.

**Epochs.** Inter-trial baseline (−500 to 0 ms), stimulus epoch (0–250 ms
post-stimulus), decision epoch (250 ms pre-choice), feedback epoch
(0–500 ms post-feedback).

## 2. Population activity and avalanches

- Binning: population count vectors at Δ ∈ {1, 2, 4, 8} ms
  (`bin_population`), region-resolved (per structure) and pooled.
- Avalanches: runs above threshold θ = 0 (`extract_avalanches`); size =
  summed counts (`avalanche_sizes`).
- Exponents: power-law MLE `powerlaw_alpha_mle` with x_min chosen by the
  Clauset–Shalizi–Newman KS-minimisation sweep; comparison against
  log-normal and exponential alternatives by likelihood ratio.

## 3. Branching-ratio estimation

- Baseline: multistep regression estimator `branching_ratio`
  (slope of $A_{t+1}$ on $A_t$, active bins only), bootstrap CI (1000
  resamples over trials).
- Subsampling awareness: re-estimate at random unit fractions
  {100%, 50%, 25%, 12.5%} per session; retain sessions whose estimate is
  stable within CI overlap across fractions.
- Distance-to-criticality: $\delta = |1 - \hat\sigma|$ per session/epoch.

## 4. Information measures

- **O-information (discrete).** Spike counts per epoch discretised to
  {0, 1, ≥2} per unit; groups of 3–12 units sampled region-balanced;
  `o_information_discrete` (bias-corrected by subtracting the mean of 50
  trial-shuffled estimates).
- **O-information (Gaussian).** z-scored smoothed rates;
  `o_information_gaussian` on shrinkage-regularised covariance
  (Ledoit–Wolf).
- **Transfer entropy.** Stimulus/choice binary channels → population
  first principal component (discretised), `transfer_entropy_discrete`,
  lag 1 epoch-bin.
- **Total correlation.** `total_correlation_gaussian` as a scalable
  overall-dependence covariate.

## 5. Statistical model

Preregistered GAMM (Gaussian family):

$$
\Omega_{s,e} \sim f(\delta_{s,e}) + f(\bar r_{s,e}) +
f(\text{bin}) + \text{region} + (1 \mid s),
$$

with $f$ penalised splines, $\bar r$ mean rate. Interior-optimum test:
the location of the spline minimum within the observed δ-range with a
bootstrap 95% CI (resample sessions). Flat/monotone null accepted if the
CI excludes the interior or the curve's range < surrogate 99th
percentile. Behavioural link: trial accuracy/RT regressed on epoch δ
with the same controls (logistic/linear mixed models).

## 6. Surrogate battery

1. **Rate-preserving dithering** (uniform ±25 ms spike jitter): destroys
   correlations, keeps rates — must yield flat Ω–δ.
2. **Trial shuffling**: keeps within-trial statistics, breaks epoch
   alignment.
3. **Galton–Watson mock sessions** with known σ: pipeline must recover
   both σ and the analytic mean size $1/(1-\sigma)$ (unit-tested).
4. **Poisson null**: independent inhomogeneous Poisson units at matched
   rates — defines the no-coupling floor for Ω.

## 7. Validation battery (synthetic)

`tests/test_critsynergy.py` (16 tests, offline): binning correctness;
avalanche segmentation; σ recovery on Galton–Watson (±0.05); mean size
vs $1/(1-\sigma)$; supercritical divergence; power-law MLE recovery;
Ω anchors (independent 0, clones > 0, XOR = −1 bit); Gaussian entropy
closed form; Gaussian Ω of identity = 0; TE direction on a copy chain.

## 8. Software

Python 3.10+; NumPy/SciPy/scikit-learn/pandas. Optional: ONE-api (IBL),
allensdk (Allen), pyinform (TE cross-check). Google-style docstrings
throughout; tests run with the stdlib `unittest` runner, no data needed.

## 9. Limitations

- The regression σ-estimator is not the subsampling-scaling estimator of
  Levina & Priesemann; the fraction-stability check is a pragmatic
  substitute and is reported transparently.
- Discrete Ω estimation is sample-hungry; group sizes are capped where
  bias correction exceeds 20% of the raw value (reported per session).
- δ is measured per epoch, but epoch-level σ is noisy; δ enters the GAMM
  with measurement-error weights (inverse bootstrap variance).
- fMRI-scale claims are out of scope: this is a spiking-level project.
