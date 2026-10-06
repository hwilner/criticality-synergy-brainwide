# Introduction

## Criticality–Synergy Coupling in Brain-Wide Spiking Activity

> **Core hypothesis.** Epochs in which the population branching ratio is
> closest to its critical value exhibit maximal synergistic (multivariate)
> information about stimuli and choices; redundancy dominates away from
> criticality, and behavioural performance follows an inverted-U in
> distance-from-criticality.
>
> **Falsifier.** Synergy shows no interior optimum as a function of
> distance-to-criticality (flat or monotonic relationship) once
> subsampling, firing-rate, and bin-size confounds are controlled with
> surrogate ensembles.

---

## 1. Motivation

The **neuronal criticality hypothesis** holds that cortical networks
self-organise near a continuous dynamical transition — a branching
process with mean offspring (branching ratio) $\sigma = 1$ — because
model studies show that dynamic range, information capacity, and
transmission are optimised there (Beggs & Plenz, 2003; Kinouchi &
Copelli, 2006). Two decades of avalanche analyses report approximately
power-law size distributions, though inference is technically delicate:
subsampling alone can manufacture or destroy apparent criticality
(Priesemann et al., 2009; Levina & Priesemann, 2017).

Independently, **multivariate information theory** has reshaped how we
think about population codes. Information in a neural population
decomposes into *redundant* components (many neurons saying the same
thing), *unique* components, and *synergistic* components (information
available only from joint observation — the XOR pattern being the
canonical example; Williams & Beer, 2010). The O-information (Rosas et
al., 2019) gives a scalable signed summary: negative means
synergy-dominated, positive means redundancy-dominated.

These two research programmes have largely passed each other by.
Criticality work rarely asks *what kind* of information the network
carries; information-theoretic work rarely asks *which dynamical regime*
produces its decompositions. Yet the theoretical link is natural: near
critical branching, activity is maximally collective — long-range
correlations, large fluctuations, modes with slow decay — which is
precisely the regime where joint (synergistic) codes should be
cheapest to maintain. The IBL Brain-Wide Map (IBL et al., 2023) makes
the joint question testable for the first time: the same task, the same
preprocessing, hundreds of thousands of neurons across 279 regions.

## 2. Named concepts

- **Neural avalanche.** A run of consecutive time bins with
  above-threshold population activity; its *size* is the total activity
  in the run. At criticality, sizes distribute as $P(s) \sim s^{-3/2}$
  (the mean-field branching-process exponent).
- **Branching ratio $\sigma$.** The expected number of descendant spikes
  per ancestor spike. $\sigma < 1$: activity dies out (subcritical);
  $\sigma > 1$: runaway (supercritical); $\sigma = 1$: critical.
- **Subsampling problem.** Recordings see a tiny fraction of all neurons;
  naive avalanche statistics of a subsampled critical system look
  subcritical, and vice versa (Levina & Priesemann, 2017). Estimators
  must be subsampling-aware.
- **Partial information decomposition (PID).** Williams & Beer's (2010)
  framework splitting the information that sources provide about a
  target into redundant, unique, and synergistic atoms.
- **O-information.** A computationally cheap signed measure of whether a
  system's high-order dependencies are synergy- or redundancy-dominated
  (Rosas et al., 2019); scales to tens of variables where full PID
  cannot.
- **Transfer entropy (TE).** Directed, history-conditioned information
  flow $X \to Y$ (Schreiber, 2000): the information $X$'s past adds
  about $Y$'s future beyond $Y$'s own past.
- **IBL Brain-Wide Map.** The International Brain Laboratory's
  standardised Neuropixels dataset: ~620k neurons, 139 mice, 279 brain
  regions, one decision task (IBL et al., 2023).

## 3. Mathematical bridge

### 3.1 Branching processes: the model under the metaphor

A Galton–Watson branching process gives each active neuron $i$ at time
$t$ a Poisson number of "descendant" activations at $t+1$ with mean
$\sigma$. *Step-by-step:* let $Z_t$ be the activity count and
$Z_{t+1} = \sum_{i=1}^{Z_t} \xi_i$ with $\mathbb E \xi_i = \sigma$. Then
$\mathbb E[Z_{t+1} \mid Z_t] = \sigma Z_t$ — so the regression slope of
$A_{t+1}$ on $A_t$ estimates $\sigma$ directly (our `branching_ratio`).
Total avalanche size $S = \sum_t Z_t$ satisfies
$\mathbb E[S] = 1 + \sigma + \sigma^2 + \cdots = 1/(1-\sigma)$ for
$\sigma < 1$ — diverging as $\sigma \to 1^-$: the mean avalanche blows
up at criticality (our tests verify both relations against simulation).

The size *distribution* of a critical Galton–Watson avalanche with
Poisson offspring is exactly the Borel distribution,
$P(S{=}s) = (s^{s-1} e^{-s}) / s! \sim s^{-3/2}$ for large $s$ — the
origin of the famous $-3/2$ exponent. We fit exponents by the
maximum-likelihood estimator
$\hat\alpha = 1 + n / \sum_i \ln(s_i / s_{\min})$ and benchmark it on
exact power-law samples.

### 3.2 O-information: synergy minus redundancy in one signed number

For $n$ variables the O-information is

$$
\Omega_n = (n-2)\, H(X_1,\dots,X_n) + \sum_{i=1}^n
\big[ H(X_i) - H(X_{-i}) \big],
$$

with $H$ the Shannon entropy and $X_{-i}$ all variables except $X_i$.
*Why this combination, step by step.* Total correlation
$\mathrm{TC} = \sum_i H(X_i) - H(X)$ measures *all* dependence.
Dual total correlation $\mathrm{DTC} = H(X) - \sum_i
[H(X) - H(X_{-i})]$ measures dependence *that survives removing any
single variable* — the redundant part. Then
$\Omega = \mathrm{TC} - \mathrm{DTC} = (n-2)H(X) + \sum_i[H(X_i) -
H(X_{-i})]$: synergy minus redundancy. Sanity anchors (all in our test
suite):

- **Independent** variables: every $H$ term decomposes, $\Omega = 0$.
- **Identical copies** $X_i = X$: $H(X) = H(X_i) = H(X_{-i}) = h$, so
  $\Omega = (n-2)h > 0$ — redundancy-dominated, as expected.
- **XOR** $Z = X \oplus Y$ with independent bits $X, Y$:
  $H(X,Y,Z) = 2$ bits, each $H(\cdot) - H(\cdot,\cdot) = 1 - 2 = -1$,
  so $\Omega = (3-2)\cdot 2 + 3(-1) = -1$ bit — synergy-dominated,
  exactly the canonical synergistic system.

For Gaussian variables all entropies are closed form,
$H = \tfrac12 \ln\big[(2\pi e)^d \det\Sigma\big]$, giving a fast
continuous estimator (`o_information_gaussian`); the discrete
histogram estimator (`o_information_discrete`) handles spike counts.

### 3.3 The coupling statistic

Per session $s$ and epoch $e$ we compute the distance to criticality
$\delta_{s,e} = |1 - \hat\sigma_{s,e}|$ (with the subsampling-aware
confidence interval of the estimator) and the population O-information
$\Omega_{s,e}$ (over region-balanced neuron groups), plus TE between
stimulus/choice channels and population latents. The preregistered model
is the generalised additive mixed model

$$
\Omega_{s,e} \sim f(\delta_{s,e}) + (\text{rate, region, bin size})
+ (1 \mid \text{session}),
$$

and the hypothesis is the existence of an interior minimum of $f$
(minimum of redundancy / maximum of synergy at $\delta \to 0$).
The falsifier is $f$ flat or monotonic across the observed range.

## 4. Why the confounds are the heart of the project

1. **Subsampling** biases both axes: fewer neurons → lower apparent
   $\sigma$ *and* altered $\Omega$. Control: subsample the same session
   at multiple fractions; require the $\Omega$–$\delta$ relation to be
   invariant.
2. **Firing rate** drives both apparent criticality and entropy
   estimates. Control: rate-matched surrogate spike trains
   (spike-time dithering preserving rates, destroying correlations)
   must show no relation.
3. **Bin size** changes avalanche segmentation. Control: full analysis
   at 1, 2, 4, 8 ms bins; report the relation's stability.
4. **Estimator bias** at finite samples. Control: analytic anchors
   (this repo's tests) plus simulation-based calibration on
   Galton–Watson data with matched rates and counts.

## 5. Empirical pipeline

1. **Data.** IBL Brain-Wide Map spike times + behaviour (primary);
   Allen Visual Coding Neuropixels (replication).
2. **Per session:** bin population activity (region-resolved), extract
   avalanches, estimate $\hat\sigma$ with bootstrap CI, fit avalanche
   exponents.
3. **Per epoch:** compute $\Omega$ (discrete on count vectors; Gaussian
   on z-scored rates) and TE from stimulus/choice to population.
4. **Across sessions:** fit the mixed model above; test the inverted-U
   against flat/monotonic alternatives by cross-validated likelihood.
5. **Replication:** repeat end-to-end in Allen data before any claim.

## 6. Relation to the programme

This repo is the *criticality/information* member of the physics-brain
programme (thermodynamic EPR, Kibble–Zurek quenches, Berry holonomy,
active nematics, and quantum-inspired geometry repos). It shares the
house rule: every hypothesis ships with a falsifier and a surrogate
battery, and all estimators are pinned to analytic ground truth before
touching data.

## References

1. Beggs, J. M. & Plenz, D. (2003). Neuronal avalanches in neocortical
   circuits. *J. Neurosci.* 23, 11167–11177.
2. Kinouchi, O. & Copelli, M. (2006). Optimal dynamical range of
   excitable networks at criticality. *Nat. Phys.* 2, 348–351.
3. Priesemann, V. et al. (2009). Neuronal avalanches differ from
   wakefulness to deep sleep. *PLoS ONE* 4, e8985.
4. Levina, A. & Priesemann, V. (2017). Subsampling scaling.
   *Nat. Commun.* 8, 15140.
5. Williams, P. L. & Beer, R. D. (2010). Nonnegative decomposition of
   multivariate information. arXiv:1004.2515.
6. Rosas, F. E., Mediano, P. A. M., Gastpar, M. & Jensen, H. J. (2019).
   Quantifying high-order interdependencies via multivariate extensions
   of the mutual information. *Phys. Rev. E* 100, 032305.
7. Schreiber, T. (2000). Measuring information transfer.
   *Phys. Rev. Lett.* 85, 461.
8. International Brain Laboratory et al. (2023). A brain-wide map of
   neural activity during complex behaviour. bioRxiv
   2023.07.04.547681; Nature (2025).
9. de Vries, S. E. J. et al. (2020). A large-scale standardized
   physiological survey. *Nat. Neurosci.* 23, 138–151.
10. Harris, T. E. (1963). *The Theory of Branching Processes*. Springer.
