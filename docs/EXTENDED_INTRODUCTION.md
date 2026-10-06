# Extended Introduction (start here if this is all new)

This page explains the project assuming **no college maths or physics**.
Every technical term gets a plain definition, an analogy, and a link. The
math is built step by step with the *why* stated at each step.

---

## 1. The one-sentence version

> Scientists have two separate guesses about why brains compute well —
> "they sit at a special balance point between fizzing out and exploding"
> and "they combine information across neurons in clever joint codes".
> We test whether those two guesses are actually **the same fact**, using
> the biggest standard brain-recording dataset in the world.

## 2. Two analogies that carry the whole project

### 2.1 The forest-fire balance (criticality)

Imagine lighting fires in a forest. Each fire can spark new fires in
neighbouring trees. The key number is the **branching ratio** σ: on
average, how many new fires does one fire start?

* σ < 1: fires fizzle out. Boring, safe — but nothing spreads.
* σ > 1: the whole forest burns down. Activity explodes — useless.
* σ = 1: the knife edge. Fires of *every size* happen — from one tree to
  huge patches — and the forest is maximally "alive": a spark can travel
  any distance, but the system never fully explodes.

Physicists call the knife edge a **critical point**. In 2003, John Beggs
and Dietmar Plenz measured bursts of activity in slices of cortex and
found exactly this signature: bursts (**avalanches**) of every size,
distributed as a power law, like the forest at σ = 1. The **criticality
hypothesis** says real brains park themselves there because models show
that at σ = 1 the network has the largest possible dynamic range and
carries the most information.
([Beggs & Plenz 2003](https://doi.org/10.1523/JNEUROSCI.23-35-11167.2003);
nice explainer: [Quanta, "On the Edge of
Chaos"](https://www.quantamagazine.org/))

### 2.2 The three-friends riddle (synergy)

Three friends each whisper a clue to you. Sometimes any one clue alone
solves the riddle — the clues are **redundant**. Sometimes clue A alone
says nothing and clue B alone says nothing, but **A and B together** solve
it completely — that's **synergy**: information that exists only in the
*combination*.

The textbook example: the XOR rule. Two coins, X and Y; the answer is
"are they different?" (XOR). Looking at either coin alone tells you zero
about the answer. Looking at both tells you everything. The information
is 100% synergistic.

Brains might use both styles. A **redundant** code is robust (lose a
neuron, keep the message). A **synergistic** code is compact and
powerful (joint patterns carry what singles cannot). Information theory
can measure both — the number we use most is the **O-information**: a
single signed score where *negative = synergy-dominated*,
*positive = redundancy-dominated*, zero = independent.
([Rosas et al. 2019](https://doi.org/10.1103/PhysRevE.100.032305))

### 2.3 The project in one line

Is the forest-fire knife edge (σ = 1) **also** where the three-friends
riddle is solved best — does **synergy peak at criticality**? Two famous
ideas would turn out to be one.

## 3. The brain data

The [International Brain Laboratory](https://www.internationalbrainlab.com/)
(IBL) recorded **over 600,000 neurons** from 139 mice across 279 brain
regions, all doing the same decision task (spot a visual stimulus, turn a
wheel). Everything is standardised and public — the perfect arena: we can
measure the fire-spreading number σ and the synergy score Ω in the *same
neurons, same epochs*, and ask whether the best-synergy epochs are
exactly the near-critical ones. We then double-check any finding in the
independent [Allen Brain
Observatory](https://brain-map.org/our-research/circuits-behavior/visual-coding)
dataset.

## 4. The math, step by step

### Step 1 — counting avalanches

Bin time into tiny windows (say 4 ms) and count spikes across the
population. Any run of non-empty bins is an **avalanche**; its size is
the total spikes in the run. (`extract_avalanches` — ten lines of code.)

### Step 2 — measuring the fire number σ

If every spike at time $t$ spawns on average σ spikes at $t+1$, then
expected activity tomorrow is σ times activity today:
$\mathbb E[A_{t+1}] = \sigma A_t$. So **fit a line** of tomorrow vs
today; its slope *is* σ. That is the whole estimator
(`branching_ratio`). Our tests simulate fake brains with a known σ and
check we get it back.

A bonus exact fact: the average avalanche size in such a process is
$\langle s \rangle = 1 + \sigma + \sigma^2 + \dots = 1/(1-\sigma)$.
At σ = 0.9 that's 10; as σ → 1 it explodes to infinity — mathematically,
that's *why* criticality feels special. (The test suite checks simulated
avalanches average 10 at σ = 0.9.)

### Step 3 — entropy, the ruler of information

**Entropy** $H$ counts how many yes/no questions you need to pin down a
signal: a fair coin has $H = 1$ bit, a double-headed coin has $H = 0$.
Joint entropy $H(X, Y)$ is the same for two signals together. Everything
else is addition and subtraction of entropies.

### Step 4 — the O-information formula, demystified

$$
\Omega_n = (n-2)\,H(\text{all}) + \sum_i \big[H(X_i) - H(\text{all but } i)\big]
$$

Don't memorise it — check the three anchor cases (our unit tests do
exactly this):

* **All independent:** every bracket and the first term cancel → Ω = 0.
* **All copies of each other** (4 clones): Ω = (4−2)·h > 0 → redundancy
  wins, as it should.
* **XOR** (the riddle): $H(\text{all}) = 2$ bits, each bracket = 1 − 2 =
  −1, so Ω = 2 − 3 = **−1 bit** → pure synergy detected.

If a formula gives the right answer on all three, it earns trust.
([O-information paper](https://doi.org/10.1103/PhysRevE.100.032305);
[Williams & Beer PID](https://arxiv.org/abs/1004.2515))

### Step 5 — the actual test

For each recording epoch: compute σ̂ (how close to 1?) and Ω (how
synergistic?). Then across hundreds of epochs, plot Ω against
"distance from criticality" $|1 - \hat\sigma|$ and ask:

* **Hypothesis:** the best synergy sits at the smallest distance — an
  interior optimum (inverted-U).
* **Falsifier:** the curve is flat or monotone — criticality has nothing
  special to do with synergy.

## 5. Why this is hard (and how we cheat-proof it)

1. **We see only some neurons.** Subsampling fakes both criticality and
   synergy changes. → We re-run everything at several subsampling
   fractions; the conclusion must not move.
2. **Firing rate changes look like everything.** → We make fake spike
   trains with identical rates but shuffled timing; they must show no
   effect.
3. **Bin size is arbitrary.** → Run at 1/2/4/8 ms; the effect must
   survive all.
4. **Estimators lie on small data.** → Every estimator is unit-tested
   against fake data where the truth is known (16 tests, all offline).

## 6. Keyword table

| Term | Plain meaning | Where used | Learn more |
|---|---|---|---|
| Avalanche | A burst of population activity from ignition to silence | Step 1 | [Beggs & Plenz 2003](https://doi.org/10.1523/JNEUROSCI.23-35-11167.2003) |
| Branching ratio σ | Average "children" per spike; 1 = critical | Step 2 | [Criticality review](https://doi.org/10.3389/fphys.2010.00024) |
| Critical point | The knife-edge regime between dying and exploding | §2.1 | [Kinouchi & Copelli 2006](https://doi.org/10.1038/nphys289) |
| Power law | "No typical size": small common, huge rare | Step 1 | [Power law primer](https://en.wikipedia.org/wiki/Power_law) |
| Subsampling | Seeing a fraction of neurons; distorts stats | §5 | [Levina & Priesemann 2017](https://doi.org/10.1038/ncomms15140) |
| Entropy | Expected surprise, in bits | Step 3 | [Khan Academy: information entropy](https://www.khanacademy.org/computing/computer-science/informationtheory) |
| Synergy | Info only in the *combination* of signals | §2.2 | [Williams & Beer](https://arxiv.org/abs/1004.2515) |
| Redundancy | Same info carried by many signals | §2.2 | same as above |
| O-information | Signed synergy-vs-redundancy score | Step 4 | [Rosas et al. 2019](https://doi.org/10.1103/PhysRevE.100.032305) |
| PID | Full decomposition into redundant/unique/synergistic | Step 4 | [PID review](https://arxiv.org/abs/1004.2515) |
| Transfer entropy | Directed info flow with history removed | Pipeline | [Schreiber 2000](https://doi.org/10.1103/PhysRevLett.85.461) |
| IBL Brain-Wide Map | 620k-neuron standard dataset | §3 | [IBL data](https://www.internationalbrainlab.com/data) |
| Neuropixels | Probe recording 100s of neurons at once | §3 | [neuropixels.org](https://www.neuropixels.org/) |
| Surrogate data | Fake data keeping some properties, breaking others | §5 | [Surrogate testing](https://en.wikipedia.org/wiki/Surrogate_data_testing) |
| GAMM | Regression with smooth curves + group effects | Step 5 | [mgcv docs](https://cran.r-project.org/package=mgcv) |

## 7. What would change our minds

We drop the coupling hypothesis if, with adequate power: Ω vs
distance-to-criticality is flat or monotone; the relation flips sign
under subsampling or bin-size changes; or rate-matched surrogates
reproduce the "effect". Any of these means criticality and synergy are
separate stories — a clean, useful negative result.

## 8. Try it yourself

```bash
pip install -r requirements.txt
export PYTHONPATH=src
python -m unittest discover -s tests -v   # 16 tests, fully offline
```

The suite simulates branching-process "brains" with a known σ and checks
we recover it; computes O-information on XOR (must be −1 bit) and clones
(must be positive); and verifies transfer entropy points the right way on
a copy-chain. One second of compute, all the core machinery proven.
