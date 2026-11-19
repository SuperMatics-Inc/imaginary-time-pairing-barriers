# STAGE3 | New method barriers, not a linear pairing ceiling

**Verdict: PROVED analytical method barriers and positive lemmas; the fixed-density linear or linear-times-polylog pairing ceiling remains OPEN.** Averaging the phonon histories does not preserve Wick's rule: the exact Holstein atomic connected cross-spin occupation cumulant is `tanh(βU/4)/4`, whereas Wick applied to the averaged covariance predicts zero. The full averaged grand-canonical two-point kernel nevertheless has a chronological determinant bound **one**; it describes an auxiliary Gaussian state and misses the physical interacting vertices. On a side-four square SSH torus at balanced density **1/2**, the exact fermion-integrated action has a **negative nonzero-Matsubara Hessian** whenever `τ≤U/32` and `τ³≤UΩ²/1024`. Hence a uniform convexity/Gaussian-domination argument for that action cannot work at any `τ=CΩ[log(1/Ω)]^p`, with fixed C>0 and p≥0. These are method-specific obstructions, not physical counterexamples. Independent mathematical review of the new analytic proofs is pending.

## 1. New statements and their precise scope

All proofs are self-contained in [RIGOR.md](RIGOR.md).

| Statement | Label | Proof |
|---|---|---|
| Exact interacting atomic cumulant `tanh(βU/4)/4`; history averaging and determinants do not commute | PROVED method barrier | Proposition 1 |
| Full grand-canonical chronological electron moments bounded by product of insertion norms | PROVED positive lemma | Proposition 2 |
| Lehmann auxiliary-CAR realization of the averaged two-point kernel; all chronological covariance determinants bounded by one, including PSD vertex interpolation | PROVED positive lemma; **not** physical Wick factorization or a claimed universal time-independent Gram representation | Proposition 3 |
| `log E exp(c∫|v|)≥c²U/[164Ω]` at βΩ=1 in the exact atomic dynamic process; more generally `≥c²UΩ²/(320τ³)` for βΩ≤2π | PROVED barrier to paying this particular absolute moment as the forest/activity majorant | Proposition 4 |
| Negative SSH slow-field Hessian in a nonzero Matsubara direction under the explicit sufficient conditions above | PROVED barrier to uniform convexity of the exact phonon-only action, including bounded frozen fields and fixed bare t>0 | Proposition 5 |
| Frozen two-electron resolvent bound `(2/τ)exp[−d log(1+τ/(16J))]`, arbitrary static diagonal, no gap required | PROVED positive spectral lemma | Proposition 6 |
| Smooth natural-amplitude good SSH histories can have exponentially amplified intermediate charged two-time pair insertions, with constant electron history weight | PROVED barrier to **pointwise** good-block temporal gluing; not to an annealed insertion estimate | Proposition 7 |
| Exact atomic pair ground-state overlap `exp(−U/Ω)` and low-sideband weight `exp(−2U/Ω) Σ_(m≤floor A)(2U/Ω)^m/m!` up to energy AΩ | PROVED atomic spectral lemma; not a lattice tunnelling or stiffness bound | Proposition 8 |

The conditional Hessian obstruction has a two-thirds-power scale `(UΩ²)^(1/3)`. That scale is **not a new physical temperature ceiling**. None of these results changes the released square-lattice square-root upper bound.

## 2. The focused annealed-block pivot

The gap-free frozen-field spectral estimate was completed. It works over distances of order J/τ, retaining arbitrarily low electron energies. It is not enough to estimate exact interacting charged transport across blocks: the interacting four-point function is not a determinant of the averaged two-point kernel, and pair susceptibility is a bosonic zero-frequency quantity rather than automatically the two-electron resolvent at iτ.

A smooth SSH sine history, with static hopping canceled and amplitude `ρsqrt(UΩ)`, is good for the natural slab cutoff. Its full-period matrix is I and its electron history weight is constant, but an intermediate charged pair insertion in the balanced sector is

```
(1−k/V)² sinh²[(ρa/π)sqrt(U/Ω)],   β=a/Ω.
```
This defeats the proposed pointwise good-set temporal gluing substitute. It does not defeat an annealed estimate that integrates the exact electron weight and boundary data before absolute values. The source-free joint excursion probability bound does not provide that insertion-weighted good-set estimate.

**Precise remaining lemma:** an insertion-sensitive annealed spatial/temporal block inequality, retaining low-energy fermions and their open boundary data, with a summable charged attachment cost. No such lemma was proved in this focused attempt. The unconditional absolute integrated-field moment would restore an exponentially large cost at the proposed linear-times-polylog scale, so it cannot simply be substituted for the missing inequality.

## 3. Physics sanity check

The atomic Franck–Condon calculation genuinely suppresses the low-energy pair spectral weight exponentially in U/Ω. **[INFERENCE]** A similar coherent-transport suppression could make the true coherence scale smaller than linear in a regime with separated adiabatic distortion minima and persistent tunnelling barriers. Those regime hypotheses and their finite-density, fixed-t transport consequence remain unproved.

The parent LINEAR2 result forbids promoting a vacuum overlap to an exact kinetic or stiffness bound: at regular density 1/2 its weak-coupling kinetic expectation stays of order t uniformly in Ω, and its exponentially small pair-transfer quantity is only a t=0 Taylor coefficient with soft-phonon gaps in the complement. No universal exponential pairing ceiling, polaron mass estimate, or physical violation of the proposed linear bound is claimed here. All numerical values in this stage are calculated, not measured.

## 4. Executed certificate and mutants

From `ADIABATIC/LINEAR5/STAGE3/`, the prescribed offline laptop command was executed:

```
UV_OFFLINE=1 UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with numpy --with scipy --with mpmath python -B replay.py --all-mutants
```

**Exit 0; reported wall time 3.29 seconds.** It ran the base certificate and all twelve child commands `<offline_uv_python> -B replay.py --mutant NAME`. Every child exited **1** at its named guard:

| Mutant | Exit | Rejection guard |
|---|---:|---|
| `wick-average` | 1 | averaging-does-not-preserve-wick |
| `canonical-wick` | 1 | canonical-fixed-number-is-not-gaussian |
| `lehmann-weight` | 1 | lehmann-sum-not-difference |
| `thermal-normalization` | 1 | thermal-normalization-bounds-negative-energy-branch |
| `chronological-sign` | 1 | chronological-antiperiodic-sign |
| `fourier-variance` | 1 | two-conjugate-fourier-modes-and-beta-factor |
| `ssh-canonical` | 1 | sixteen-site-canonical-binomial-weights |
| `ssh-spin` | 1 | balanced-two-spin-susceptibility-factor |
| `ssh-kinetic` | 1 | single-bond-cosine-kinetic-half-factor |
| `hessian-sign` | 1 | negative-nonzero-matsubara-hessian |
| `good-block` | 1 | good-amplitude-flag-does-not-control-charged-gluing |
| `franck-condon` | 1 | franck-condon-amplitude-versus-spectral-weight |

**CERTIFIED finite evidence:** exact rational Machin/π bounds and sufficient-condition constants, canonical binomial coefficients, kinetic and spin factors, and 70-digit mpmath interval spot certificates. Among the observed intervals:

- Atomic cumulant at βU=4: `0.1903985389889412220... > 0`.
- SSH `S_eff''/β` at U=1, τ=1/64, Ω=1/16: `−3.3308489377064198424... < 0`; this point satisfies both sufficient conditions, including the cubic boundary equality.
- Atomic pair sideband weight for m≤2 at U/Ω=8: `0.00001631760033429257...`; the exact expression is `145 exp(−16)`.

The replay does not machine-prove the all-order chronological determinant argument, the infinite-dimensional Hölder proof, every Hessian parameter, or any annealed pairing bound. Independent mathematical review, Lean, TeX, repository builds/lint/tests and broader checks were **not run**. The requested small offline replay is the only executed verification. No network, publishing, Git operation, branch, commit, or paid/cloud compute occurred; spend $0.

## 5. Files and next decision

Created only under `ADIABATIC/LINEAR5/STAGE3/`:

- [RESULT.md](RESULT.md): verdict and observed replay evidence.
- [RIGOR.md](RIGOR.md): full proofs of all eight propositions and explicit limitations.
- [replay.py](replay.py): finite certificate and all retained mutant modes.
- [SUMMARY.md](SUMMARY.md): consolidated status of parent §§12–18 and this stage.

Protected parent §§1–18 and the parent replay were not edited. No throwaway files were created. The positive lemmas and barriers are ready for independent review; an end-to-end linear, logarithmically corrected linear, or tunnelling-based fixed-density pairing ceiling remains OPEN.
