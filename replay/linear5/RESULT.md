# Linear adiabatic target: positive machinery and an exact reference barrier

**Verdict: PROVED interacting Gaussian domination and a natural-amplitude slab budget, with PROVED method barriers; the linear ceiling remains OPEN (analytical proofs, independent review pending).** For both finite-volume Holstein and bond-SSH models, including every balanced canonical sector `N_up=N_down=k`, `E exp[η(b_e(s)−b_e(0))]≤exp[UΩ²η²s(β−s)/(2β)]` under the full interacting positive measure. Disjoint same-oscillator slabs of width `ℓ/Ω` have joint excursion probability at most `[4exp(−ρ²/(32ℓ))]^m` at amplitude `ρsqrt(UΩ)`, for `ρ²≥512ℓ`, uniformly in volume and external static parameters, but not conditional on a prescribed zero mode. This positive linear-time budget does not close the LINEAR4 forest: a continuous palindromic SSH reference has normalized unequal-time attachment norm `exp[(ρa/8)sqrt(U/Ω)]/2` at `β=a/Ω`. Small field increments against fixed t also do not justify instantaneous Fermi-equilibrium tracking; the occupation error tends to `1/2`. Ordinary spatial phonon-coordinate reflection positivity fails in fixed-density balanced Holstein sectors, at any t, because the reflected zero-mode form is `−U²Var_D(n_A)<0`. These are precise auxiliary-method barriers, not physical pairing-temperature counterexamples. The earlier imaginary-time RP, fugacity comparison, and atomic order-Ω source lemma remain valid; no `k_B T_pair≤CΩ` or Ω-times-log theorem has been obtained.

## 1. Context and conventions

The open target is the pure spinful square-lattice Holstein and bond-SSH problem at fixed spring constant, force, bare hopping and positive regular limiting density. Here `U=F²/K` is a phonon coupling energy, not a Hubbard repulsion; `Ω=ℏsqrt(K/m)` and `τ=k_B T` are energies; `β=1/τ`. The asymptotic parameter is `Ω→0` with `t,U>0` fixed. No phonon-vacuum projection is made.

The released parent work proves the annealed static no-pair-order endpoint and a square-root upper scale at positive ion mass. LINEAR2 rules out using vacuum hopping as equilibrium kinetic suppression, simultaneous canonical SSH displacements, and a reference-uniform positive tree estimate as a route to a linear scale. LINEAR3 proves the fixed-width cylinder result and a signed bulk quantum/static free-energy sandwich, but explains why it does not control the square-lattice winding response or a sourced quadratic interpolation coefficient. These are method obstructions, not physical lower bounds on a transition temperature.

The two-thirds-power claim in LINEAR4 is under independent review, particularly its fermion-weighted tail replacement. Nothing in LINEAR4 was edited in this pass, and its theorem is not assumed here. The present exact counterexample uses temporal amplitudes outside its small-reference cutoff, so it does not itself disprove that manuscript.

No new third-party literature or paper content was used. All new statements below are derived in [RIGOR.md](RIGOR.md). Parent context is [ADIABATIC](../RESULT.md), [LINEAR2](../LINEAR2/RESULT.md), [LINEAR3](../LINEAR3/RESULT.md), and the released companion proof under `PUBLISH/adiabatic-merged/repo/RIGOR.md`.

## 2. Positive direction: reflection in imaginary time

**PROVED, analytical, independent review not yet performed.** Let `s=β/2`, `K_s(x,y)` be the positive harmonic heat kernel and `ν_xy` the corresponding normalized half-period bridge law. For a half path `w`, let `M_w` be its time-ordered electronic evolution, in the full Fock space or one fixed-number sector. Define `Q_f(x,y)=∫ν_xy(dw) f(w)M_w`. Then

```
Z⟨(Θf)g⟩ = ∫ dx dy K_s(x,y)² Tr[Q_f(x,y)†Q_g(x,y)].
```

Thus the full interacting functional satisfies reflection positivity and the associated Schwarz inequality. For a grand-canonical half-path event `E`, where the determinant-square measure is positive,

```
⟨1_E⟩² ≤ ⟨1_E (Θ1_E)⟩.
```

The right side retains the interacting electronic weight. No claim that a rare event under the free Gaussian is rare under the interacting measure enters the proof. Noncommuting SSH bond Hamiltonians pose no obstruction: reversal takes the adjoint of the ordered product. The same functional identity is valid in canonical sectors, without assuming pointwise canonical path positivity.

Every mirrored one-particle reference has period matrix `m_w†m_w>0`, so `det(I+q m_w†m_w)>0` for every real `q>0`. This meets the request to construct references that remain invertible. The reflection comparison is established; an Ω-dependent estimate of its reflected restrictions is not.

**Precise stopping lemma.** Bound the reflected **interacting** event or insertion functional at `βΩ` of order one strongly enough to make a spatial charged-correlation expansion summable, or bound the mirrored/original connected source-interpolation difference at order Ω. The Schwarz inequality itself contains no small Ω factor. A chessboard iteration with a sufficiently small weighted restriction budget remains OPEN.

## 3. Positive direction: fugacity averaging with pressure equivalence

**PROVED, analytical.** Write the one-spin history polynomial as `D(z)=Σ_(n=0)^V a_nz^n`, with real coefficients and `a_0=1`. Opposite-spin angular twists yield

```
W_bal = (1/(2π))∫ |D(ze^(iθ))|²dθ
      = Σ_n a_n²z^(2n) ≥ 1.
```

After integrating the phonons, this is exactly the physical partition function restricted to `N_↑=N_↓`. Pointwise Cauchy–Schwarz proves `Z≤(V+1)Z_bal`; positivity of the physical spin-sector quantum traces proves `Z_bal≤Z`. Their pressure difference is at most `τ log(V+1)/V`.

This is a genuinely positive, nonzero unnormalized history reference with the correct thermodynamic grand-canonical pressure. It is not a free-Gaussian substitute and is not merely a chosen invertible twist. It does **not** establish equality of fixed-density sourced correlations, which would require additional work.

**PROVED limitation.** The same exact SSH loop described below has `W_bal=4`. Its conditional angular-averaged onsite pair transfer on one edge is `200/9`, already exceeding a physical state norm bound. Conjugating that loop by extra allowed bond pulses makes some pair transfers arbitrarily large while leaving its spectrum and `W_bal` unchanged. The conditional reference is not a positive electronic density matrix. Weighted spatial determinant/minor estimates are still needed; strict positivity of the averaged denominator alone does not control nonnormality.

## 4. Exact unrestricted-reference barrier

**PROVED, analytical; finite algebra CERTIFIED.** On a three-site nearest-neighbor path, let `K_12,K_23` be its real symmetric bond matrices. Set

```
A=exp[(log 3)K_12],       B=exp[log(4+sqrt(15))K_23].
```

Four continuous pulses, with equal positive/negative time integrals on each bond, realize

```
Y=A B A^(-1) B^(-1),
Y²=I,      tr Y=−1,      spec(Y)={1,−1,−1},
det(I+qY)=(1+q)(1−q)².
```

The pulse profiles are explicit continuously differentiable periodic functions. They have zero temporal mean, and their peak coefficients are `8τ log 3` and `8τ log(4+sqrt(15))`. Embed the path into the square SSH model by setting all static bond coordinates to `b_0=t`, and all other temporal coefficients to zero.

At `q=1` the conditional electronic trace vanishes. At `q=1+1/n` a unit-norm one-particle Wick contraction has magnitude `n+1`. At any other fixed positive `q≠1`, conjugating the loop by increasingly large continuous bond pulses creates unbounded contractions through nonnormality. This kills an unrestricted pointwise-reference Wick/minor estimate even before a multiscale spatial argument starts.

The barrier evades the earlier positive-tree obstruction by being an exact electronic counterexample to its putative unrestricted reference premise. It is not an argument based on positive absolute tree weights. It neither deletes the zero Matsubara mode nor identifies a vacuum overlap with equilibrium kinetic energy. The coordinate construction stays entirely within actual SSH bonds and does not invoke a nonexistent simultaneous canonical transformation.

**Not asserted:** that the exact path has positive Gaussian measure, that nearby large contractions survive the interacting weighted integral, or that pair order occurs at this scale. Vanishing determinant weights might cancel large normalized contractions. Such weighted cancellation is exactly what a viable replacement method must preserve.

## 5. Atomic test of the linear scale

**PROVED partial lemma, not a lattice ceiling.** For the isolated Holstein site (`t=0`), use `H_0=Ωa†a+g(a+a†)(n−1)`, `g²=UΩ/2`, and source `−h(c†_↑c†_↓+c_↓c_↑)`. Define `χ=∂_h²[β^(-1)log Z(h)]` at zero source. At each fixed `a=βΩ>0` and `U>0`,

```
χ_qu = 1/U + Ω coth(a/2)/(2U²) + o(Ω),
χ_st = 1/U + Ω/(aU²) + o(Ω).
```

The full exact two-insertion formulas, oscillator trace derivation, static Gaussian integral and controlled endpoint asymptotics are in RIGOR §4. No external paper formula is imported. The coefficient of `χ_qu−χ_st` is strictly positive and finite for fixed `a`.

This calculation tests the positive source route at `τ=Ω/a`. It shows that an order-Ω **connected source correction** is compatible with the exact dynamics in a nontrivial soluble limit; the bulk pressure sandwich does not force a square-root or two-thirds-power susceptibility loss. It does not establish that the linear ceiling is the true answer on a square lattice. An isolated atom has no spatial pairing transition, the limiting susceptibility `1/U` is not a superfluid stiffness, and the expansion is not uniform when `a` is allowed to tend to zero or infinity with Ω.

## 6. What remains OPEN

1. An order-Ω connected lattice source-interpolation estimate at fixed `βΩ`, with uniform spatial/momentum and volume control.
2. A pairing-sensitive summability estimate that can use the now-proved interacting natural-amplitude slab budget. Ordinary fixed-density spatial RP and direct substitution into the normalized Wick forest are blocked by §§14 and 17.
3. A weighted spatial insertion-minor estimate in the balanced fugacity-averaged reference that handles nonnormal histories.
4. Transfer of a successful sourced or charged-correlation estimate to every required fixed-density thermodynamic state.
5. A physical counterexample to the proposed linear ceiling, if such a regime exists. None was constructed.

These are mathematical gaps, not inaccessible files or missing literature. The new §§12–15 settle the positive interacting slab budget and its balanced fixed-density transfer, in addition to the earlier invertibility and imaginary-time reflection lemmas; the pairing-sensitive spatial step remains OPEN.

## 7. Earlier reproducible evidence

The prescribed laptop-only offline command was run from `ADIABATIC/LINEAR5/`:

```
UV_OFFLINE=1 UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with numpy --with scipy --with mpmath python -B replay.py --all-mutants
```

**Observed exit code: 0.** Runtime reported 1.32 seconds. The exact certificate uses `fractions.Fraction` and arithmetic in `Q(sqrt(15))`; numpy, scipy and mpmath were not needed by the script. It verifies the boost identities, monodromy, fugacity polynomial, negative eigenspace projector, occupation divergence, spin-sector coefficients, atomic endpoint factors, opposite-spin Fourier projection, the conditional pair-transfer factor and finite Hilbert–Schmidt Gram identities. The infinite analytic steps in RIGOR are not machine-formalized by these finite checks.

The same invocation executed these child commands, each through the offline environment's Python executable, and read their exits:

| Child command suffix | Observed exit | Rejected guard |
|---|---:|---|
| `python -B replay.py --mutant field-radicand` | 1 | second boost hyperbolic identity |
| `python -B replay.py --mutant inverse-sign` | 1 | monodromy trace |
| `python -B replay.py --mutant commuting-order` | 1 | monodromy trace |
| `python -B replay.py --mutant thermal-occupation` | 1 | divergent reference occupation |
| `python -B replay.py --mutant canonical-positivity` | 1 | conditional canonical sign |
| `python -B replay.py --mutant one-spin-weight` | 1 | two-spin square weight |
| `python -B replay.py --mutant endpoint-factor` | 1 | atomic first-order factor |
| `python -B replay.py --mutant same-spin-twist` | 1 | opposite-spin Fourier projection |
| `python -B replay.py --mutant mirror-without-adjoint` | 1 | reflected Hilbert–Schmidt positivity |

The earlier exploratory exact-rational boost calculation used the same offline UV command with a Python heredoc and exited 0. Directory creation via `mkdir -p .../ADIABATIC/LINEAR5` exited 0. Those operations did not simulate a many-body system.

**Files written:** `ADIABATIC/LINEAR5/RESULT.md`, `RIGOR.md`, `replay.py`, and the requested self-contained `SHARPNESS.md`. All twenty-four mutants are retained as replay modes. No throwaway files were left, and LINEAR4 was not edited. No Git operation, branch or commit was created. No network publishing, external literature intake, many-body simulation, cloud compute or paid compute occurred; spend $0. Skipped: independent review of the new §§12–17, Lean formalization, TeX build, repository builds, lint and broader tests. SHARPNESS.md has separately been submitted to independent review by Main. The requested exact replay and mutants are the only local verification.

## 8. New positive route and its exact stopping point

### Full interacting temporal Gaussian domination

RIGOR §§12–13 prove the source comparison by conjugating only the phonon kinetic energy. For a mean-zero temporal coordinate source with piecewise-constant gauge q, its exponential moment is bounded by `exp[(UΩ²/2)∫|q−q̄|²]`. Schatten Hölder controls the exact interacting semigroups; no positive conditional fermion covariance or free-Gaussian replacement is assumed. The continuous-slab estimate is obtained by a dyadic union bound and applies jointly, without independence, to slabs disjoint for each oscillator.

At width `ℓ/Ω` and excursion amplitude `ρsqrt(UΩ)`, the budget is an Ω-independent cost per actual slab. This escapes the `κτ` small-cutoff hypothesis in SHARPNESS.md. It does not assert a doubled cost `βΩ` per individual slab. For the balanced sector, the history coefficient is the square of the real one-spin canonical coefficient, so the measure remains positive at prescribed density. RIGOR §15 proves both the balanced fixed-density budget and `Z_(N/2,N/2)≤Z_N≤(N+1)Z_(N/2,N/2)` for even N. The latter is a zero-source pressure comparison, not by itself a pair-correlation transfer.

### The forest attachment cost is still uncontrolled

RIGOR §14 supplies a natural-amplitude, zero-mean, palindromic one-bond SSH reference with period matrix I but a large unequal-time contraction. Thus reflected reference invertibility is insufficient. A factorized slab rarity bound cannot simply replace the small-history estimate in the old all-insertion-time normalized Wick forest. An unnormalized insertion estimate or a reflection-paired operator comparison may still evade this barrier.

### Two further routes tested

**Adiabatic multiscale: PROVED barrier to uniform instantaneous-equilibrium tracking, not all multiscale methods.** The same reference has `|v'|≤4πr/β`, so its increments across windows of width `1/t` are of order `Ω^(3/2)`. Nevertheless its history occupation is `I/2` whereas the instantaneous Fermi occupation differs by `tanh(βr/2)/2→1/2`. Static fields can cancel the bare hopping, so bare t supplies no uniform instantaneous gap. A tracking scheme must retain imaginary-time memory and handle gap crossings, or justify an annealed exception estimate.

**Chessboard on positive balanced sectors: PROVED barrier to ordinary spatial RP in Holstein.** Positivity alone does not imply RP. RIGOR §17 proves the exact source identity `E exp(ξ·b_0)=exp[(Uτ/2)|ξ|²]Z(H+Uτξ·n)/Z(H)` and the resulting negative reflected covariance for spatially reflected halves. This failure holds for arbitrary fixed t>0 and nonempty, nonfull balanced canonical density, not merely in an isolated atom. It does not exclude modified or grand-canonical spatial comparisons, or imaginary-time Gaussian domination. No phonon-only canonical SSH spatial-RP conclusion is claimed.

### What remains OPEN

1. A pairing-sensitive insertion/transport estimate that survives natural-amplitude histories without the exponentially large arbitrary-time normalized Wick factors.
2. A successful adiabatic tracking estimate through static gap closings, or an integrated estimate that does not require uniform local thermal equilibrium.
3. A spatial comparison bypassing the fixed-density Holstein RP obstruction.
4. The source/ensemble argument needed to turn any completed spatial bound into a full fixed-total-density electronic pairing ceiling. Equality of zero-source free energies does not supply it.

Uniformity of the proved gauge bound is in external static parameters of H, not in a conditional law fixing the zero Matsubara coordinate. No new third-party references, measured quantities, or physical ordering-temperature lower bounds were introduced.

## 9. Observed replay evidence for the new results

After the new proofs and replay modes were written, the prescribed command from §7 was executed from `ADIABATIC/LINEAR5/` and exited **0**, with reported wall time **6.68 seconds**. It ran all twenty-four retained child mutants; **every child exited 1** at its intended rejection guard. The first fourteen include the complete SHARPNESS.md replay. The ten new modes were:

| Mutant command suffix | Observed exit | Rejected guard |
|---|---:|---|
| `--mutant kinetic-gauge-sign` | 1 | imaginary coordinate gauge square sign |
| `--mutant gauge-centering` | 1 | optimal circle gauge centering |
| `--mutant dyadic-allocation` | 1 | dyadic telescoping threshold budget |
| `--mutant natural-amplitude-power` | 1 | linear-scale natural-amplitude budget |
| `--mutant mirror-equal-time` | 1 | mirrored unequal-time minor growth |
| `--mutant drift-unitarity` | 1 | anti-Hermitian drift Schatten point |
| `--mutant balanced-density-weight` | 1 | fixed balanced positive history weight |
| `--mutant canonical-multiplicity` | 1 | fixed-density SU(2) multiplicity bound |
| `--mutant adiabatic-occupation` | 1 | instantaneous-equilibrium tracking error |
| `--mutant spatial-RP-sign` | 1 | balanced canonical spatial RP counterexample |

Each child was `python -B replay.py --mutant NAME` under the offline UV environment's Python executable. Exact-rational checks certify the gauge constants, zero-mode point factor, dyadic thresholds and series ratio, natural-amplitude power, mirrored contraction at `η=log 4`, balanced-sector spin factors, tracking error at `βr=log 4`, and the atomic negative spatial form at `βU=log 2`. A 60-digit mpmath interval check certifies a finite `2×2` anti-Hermitian-drift Schatten `S_2` example. The seven earlier interval covariance points were also retained.

**Certification scope:** these finite checks are not machine proofs of the unbounded semigroup comparison, continuous-time joint-event chaining, infinite-volume pairing statement, or arbitrary-t canonical spatial-RP lemma. The complete analytical proofs are in RIGOR §§12–17; independent mathematical review remains pending. No additional test suites, Lean checks or TeX builds were run.
