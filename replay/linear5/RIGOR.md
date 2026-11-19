# Unrestricted dynamic references and the linear adiabatic target

## 1. Verdict and scope

**PROVED method barrier.** A pointwise, history-uniform normalized Wick bound does not exist for arbitrary continuous, real, periodic, spin-independent bond-SSH histories, even on a three-site path and with zero temporal mean. Four allowed bond pulses give an electronic one-period monodromy with eigenvalues `1,−1,−1`. Its one-spin grand-canonical trace at chemical potential zero is zero. Nearby real chemical potentials give normalized one-particle contractions of arbitrarily large magnitude. This refutes a multiscale extension that simply admits unrestricted temporal histories into the quadratic reference and retains the normalized pointwise Wick estimates. It is not a counterexample to a thermodynamic pairing-temperature bound.

**PROVED partial lemma.** In the isolated Holstein atom, with source coupling the empty and doubly occupied states, the quantum correction to the zero-source pair susceptibility is of order `Ω` at every fixed positive ratio `βΩ`. Its exact coefficient is given in §4. This supplies an exactly soluble check for a connected source-interpolation estimate, not such an estimate on a lattice with fixed hopping.

**PROVED positive lemmas.** The full finite-volume functional is imaginary-time reflection positive in both models, grand canonically and in each fixed particle-number sector (§5). Opposite-spin fugacity averaging gives a strictly positive unnormalized history weight and a balanced-spin partition function satisfying `Z_bal≤Z≤(V+1)Z_bal` (§6). Neither lemma alone contains the required Ω-dependent spatial estimate. Mirrored references stay invertible; twist averaging alone still leaves unbounded conditional pair transfers.

**PROVED sharper restriction barriers (§§8–11).** Removing or conditioning on the zero Matsubara coordinate destroys imaginary-time reflection positivity, already in the exact interacting Holstein atom. If a positive rare-history comparison retains the small-amplitude reference cutoff `sup|v|≤κτ`, its necessary and sufficient atomic restriction scale is `τ` of order `(UΩ²)^(1/3)`, not Ω. This is a sharp cutoff-probability barrier for the stated method class, not a square-lattice physical ceiling or a no-go for every reflected reference.

**PROVED positive interacting estimate and further precise barriers (§§12–17; independent review pending).** A kinetic coordinate gauge and Schatten Hölder give `E exp[η(b_e(s)−b_e(0))]≤exp[UΩ²η²s(β−s)/(2β)]`, and a factorized natural-amplitude slab budget, under the full grand-canonical or fixed balanced-spin interacting measure. The constants are uniform in volume and external static parameters, not conditional on a prescribed zero Matsubara coordinate. Balanced fixed-density transfer and its zero-source pressure comparison are proved. A palindromic SSH reference still has arbitrary-time normalized Wick contraction of order `exp[c sqrt(U/Ω)]`, and its instantaneous thermal occupation need not track the slow field. Ordinary spatial phonon RP also fails in fixed-density balanced Holstein sectors, at arbitrary t, with reflected form `−U²Var_D(n_A)<0`. These are method-specific barriers, not physical counterexamples; the linear pairing ceiling remains OPEN.

**OPEN.** A square-lattice upper bound `k_B T_pair≤CΩ`, or `Ω` times a logarithm, for the pure spinful Holstein and bond-SSH models at fixed positive regular density and fixed `t,U>0`. No physical family violating that bound is exhibited. The LINEAR4 manuscript remains under independent review and is not edited or relied on as a released theorem here.

## 2. Exact continuous SSH histories with a zero conditional trace

Write `K_12` for the real symmetric `3×3` matrix with unit entries at `12` and `21`, and `K_23` for the corresponding `23` and `32` matrix. They are allowed one-particle nearest-neighbor bond bilinears on a path. If a history on one time interval has Hamiltonian `−f(u)K_ij`, its electron evolution on that interval is `exp(K_ij∫f)`, since all generators on that interval are multiples of the same matrix. No commutation between distinct intervals is assumed.

Set

```
η_1 = log 3,               η_2 = log(4+sqrt(15)),
c_1 = cosh η_1 = 5/3,      s_1 = sinh η_1 = 4/3,
c_2 = cosh η_2 = 4,        s_2 = sinh η_2 = sqrt(15).
```

The identities follow immediately from `3^(-1)=1/3` and `(4+sqrt(15))^(-1)=4−sqrt(15)`. Put `A=exp(η_1 K_12)`, `B=exp(η_2 K_23)`. Thus

```
A = [5/3  4/3  0;  4/3  5/3  0;  0  0  1],
B = [1  0  0;  0  4  sqrt(15);  0  sqrt(15)  4].
```

Both are positive definite symmetric matrices. Their inverses negate the two off-diagonal entries in the active `2×2` block.

To realize the products by continuous zero-mean histories, divide `[0,β]` into four intervals of length `β/4`. On interval `j`, with left endpoint `l_j`, let

```
f_j(u) = (8/β) sin²(4π(u−l_j)/β),
```

and set it to zero outside the interval. Each `f_j` is continuously differentiable on the period, vanishes with its derivative at its endpoints, and has integral one. Assign chronological Hamiltonians

```
+η_2 f_1 K_23,  +η_1 f_2 K_12,
−η_2 f_3 K_23,  −η_1 f_4 K_12.
```

Each bond coefficient has zero time average. The one-period evolution, with later factors to the left, is exactly

```
Y = A B A^(-1) B^(-1)
  = [−13/3   20/3   −4sqrt(15)/3;
     −20/3   37/3   −8sqrt(15)/3;
     −4sqrt(15)/3   8sqrt(15)/3   −9].                 (1)
```

Multiplication in `Q(sqrt(15))` gives `Y²=I`, `tr Y=−1`, and `det Y=1`. For example the `11` entry of the square is `(169−400+240)/9=1`; the off-diagonal `12` entry is `(−260+740−480)/9=0`; all other entries follow by the same multiplication. The minimal polynomial divides `(z−1)(z+1)` with distinct roots, so `Y` is diagonalizable. Its trace implies that the `+1` eigenspace has dimension one and the `−1` eigenspace dimension two. In particular

```
det(I+qY) = (1+q)(1−q)²                              (2)
```

for every scalar `q`.

The conditional one-spin electronic trace is `det(I+Y)`. To see this without any positivity assumption, second quantization of the number-conserving one-particle evolution is its exterior lift: on the `N`-particle space it is `∧^N Y`. The trace of its direct sum is the sum of elementary symmetric functions of the eigenvalues, equal to `det(I+Y)`. This holds for products of quadratic evolutions by multiplicativity of the exterior lift, not merely for a static Hermitian Hamiltonian. Two identical spins give the square of this real determinant. Equation (2) therefore gives a zero conditional integrand at `q=1`.

These histories are allowed in the square bond-SSH model. Choose a three-site path in a torus with side at least three. Set every static bond energy coordinate to `b_0=t`, so all reference hoppings `−t+b_0` vanish. On the two active edges take the continuous temporal coefficients above, with the appropriate sign relative to the one-particle Hamiltonian; all other temporal coefficients are zero. The remaining sites are spectators. Their one-spin trace is `(1+q)^(V−3)`, and they do not remove the zero in (2). The physical coordinates remain finite, continuous, and periodic. There is no phonon-vacuum projection, ionic-mass limit, or approximation in the monodromy. For sufficiently small `τ` at fixed `t,U>0` these entire histories are inside the auxiliary static-size scale `4096(t+U+τ)` if that scale is retained. Their temporal amplitude is at most `8τ max(η_1,η_2)`, much larger than `τ/64`; the construction is therefore outside, and does not invalidate, a small-amplitude dynamic reference.

The argument is a counterexample to a **pointwise uniform estimate over the histories**. It does not require the exactly specified path to have positive Gaussian measure, and makes no claim that it has typical interacting weight.

## 3. Divergent Wick contractions, weighting, and the method boundary

Include a real chemical potential in the one-particle Hamiltonian. Since it multiplies identity, the period matrix becomes `qY`, `q=exp(βμ)>0`. Whenever `q≠1` the conditional trace is nonzero. Its normalized one-particle contraction at the beginning/end of the period is

```
N_q = qY(I+qY)^(-1)
    = [q/(1+q)] P_+ + [q/(q−1)] P_-,
P_+ = (I+Y)/2,           P_- = (I−Y)/2.                (3)
```

For completeness, the normalized exterior trace of `c†(v)c(w)` is `w* qY(I+qY)^(-1)v`, with the transpose or complementary contraction depending on the creator/annihilator convention. Diagonalize `Y` and sum the independent orbital occupation weights `1, qλ`; this gives `qλ/(1+qλ)` for a diagonal orbital. Under a one-particle similarity the creation and dual annihilation vectors transform contragrediently. Continuity extends the formula to all nonsingular `I+qY`. Thus the large eigenvalue in (3) occurs in an actual one-body Wick contraction, not in an unrelated operator inverse.

For a unit right eigenvector of `P_-`, the quadratic form of `N_q` has magnitude `|q/(q−1)|`, since `N_q v=[q/(q−1)]v`. In particular `q=1+1/n` gives the value `n+1`. Unit-norm creation and annihilation operators are being used, but their conditional evolution is not a positive normalized density matrix, so the ordinary thermal bound by one is inapplicable. Consequently even a `1×1` minor has no finite history-uniform bound. Finite-volume matrix differential equations, their evolution matrices, and their inverses where defined depend continuously on the coefficients. A fixed `q≠1` sufficiently close to one therefore has a neighborhood of this continuous history with contractions above any prescribed finite bound.

There is also no history-uniform repair at a fixed `q≠1` merely by removing the singular point. Add two bond pulses to conjugate the monodromy by `T_γ=exp(γK_12)`. The conjugated history still has zero integral in every bond, because the added pulse and its inverse cancel in the time average. Its period matrix is `T_γ Y T_γ^(-1)`, and the corresponding contraction is `T_γ N_q T_γ^(-1)`. The entries of this matrix grow without bound with `γ`. To check this, decompose into the `±1` eigenspaces of `K_12` and its zero eigenspace. The `+` to `−` matrix element of `Y` in the orthonormal vectors `(1,1,0)/sqrt(2)` and `(1,−1,0)/sqrt(2)` is

```
(Y_11−Y_12+Y_21−Y_22)/2 = −15.
```

After conjugation it is multiplied by `exp(2γ)`. In `N_q` its coefficient relative to `Y` is `−q/(q²−1)`, nonzero for every fixed positive `q≠1`; the identity part contributes zero to this matrix element. Hence its magnitude is `15 |q/(q²−1)| exp(2γ)`. These larger pulses remain finite and continuous for each `γ`, but are not amplitude-bounded uniformly in `γ`. This establishes the barrier for unrestricted histories also at fixed chemical potential away from zero.

The free or electronically weighted path integral is not contradicted. At `q=1` its two-spin zero-source electronic weight vanishes. Near this path the normalized Wick contractions can be large while the squared determinant weight is small. An estimate for **unnormalized** insertion determinants, or an integrated estimate retaining this weight, might exploit exactly that compensation. The present counterexample does not rule out such estimates, a multiscale treatment controlling the dangerous history subsets, or cancellations in connected correlations. It rules out passing from the static reference to an unrestricted history while keeping pointwise normalized covariance/minor bounds.

A related sign warning is exact on the isolated three-site block. Its one-spin fugacity polynomial is `1−q−q²+q³`. For two spins it is

```
(1−q−q²+q³)² = 1−2q−q²+4q³−q⁴−2q⁵+q⁶.             (4)
```

Some conditional fixed-number coefficients are negative. Thus the spin-square positivity of the **grand-canonical** history integrand does not establish positive history weights in each canonical sector. Equation (4) is a finite-block fact, not a statement about the sign of a full square-torus sector after adding spectators, nor about the sign of a physical partition function after integrating the phonons. Fixed-density conclusions need an argument that does not assume conditional canonical positivity.

## 4. An exact atomic Holstein test of the quadratic source route

This subsection concerns `t=0`, an isolated spinful Holstein site, not the fixed-positive-hopping square-lattice target. Use coupling to `n−1`, which is a coordinate shift convention, and subtract the free oscillator free energy. The Hamiltonian is

```
H_0 = Ω a†a + g(a+a†)(n−1),   g²=UΩ/2.
```

Add `−h A`, with `A=c†_↑ c†_↓+c_↓c_↑` in a phase convention where the empty/double matrix element is one. Define pressure `p(h)=β^(-1) log Z(h)` and susceptibility `χ=p''(0)`. Both quantum and static systems at zero source have

```
Z(0)=2[exp(βU/2)+1],
```

because empty/double states have `(n−1)²=1` and displaced energy `−U/2`, whereas the two singles have zero displacement.

Put `a=βΩ>0` and

```
R_a(s) = [cosh(a/2)−cosh(Ω(s−β/2))]/sinh(a/2).
```

The exact susceptibilities are

```
χ_qu = exp(βU/2)/(exp(βU/2)+1)
       × ∫_0^β exp[−(2U/Ω) R_a(s)] ds,                (5)
χ_st = exp(βU/2)/(exp(βU/2)+1)
       × ∫_0^β exp[−2Us(1−s/β)] ds.                 (6)
```

To derive (5), complete the square separately for empty and doubly occupied states. The two displacement parameters are opposite, of magnitude `g/Ω`; flipping the pair therefore inserts a displacement of magnitude `d=2g/Ω`, `d²=2U/Ω`. The free thermal expectation of the two opposite displacement insertions at imaginary-time separation `s` is

```
C(s)=exp[−d² (1−e^(−Ωs))(1−e^(−Ω(β−s)))/(1−e^(−βΩ))].
```

For a self-contained oscillator trace calculation, put `r=exp(−a)` and `n̄=r/(1−r)`. The number-basis thermal probabilities are `(1−r)r^n`, so `⟨(a†)^k a^k⟩=k!n̄^k`, by differentiating the geometric series `k` times. Normal ordering and the scalar commutator `[a,a†]=1` then give

```
⟨exp(v a†+w a)⟩ = exp[(n̄+1/2)vw],
exp(v a†+w a) exp(v' a†+w' a)
 = exp[(wv'−vw')/2] exp[(v+v')a†+(w+w')a].
```

Indeed `exp(v a†+w a)=exp(vw/2)exp(v a†)exp(w a)`, and the normal-ordered thermal series is `Σ_k (vw)^k n̄^k/k!=exp(vwn̄)`. Apply the two formulas to `v=d exp(Ωs)`, `w=−d exp(−Ωs)`, `v'=−d`, `w'=d`, obtained by conjugating the first displacement by `exp(sΩa†a)`. The logarithm of its thermal two-insertion expectation is

```
−d² sinh(Ωs)+d² coth(a/2)[cosh(Ωs)−1] = −d² R_a(s).
```

This proves the displayed `C(s)`, including its equivalent exponential-product form. It checks `C(0)=C(β)=1` and `C'(0)=−2U`. The second derivative of the partition at zero source is `β∫_0^β Tr(e^(−(β−s)H_0) A e^(−sH_0) A)ds`; the even-sector trace is `2exp(βU/2)C(s)`. Dividing by `βZ(0)` gives (5).

For (6), the static coordinate `b` has Gaussian variance `U/β`. A flip from `n−1=1` to `−1` for a duration `s` changes its integrated coefficient from `β` to `β−2s`. The normalized Gaussian moment is `exp[U(β−2s)²/(2β)]`; divide by the zero-source even-sector moment `exp(βU/2)` to obtain the integrand in (6). An equivalent direct calculation uses the even-sector matrix eigenvalues `±sqrt(b²+h²)`: twice differentiating `2cosh(βsqrt(b²+h²))` yields `2β sinh(βb)/b`, with its continuous value at zero. This equals (6) after Gaussian integration and the substitution `u=2s−β`, confirming the source factors independently.

At fixed `U>0` and fixed `a=βΩ>0`, as `Ω→0`,

```
χ_qu = 1/U + [Ω/(2U²)] coth(a/2) + o(Ω),
χ_st = 1/U + Ω/(aU²) + o(Ω),
χ_qu−χ_st = [Ω/(2U²)] [coth(a/2)−2/a] + o(Ω).         (7)
```

Here is a full endpoint justification. Both integrands are symmetric about `β/2`, so integrate twice on `[0,β/2]`. For fixed `s`, Taylor expansion gives

```
(2U/Ω)R_a(s) = 2Us−UΩ coth(a/2)s²+O_a(UΩ²s³),
2Us(1−s/β) = 2Us−(2UΩ/a)s².
```

For `0≤s≤β/2`, the function `R_a(s)/(Ωs)` decreases to `2tanh(a/4)/a>0`: write the numerator as the integral of `sinh(a/2−v)` from `v=0` to `Ωs`, and use the decreasing integrand to see that its average decreases. Its endpoint average is `(cosh(a/2)−1)/[(a/2)sinh(a/2)]=2tanh(a/4)/a`. Thus the quantum exponent is bounded below by a positive fixed multiple of `Us`; the static exponent is bounded below by `Us`. Choose a truncation `L(Ω)=C log(1/Ω)` with a fixed sufficiently large `C`, interpreting the logarithm with a fixed energy unit. The integrals beyond `L` are `o(Ω)` by these exponential bounds. On `[0,L]`, the Taylor remainder bound is uniform because `ΩL→0`; expansion of the exponential has integrated error bounded by a constant times `Ω²` times a polynomial in `L`, which is `o(Ω)`. The first-order integrals are therefore obtained by integrating the expansions against `exp(−2Us)` on `[0,∞)`. The moments are `∫ exp(−2Us)ds=1/(2U)` and `∫s²exp(−2Us)ds=2/(2U)³`, by two integrations by parts. Multiply by the two endpoints. The prefactor in (5)–(6) differs from one by `O(exp(−aU/(2Ω)))`, which does not change the coefficients. This proves (7).

The coefficient of the difference is positive: `tanh x<x` for `x>0`, since its derivative is `sech²x<1` and both vanish at zero, hence `coth x>1/x`. It is finite for every fixed `a>0`. This exactly soluble case has an order-`Ω` connected susceptibility correction, despite the much weaker bulk pressure sandwich. It therefore does not support a universal no-go for the quadratic source route. Neither the value `1/U` nor the correction in (7) is a square-lattice susceptibility estimate at `t>0`.

## 5. Positive route: imaginary-time reflection positivity

**PROVED.** In either finite-volume spinful Holstein or bond-SSH model, at every `β>0` and positive harmonic oscillator frequency, the Euclidean functional is reflection positive across the two time slices `0` and `β/2`. This is true grand canonically and in every fixed total-number sector. Noncommuting SSH bond terms are allowed. The statement is about a reflected functional, not positivity of each canonical history weight.

Here is an explicit proof. Let `K_s(x,y)>0` be the product harmonic-oscillator heat kernel at time `s=β/2`. It is symmetric in its endpoints. Write `ν_xy` for its normalized oscillator bridge law on this half period, from coordinate vector `x` to `y`. For a half path `w` let

```
M_w = T exp[−∫_0^s H_el(w(u))du],
Q_f(x,y) = ∫ ν_xy(dw) f(w) M_w,
```

where the real chemical potential can be included in `H_el`. Choose the electron Hilbert space to be the full Fock space, or its fixed total-number sector. The zero-source electronic Hamiltonian preserves that sector and is Hermitian for each real coordinate vector. An imaginary-time partition function is

```
Z = ∫ dx dy K_s(x,y)² Tr[Q_1(x,y)† Q_1(x,y)].        (8)
```

To prove the formula, condition a full oscillator loop on its two endpoints `x,y`. Its halves are independent bridges. Reflecting the second half reverses its order, so its conditional electronic evolution is `M_v†`, since all individual infinitesimal factors `exp[−du H_el]` are Hermitian. The unnormalized integrand is `Tr[M_v†M_w]`. This identity holds for the time-ordered evolution as the limit of time-sliced products; no commutation of different `H_el` is used. Integrating the two independent half bridges gives (8). Any free-phonon normalization or subtraction appears equally in numerator and denominator and cancels.

For bounded complex half-path functions `f,g`, define reflection `Θf` as complex conjugation of `f` evaluated on the reflected second half. The same calculation gives

```
Z ⟨(Θf)g⟩ = ∫ dx dy K_s(x,y)² Tr[Q_f(x,y)†Q_g(x,y)]. (9)
```

The right side is a Hilbert-space inner product of matrix-valued endpoint functions, weighted by the positive heat kernel square. Therefore

```
⟨(Θf)f⟩ ≥ 0,
|⟨(Θf)g⟩|² ≤ ⟨(Θf)f⟩⟨(Θg)g⟩.                   (10)
```

All integrals exist in finite volume. The norm of a path evolution is bounded by `exp[∫(C+Σ_e C_e|b_e(u)|)du]`, with finite constants from the finite-dimensional electron space. The Gaussian bridges and harmonic endpoint integration have finite exponential moments of these linear-growth functionals. One may also dominate `|b|` by `εb²+1/(4ε)` and choose sufficiently small `ε` to obtain a Gaussian quadratic domination. Bounded `f,g` do not change this. Time-sliced products converge under this domination to the evolution of every continuous path, which proves the path formula as well as the inequalities. The fixed-number proof uses exactly the same finite-dimensional matrix construction. In particular it does not assume the misleading conditional canonical positivity discussed in §3.

For the grand-canonical two-spin model the full path weight is a determinant square and hence is a positive measure. With `f=1_E` for an event on one half and `g=1`, (10) gives the concrete rare-event inequality

```
⟨1_E⟩² ≤ ⟨1_E (Θ1_E)⟩.                            (11)
```

Its right side is the **interacting** event weight of the reflected two-half restriction. No free-Gaussian rarity has been substituted for interacting rarity. In a canonical sector equations (9)–(10) remain valid as functional inequalities, but a general path event is not automatically a physical probability.

This gives an invertible reference in the reflected comparison: each mirrored half path has full electronic period matrix `M_w†M_w`, positive definite because finite-time matrix evolutions are invertible. In the number-conserving quadratic representation write `m_w` for the corresponding one-particle evolution, excluding the scalar fugacity; its mirrored one-particle period matrix is `m_w†m_w`. Thus `det(I+q m_w†m_w)` is strictly positive for every real positive `q`. One can use these mirrored references in the Schwarz comparison without the zero-trace defect of unrestricted histories. The mirrored quantity is an upper comparison, not the original history law; obtaining a small error on replacing one by the other is a separate issue. Equation (11) alone contains no factor of `Ω` and proves no linear temperature ceiling.

## 6. Positive route: opposite-spin fugacity averaging

**PROVED.** On a finite lattice with `V` sites, let `Z` be the original grand-canonical partition function and let `Z_bal` be its restriction to `N_↑=N_↓`. For either pure model,

```
Z_bal ≤ Z ≤ (V+1)Z_bal.                               (12)
```

They consequently have the same thermodynamic pressure whenever the limit exists, since `0≤(τ/V)log(Z/Z_bal)≤(τ/V)log(V+1)→0`. This is not yet a fixed-density source or correlation equivalence.

Proof: for a real phonon history and real fugacity `z=exp(βμ)`, put

```
D(z)=det(I+zY)=Σ_(n=0)^V a_n z^n,
```

where all `a_n` are real and `a_0=1`. Opposite-spin twists give the positive integrand

```
W_bal = (1/(2π))∫_0^(2π) |D(ze^(iθ))|²dθ
      = Σ_n a_n² z^(2n) ≥ 1.                          (13)
```

The identity is finite Fourier orthogonality. After integrating the phonon histories, angular averaging is the projection onto equal spin numbers, because both spin numbers are conserved by the full Hamiltonian. Thus its partition function is exactly `Z_bal`, not an invented modified ensemble. Cauchy–Schwarz in `R^(V+1)` gives `D(z)²≤(V+1)Σ_n a_n²z^(2n)` pointwise, proving the upper bound in (12) after Gaussian integration. The lower bound follows from the positivity of the physical quantum traces in all conserved spin-number sectors: `Z` is their sum and `Z_bal` is a subset of that sum. It does not follow from a pointwise inequality between the two history integrands, which need not hold.

The averaged **unnormalized** functional has no zero denominator, so it is one precise way to keep electronic weights meaningful near determinant zeros. However it does not alone repair normalized Wick or pair bounds. For the three-site monodromy (1) at `z=1`, (13) is `W_bal=4`. The unnormalized one-spin number-conserving insertion matrix is

```
D(q) N_q = q(1−q)(Y−qI),                              (14)
```

by (2)–(3), as a polynomial identity including the originally singular values. An off-diagonal insertion at `x≠y` is therefore `q(1−q)Y_xy`. Opposite twists for the two spins give the onsite pair-transfer insertion polynomial and its complex conjugate. At `q=exp(iθ)` its angular-averaged numerator is `2|Y_xy|²`, since the mean of `|1−exp(iθ)|²` is two. Its conditional normalized pair transfer is `|Y_xy|²/2`, which already equals `200/9` for `Y_12=20/3` in (1). Conjugating by the extra boost in §3 makes some coordinate entries grow without bound, and hence makes this averaged conditional pair transfer unbounded too. The full physical balanced state still has ordinary norm bounds after integrating the histories. The conditional angular-averaged object is not itself a positive electronic density matrix.

Thus fugacity averaging achieves strict positivity of an integrated reference weight and exact pressure equivalence, but an integrated **spatial** determinant/minor estimate remains essential. It cannot be replaced by an amplitude-independent pointwise covariance bound. A finite twist selection is likewise not automatically the original ensemble: it needs a proved Fourier-projection or comparison identity, not just an invertible matrix at a selected angle.

## 7. Exact remaining positive lemma targets

The reflection-positive construction (8)–(11) solves the existence of a valid imaginary-time Schwarz comparison for both models and fixed particle number, including the noncommuting SSH terms. What is missing is an **Ω-dependent bound on the reflected interacting restrictions**, or a quantitative comparison of mirrored histories to the original full functional at `βΩ` of order one. One must establish that bound with the harmonic bridge weight and electronic matrix norm both present. Free-Gaussian excursion estimates alone do not supply it. No chessboard iteration or multiscale convergence theorem has been asserted here without proof.

The balanced construction (12)–(14) gives a positive unnormalized reference with the correct grand-canonical pressure. What is missing is control of weighted spatial insertion minors, and transfer of the relevant sourced estimate to a fixed-density state. Nonnormality persists after twist averaging, so a lower bound on the averaged determinant by itself is insufficient.

For the source route, the atomic calculation suggests the precise lattice target: a quadratic bound on the source-dependent **connected** quantum/static interpolation difference, with coefficient of order `Ω` at fixed `βΩ`, uniform in volume and relevant momentum and compatible with the fixed-density limit. The additive free-energy sandwich alone supplies no such coefficient. Even obtaining that estimate would also require showing how it excludes the specified algebraic electronic order, rather than only bounding a uniform-source response with an unsuitable denominator. These steps are OPEN.

No physical lower bound on `T_pair`, false-regime claim, or violation of the linear conjecture follows from this slice. Completed contributions are the unrestricted-reference barrier, the atomic source-response lemma, the imaginary-time reflection-positive comparison, and the balanced-spin pressure equivalence, with an exact demonstration of the latter's remaining nonnormality obstruction.

## 8. Exact interacting atomic path law and loss of projected RP

This section addresses the proposed chessboard bound on rare **nonzero-Matsubara** histories. For one oscillator, the energy-coordinate Euclidean action is

```
S_ph(b) = (1/(2U))∫_0^β [b(u)²+Ω^(−2) b'(u)²]du.
```

Its stationary Gaussian covariance is

```
C_Ω(s) = (UΩ/2) cosh[Ω(β/2−d(s,0))]/sinh(βΩ/2),
```

where `d(s,0)` is the circle distance in `[0,β/2]`. This follows by solving `(-∂_s²+Ω²)C=UΩ²δ_periodic`, with continuity, the derivative jump `−UΩ²` at zero, and periodicity. Fourier inversion equivalently gives coefficients `UΩ²/[β(ω_n²+Ω²)]`, `ω_n=2πn/β`. Let `b_0=β^(−1)∫b` and `v=b−b_0`. The zero coefficient has variance `U/β=Uτ`, is independent of all nonzero coefficients, and

```
D_Ω(s) = E[v(s)v(0)] = C_Ω(s)−U/β.                  (15)
```

This is also the **exact interacting** dynamic law of the isolated Holstein atom, grand canonically or in any number sector. In the coordinate path representation its electron occupation is conserved, and the electron factor in sector `n` is `exp[−(n−1)∫b]`, multiplied by any real fugacity and spin multiplicity. It depends only on `b_0`. Thus integrating or conditioning on the electron occupation and static coordinate changes the law of `b_0` but leaves the independent law of `v` exactly (15). The same factorization holds on an atomic lattice, even at fixed total particle number, since its electronic factor is a function of the static coordinates alone. This atom has `U>0` and finite positive frequency, and no approximation to its interacting weight is made. It is an auxiliary `t=0` model, outside the fixed-`t>0` square-lattice target.

**PROVED loss of RP.** On the projected path space, use reflection `θu=β−u` and the half-path functional `f(v)=v(β/4)`. Equation (15) gives

```
E[(Θf)f] = D_Ω(β/2)
 = Uτ [a/(2sinh(a/2))−1] < 0,       a=βΩ>0.         (16)
```

The inequality follows from `sinh x>x` for every `x>0`, by integrating `cosh x>1`. The unbounded linear test can be replaced by the bounded truncation `f_R=max(−R,min(R,f))`: its reflected product converges in expectation to (16) as `R→∞`, since Gaussian second moments are finite and the truncations converge in `L²`. Consequently a sufficiently large finite truncation is a bounded RP counterexample.

This does not contradict §5. In the full `b`-path space, `v(β/4)=b(β/4)−β^(−1)∫b` uses both halves and is not a half-path functional. Conditioning on `b_0` imposes a global integral constraint and destroys the Markov bridge factorization that proved (9). Therefore neither the full-path RP inequality nor its chessboard iterations can be applied directly to the zero-mode-conditioned Gaussian law.

For orientation, the small-`a` reflected covariance is

```
D_Ω(β/2) = −UΩ²/(24τ)+O(UΩ⁴/τ³),                   (17)
```

by the Taylor expansion of `sinh`. At fixed `a` it is a strictly negative constant times `UΩ`. The obstruction exists at every positive frequency; it is not an order-of-limits defect.

## 9. A legitimate local reflected event still has the same cutoff scale

Take `Δ=b(β/4)−b(0)=v(β/4)−v(0)`. This is a half-path observable in the full path space, so its event indicators are eligible for §5. In the exact interacting atomic law it is a centered Gaussian with variance

```
σ² = 2[C_Ω(0)−C_Ω(β/4)]
   = UΩ H(a),
H(a) = [cosh(a/2)−cosh(a/4)]/sinh(a/2)
     = (1−exp(−a/4))(1−exp(−3a/4))/(1−exp(−a)).       (18)
```

In particular `H(a)=3a/16+O(a³)` as `a→0`. The following bounds are exact:

```
21a/128 ≤ H(a) ≤ 3a/14       for 0<a≤1,
H(a) ≥ 3/35                  for a≥1.                 (19)
```

Proof of the first bounds: the positive Taylor series of the cosh difference gives a numerator at least `3a²/32`. For `0≤x≤1/2`, `(2n)!≥2^n` implies `cosh x≤1/(1−x²/2)≤8/7`. Integrating `sinh x≤x cosh(1/2)` from `a/4` to `a/2` gives a numerator at most `3a²/28`. The denominator lies between `a/2` and `(a/2)cosh(a/2)≤4a/7`. Their ratios give the first two constants in (19). For `a≥1`, use the product form in (18), drop its denominator since it is at most one, and use `1−exp(−x)≥x/(1+x)`. The numerator is at least `(1/5)(3/7)=3/35`.

Define the small-reference event `G_r={sup_u |v(u)|≤r}` and its complement `B_r`. Since `G_r` implies `|Δ|≤2r`, the Gaussian density bound yields

```
P(G_r) ≤ P(|Δ|≤2r) ≤ 4r/(sqrt(2π)σ) ≤ 2r/σ.         (20)
```

The last inequality only uses `π≥2`. Fix a dimensionless `κ>0`, put `r=κτ`, and set `Q=τ³/(UΩ²)`. For `a≤1`, (19)–(20) imply

```
P(B_(κτ)) ≥ 1−2κ sqrt(128/21) sqrt(Q).                (21)
```

For `a≥1` they instead give

```
P(B_(κτ)) ≥ 1−2κ sqrt(35/3) τ/sqrt(UΩ).               (22)
```

The same lower bounds hold for the legitimate local bad event `E_r={|Δ|>2r}`, because (20) bounds its complement directly. Thus if `τ=o((UΩ²)^(1/3))` at fixed `U>0` and `Ω→0`, both `P(B_(κτ))` and `P(E_(κτ))` tend to one. In the `a≤1` branch use `Q→0`; in the `a≥1` branch use `τ≤Ω` in (22).

At the proposed linear scale `τ=Ω/a` with fixed `a>0`, the sharper direct bound is

```
P(E_(κτ)) ≥ 1−[2κ/(a sqrt(H(a)))] sqrt(Ω/U) → 1.     (23)
```

For a finite collection of nontrivial reflected increments on slabs of width proportional to `1/Ω`, at fixed `a`, each increment has variance `UΩ` times a positive fixed constant. The union bound therefore gives `P(all increments bad)≥1−C sqrt(Ω/U)→1`. Independence is not needed. In particular a uniformly positive, Ω-independent suppression `exp(−c×number of slabs)` for these small-cutoff restrictions cannot hold. A proposed extra fixed factor `βΩ=a` in its exponent does not fix the contradiction. If the full period is shorter than `1/Ω`, there is not even a slab of that proposed width to tile; at the cutoff threshold in §10, `βΩ→0`.

Under a valid RP/chessboard comparison `P(E)≤P(E_repeated)^(1/m)`, this lower bound forces the repeated-restriction root to tend to one too. RP cannot turn a genuinely common event in the interacting atomic model into a uniformly rare one.

## 10. Matching upper bound: sharp two-thirds-power cutoff threshold

For the same exact interacting atomic `v` law, there is a matching sufficient tail bound at the same Ω power:

```
P(B_r) ≤ 3 exp[−τ r²/(4UΩ²)],
P(B_(κτ)) ≤ 3 exp[−κ² Q/4].                           (24)
```

Here is a self-contained Gaussian proof. Let `v_*` be the mean-zero periodic massless Gaussian with Fourier variances `UΩ²/(βω_n²)`, `n≠0`. Its covariance dominates (15), since

```
UΩ²/(βω_n²) − UΩ²/[β(ω_n²+Ω²)] ≥ 0.
```

One can couple `v_*=v+z` with `z` an independent centered Gaussian having these difference coefficients. The Fourier standard deviations of `z` are summable (of order `n^(−2)`), so its Fourier series is almost surely uniformly convergent and defines a continuous process. For the convex functional `F(w)=exp(λ||w||_∞²)`, conditional Jensen gives `E F(v)≤E F(v_*)`.

Let `B(u)` be a standard Brownian bridge on `[0,β]`, with covariance `min(u,w)−uw/β` and endpoints zero. Then

```
v_* = sqrt(U)Ω [B−β^(−1)∫B],
||v_*||_∞ ≤ 2sqrt(U)Ω ||B||_∞.
```

To verify its covariance, double-centering the Brownian bridge covariance gives `β/12−d/2+d²/(2β)` at circle distance `d∈[0,β]`. Its Fourier coefficients are `1/(βω_n²)` for `n≠0`, as obtained by twice differentiating that quadratic away from the derivative jump at zero. This proves the representation.

The one-sided bridge maximum has tail `P(sup B≥x)=exp(−2x²/β)`. Indeed reflect a Brownian-motion trajectory at its first passage of level `x`: paths ending at zero after crossing `x` have the transition density of motion ending at `2x`. Dividing its density at time `β` by the density at zero, and conditioning the endpoint to zero by a shrinking interval, gives `exp[−(2x)²/(2β)]`. Symmetry and a union bound then give `P(||B||_∞≥x)≤2exp(−2x²/β)`. Consequently, with `A=UΩ²β`,

```
P(||v_*||_∞²>t) ≤ 2exp[−t/(2A)].
```

Integrating this tail at `λ=1/(4A)` gives

```
E exp[||v_*||_∞²/(4A)]
 = 1+∫_0^∞ (1/(4A)) exp[t/(4A)] P(||v_*||_∞²>t)dt
 ≤ 1+2(1/(4A))/(1/(2A)−1/(4A)) = 3.
```

Conditional Jensen transfers the same moment bound to `v`; Markov's inequality proves (24). This estimate is for the exact atomic interacting measure, not an assumed replacement of a general lattice interacting measure by a Gaussian.

For every fixed target probability `0<p_0<1`, let `τ_min(Ω)` denote the infimum of temperatures where `P(B_(κτ))≤p_0`. Equations (21)–(24) show

```
τ_min(Ω) = Θ_(κ,p_0,U)((UΩ²)^(1/3)) as Ω→0.         (25)
```

More explicitly, in the `a≤1` regime a necessary condition is

```
Q ≥ 21(1−p_0)²/(512κ²),
```

and the sufficient condition is `Q≥(4/κ²)log(3/p_0)`. The `a≥1` regime supplies no lower-temperature exceptions, by (22). These two bounds establish the exact **power** in (25), not an optimal numerical quantile constant.

## 11. Precise method-class verdict and remaining scope

The sharp cutoff barrier concerns positive rare-restriction methods with all of these properties:

1. Good quadratic references obey the fixed dimensionless perturbation bound `β sup|v|≤κ`, so their complement is `B_(κτ)` (or a local increment restriction necessary for that bound).
2. The returned bad-history budget is a positive probability/restriction comparison, or a domination of it, without compensating signed connected insertions.
3. The restriction estimate is claimed uniformly over the auxiliary quadratic references, including the atomic Holstein reference. This is a method-uniformity premise, not a claim that the fixed-`t>0` target is atomic.

In that class a small Ω-independent budget cannot hold at `τ=O(Ω)` or `τ=Ω` times any fixed logarithmic power. The exponent accessible from the cutoff probability is precisely `2/3`, by (25). RP/chessboarding cannot improve that exponent while the three premises remain in force. Furthermore the zero-mode-conditioned measure cannot be used as an RP input at all, by (16). These are two distinct failures: invalid reflection locality after zero-mode removal, and an unchanged fluctuation scale even for valid full-path local increment events.

The conclusion does **not** cover mirrored references permitting natural amplitude `sqrt(UΩ)` at `βΩ` of order one. Such a reference can remain invertible by §5 and does not satisfy premise 1; its temporal/spatial insertion estimates still require proof. Nor does the conclusion exclude integrated weighted cancellation, a refreshed reference with a non-pointwise transport estimate, or a connected lattice source estimate. A natural-amplitude cutoff can have a nontrivial Gaussian budget at the linear scale, but that observation alone supplies no forest summability theorem.

No square-lattice ordering-temperature inequality has been obtained by combining these restrictions with fugacity averaging or LINEAR4. LINEAR4's interacting tail-replacement review is not assumed resolved here. The full linear or Ω-times-log ceiling at fixed positive regular density and fixed `t>0` remains OPEN; the completed new result is the sharp obstruction to the specifically stated small-reference positive-RP/chessboard route.

## 12. Positive result: kinetic Gaussian domination with the interacting weight

**PROVED, analytical; independent review pending.** Consider either finite-volume pure model, with the real coordinate representation

```
H = (D/2)Σ_(e=1)^M p_e² + V(b),
D=UΩ²,    p_e=−i∂_(b_e),
V(b)=Σ_e b_e²/(2U) I + H_el(b).
```

The electron space can be full Fock space or any invariant fixed-number sector. The matrix potential is Hermitian, is quadratic scalar plus a linear matrix-valued electronic coupling, and commutes with multiplication by any scalar function of the coordinates. A real chemical potential is included in `H_el`. The constants are uniform in finite volume, hopping, and all external static parameters of that potential. This does **not** assert an estimate conditional on a prescribed zero Matsubara coordinate: such conditioning is a nonlocal-in-time constraint and is not the trace of the Hamiltonian used below.

Take times partitioning the thermal circle into positive lengths `δ_j`, `Σ_jδ_j=β`. Insert scalar multipliers `exp(λ_j·b)` at the vertices, with `Σ_jλ_j=0` coordinatewise. Choose constant vectors `q_j` on each interval such that the vertex coefficient is `λ_j=q_(j−1)−q_j`. Let `q̄=β^(−1)Σ_jδ_j q_j`. The unnormalized Euclidean insertion trace obeys

```
|Z E exp(Σ_j λ_j·b(t_j))|
 ≤ Z exp[(D/2)Σ_jδ_j |q_j−q̄|²].                    (26)
```

In the grand-canonical two-spin model this is an ordinary moment-generating-function bound for the **full interacting positive path measure**. The same is true in each fixed balanced sector `N_up=N_down=k`, as proved in §15. It makes no conditioning on the zero mode and no replacement of the interacting measure by a free Gaussian. In general fixed-total-number sectors (26) is only a trace-functional bound; positive path-event probabilities are not asserted there.

### 12.1 A Schatten estimate without positivity of conditioned electronic matrices

We first prove the needed estimate. For a real constant vector `q`, define the self-adjoint first-order operator `A_q=Dq·p`. For every `t>0` and `1≤p<∞`,

```
||exp[−t(H+iA_q)]||_(S_p) ≤ [Tr exp(−ptH)]^(1/p).    (27)
```

Finite-dimensional proof: Lie product expansion gives the limit of `[exp(−tH/n)exp(−itA_q/n)]^n`. The second factor is unitary. Schatten Hölder with exponent `np` for each of the n heat factors and operator norm for the unitaries bounds the p-norm by

```
Π_(k=1)^n [Tr exp(−ptH)]^(1/(np))
 =[Tr exp(−ptH)]^(1/p).
```

Continuity of matrix exponentiation proves the finite-dimensional statement.

For the present unbounded operators the same bound follows by spectral Galerkin approximation, whose details matter. The confining quadratic-plus-linear potential gives `H≥(D/2)p²+|b|²/(4U)−C`, with finite C. Thus H has compact resolvent and finite heat trace, and the form of `A_q` is infinitesimally bounded relative to the shifted H form: Cauchy–Schwarz and `2xy≤εx²+ε^(−1)y²` bound `|⟨ψ,A_qψ⟩|` by `ε⟨ψ,(H+C)ψ⟩+C_(ε,q)||ψ||²`. Its off-diagonal form is bounded on the same form domain. The closed sectorial form `h+i a_q` consequently defines `H+iA_q`.

Let `P_L` project onto the first L eigenvectors of H. Their span is dense in the H form domain. The finite Galerkin resolvent equations for `P_LHP_L+iP_LA_qP_L+z` converge to the full form equation: after shifting the real part of z, coercivity is supplied by h, and the form-norm error is bounded by a fixed constant times the best approximation error from `P_L`. This proves strong resolvent convergence; the same coercive sectorial bounds give strong semigroup convergence for every positive time, for instance by the resolvent contour integral with a common sector and dominated resolvent bounds. The finite-dimensional right side of (27) is the truncated heat trace and is at most the full heat trace. Schatten norms are lower semicontinuous under this strong convergence, by restricting to arbitrary finite orthonormal systems and then taking their supremum. Hence (27) holds for the full operator. This avoids applying an unproved non-Hermitian thermal-density estimate to a conditional electron history.

### 12.2 Coordinate gauge and trace comparison

On a smooth coordinate core,

```
exp(q·b) H exp(−q·b)
 = H+iDq·p−D|q|²/2.                                 (28)
```

All potential terms commute with the scalar multiplier; only `p` changes to `p+iq`. No electronic bond commutation is used.

The corresponding heat-kernel identity is valid without treating the unbounded multiplication operator as a bounded similarity. The semigroup of `H+iDq·p` has diffusion with constant drift `−Dq`. Conditional on its endpoints its coordinate bridge is the same bridge as for zero drift. The Gaussian endpoint density differs by `exp[q·(x−y)−D|q|²t/2]`, while the matrix-valued potential insertion is unchanged on that bridge. Therefore

```
exp[q·(x−y)] K_H(t;x,y)
 = exp[D|q|²t/2] K_(H+iDq·p)(t;x,y).                 (29)
```

One can first use bounded matrix potentials and time slicing. The original harmonic-plus-linear potentials follow by truncation and Gaussian domination of their linear-growth terms, exactly as for the finite-volume path integrals in §5. This establishes (29) as a kernel identity and avoids assuming trace invariance under an arbitrary unbounded similarity.

In the cyclic product of kernels the factors `exp[q_j·(x_(j+1)−x_j)]` give exactly the vertex insertions with coefficient `q_(j−1)−q_j`. Subtracting the same `q̄` from every interval leaves those vertex coefficients unchanged. Equation (29) factors out `exp[(D/2)Σ_jδ_j|q_j−q̄|²]`. Apply trace Hölder to the remaining semigroups with exponents `p_j=β/δ_j`, and then (27):

```
|Tr Π_j exp[−δ_j(H+iA_(q_j−q̄))]|
 ≤ Π_j [Tr exp(−βH)]^(δ_j/β)=Z.
```

This proves (26), including the interacting electronic trace and the phonon weight. Its constants do not involve volume, electronic matrix dimension, hopping, or an electron-conditioned Gaussian variance.

### 12.3 Explicit consequences

For a single coordinate increment over `s∈(0,β)`, choose `q=η` on that arc and zero elsewhere. Equation (26) gives

```
E exp[η(b_e(s)−b_e(0))]
 ≤ exp[Dη²s(β−s)/(2β)].                              (30)
```

For `v_e(t)=b_e(t)−β^(−1)∫b_e`, approximate the uniform integral source by finite Riemann sums. The limiting mean-zero gauge is `q(s)=η(s/β−1/2)`, with squared integral `η²β/12`. Positive-path Fatou, continuity of the paths, and (26) yield

```
E exp(ηv_e(t))≤exp(Dη²β/24),
P(|v_e(t)|≥r)≤2exp[−6r²/(Dβ)].                       (31)
```

This is a positive full-interacting Gaussian domination statement, not projected reflection positivity. In particular the negative projected RP form (16) does not obstruct it.

For increments indexed by j, suppose their time arcs are disjoint whenever they involve the same oscillator (different oscillator coordinates may have overlapping arcs). If their lengths are `L_j` and their thresholds are `r_j`, the source gauges have orthogonal/disjoint support, so subtraction of the common mean can only reduce their squared norm. Applying (26) to each sign pattern, optimizing `η_j=r_j/(DL_j)`, and summing over the `2^m` sign patterns gives

```
P(∀j: |Δ_j|≥r_j)
 ≤ Π_j [2exp(−r_j²/(2DL_j))].                        (32)
```

There is no independence assumption. Equation (32) is precisely a factorized bound on joint interacting local restrictions.

## 13. Natural-amplitude slab excursions have a linear-Ω budget

**PROVED.** The increment comparison also controls whole continuous slabs, with explicit constants. Let `I_j` be time slabs of lengths `L_j`, disjoint for a given oscillator, and define

```
F_j={sup_(t∈I_j)|b_(e_j)(t)−b_(e_j)(left(I_j))|>R_j}.
```

If `x_j=R_j²/(32DL_j)≥16` for every j, then in the full grand-canonical interacting measure

```
P(∩_j F_j)≤Π_j [4exp(−x_j)].                         (33)
```

Proof: use the dyadic grids of a slab. At level n the relevant parent-child increments have length at most `L_j/2^n` and number at most `2^n`, including the full endpoint increment at level zero. Set `a_n=(1/4)(3/4)^n`, so `Σ_n a_n=1`. Continuity and telescoping imply that a slab excursion above `R_j` forces at least one level-n increment above `a_nR_j`. Across distinct slabs, any selected increment arcs are disjoint in the sense needed for (32). A countable union bound and (32) therefore give

```
P(∩_j F_j)
 ≤ Π_j [2Σ_(n≥0) 2^n exp(−x_j(9/8)^n)].             (34)
```

Consecutive terms in each sum have ratio at most `2exp(−x_j/8)≤2exp(−2)<1/2` when `x_j≥16`. The last inequality already follows from `exp(2)>4`. Thus each bracket is at most `4exp(−x_j)`, proving (33).

In particular choose slabs of width `L=ℓ/Ω` and natural amplitude `R=ρsqrt(UΩ)`. Then

```
x=ρ²/(32ℓ),
P(∩_(j=1)^m F_j)
 ≤ exp[−m(ρ²/(32ℓ)−log 4)]
```

provided `ρ²≥512ℓ`. This is an Ω-independent cost per slab, so a tiling has a cost proportional to its actual number of slabs, of order `βΩ/ℓ` per oscillator. No extra factor of `βΩ` per individual slab is asserted. The bound includes the interacting fermion weight, holds uniformly in volume for both models, and permits simultaneous slabs on different oscillators.

For comparison, using the small cutoff `R=κτ` gives `x=κ²τ²/(32UΩℓ)` on these slabs, which tends to zero at fixed `βΩ` as Ω tends to zero. The natural-amplitude positive result therefore escapes the small-cutoff sharpness proposition by abandoning its first hypothesis; it does not contradict SHARPNESS.md.

## 14. Why this positive budget does not yet close the forest at the linear scale

The interacting restriction bound is now proved at the natural amplitude. However simply combining it with a reference-uniform, all-insertion-time Wick forest is still invalid. There is an exact counterexample even among continuous, palindromic, invertible mirrored SSH references.

Take a single two-site bond with one-particle matrix `K=[0 1;1 0]`, static coordinate `b_0=t`, and Hamiltonian `h(u)=−v(u)K`. Use the four quarter profiles `f_j` of §2 and

```
v(u)=η[f_1(u)−f_2(u)−f_3(u)+f_4(u)].
```

This is continuous and palindromic under `u→β−u`, has zero time average, and has peak amplitude `r=8η/β`. All its one-particle generators commute. The full period matrix is therefore `I`, with strictly positive normalized electronic trace, and its midpoint and endpoint thermal occupations are `I/2`. At the end of the first quarter the evolution is `exp(ηK)`. Its unequal-time Wick contraction from zero to that quarter is

```
G(β/4,0)=exp(ηK)/2,
||G(β/4,0)||=exp(η)/2.                               (35)
```

The formula follows either from the exterior trace with full period matrix I or directly in the eigenmodes of K. A unit positive-eigenmode creation/annihilation pair realizes the indicated `1×1` contraction. For `r=ρsqrt(UΩ)` and `β=a/Ω`,

```
η=(ρa/8)sqrt(U/Ω)→∞.
```

Thus invertibility, palindromicity, and the positive restriction bound do not provide a uniform all-time Wick-minor bound at the natural amplitude. The example embeds in square SSH by canceling all other static hoppings. It is an auxiliary-history counterexample, not a physical pairing-temperature counterexample. It also does not disprove estimates retaining the full unnormalized insertion weight or reflection-paired operator norms.

Accordingly priority 1 has a positive natural-amplitude interacting slab budget, but priority 2 requires a new temporal/spatial insertion or transport estimate, not a direct substitution into LINEAR4's small-reference forest. The full linear or Ω-times-log lattice ceiling remains OPEN. The immediate positive target is to exploit the reflected Hilbert–Schmidt comparison for charged insertions in a way that avoids the arbitrary-time contraction (35), or to prove a connected-source interpolation estimate directly.

## 15. Which fixed-density statement has actually transferred

**PROVED.** At zero pair source, (30)–(34) hold in every canonical sector with `N_up=N_down=k`, uniformly in k and volume. For a real number-conserving history its one-spin canonical trace is the real coefficient `a_k` of `det(1+zM)`. Its two-spin balanced trace is `a_k²≥0`. The normalized phonon/history integral is therefore a positive measure. The scalar coordinate gauges preserve both spin numbers, and the proof of (26) applies to their finite-dimensional electron sector. All moment-generating-function and probability arguments consequently apply with that sector's exact interacting partition function. For any target `0<n<2`, taking `2k/V→n` gives the same volume-uniform slab bound at that density without a grand-canonical change of ensemble.

There is also an exact zero-source **pressure** comparison with the unrestricted fixed-total-number sector when N is even:

```
Z_(N/2,N/2) ≤ Z_N ≤ (N+1) Z_(N/2,N/2).                (36)
```

Proof: the pure spinful Hamiltonian commutes with global SU(2) spin rotations, including its phonon coupling. In the sector with total number N, decompose the electron representation into integer-spin irreducibles `V_S` and multiplicity spaces; append the phonon space to each multiplicity space. The Hamiltonian is the identity on `V_S` tensor the corresponding multiplicity-space Hamiltonian. Every integer-spin irreducible has exactly one `S_z=0` vector, while its full spin dimension is `2S+1≤N+1`. Write the positive heat trace on the multiplicity space as `w_S`. Then `Z_(N/2,N/2)=Σ_S w_S` and `Z_N=Σ_S(2S+1)w_S`, proving (36). The free-energy difference per site is at most `τ log(N+1)/V` and vanishes along bounded-density thermodynamic sequences.

This is not an inequality between arbitrary event-restricted canonical path integrals, whose fixed-total-number history coefficients can have either sign. It also does not automatically transfer a decay estimate or an electronic pair-source curvature: equality of zero-source free energies alone supplies neither one. We have proved a fixed-density **balanced-sector restriction budget** and the zero-source fixed-density pressure comparison. A completed pairing-order forest plus its appropriate source/ensemble transfer is still required for a full fixed-total-density pairing ceiling. Odd total numbers and other density sequences have no new pair-correlation claim here.

## 16. Testing adiabatic multiscale: instantaneous equilibrium tracking fails uniformly

**PROVED precise barrier, not a no-go for every multiscale construction.** The palindromic commuting SSH reference of §14 also refutes an instantaneous-equilibrium approximation which replaces its time-dependent conditional occupation by the Fermi function of the instantaneous Hamiltonian at temperature τ, uniformly in external static coordinates.

Its period matrix is I, and its evolution matrices commute with one another. Therefore its conditional equal-time one-particle occupation is exactly `n_hist(u)=I/2` at every time: transporting `(I+M)^(−1)M=I/2` by any prefix evolution leaves it unchanged. At the peak of the first quarter `v=r`, the instantaneous Hamiltonian is `h_inst=−rK`. Its instantaneous Fermi occupation is `n_inst=(I+exp(βh_inst))^(−1)`, with eigenvalues `(1±tanh(βr/2))/2`. Hence

```
||n_hist−n_inst|| = (1/2)tanh(βr/2) → 1/2             (37)
```

for `r=ρsqrt(UΩ)`, `β=a/Ω`.

The field really is slowly varying, not a rapid-pulse trick. On a quarter its form is `v=±r sin²(4π(u−u_j)/β)`, so `|v'|≤4πr/β`. For any window of length δ,

```
|Δv|≤4πrδ/β.
```

At a fixed electronic-band window `δ=1/t`, this is at most `4πρsqrt(U)Ω^(3/2)/(at)`, small against fixed t; it is also small against r. Nevertheless (37) stays order one. The external static coordinate `b_0=t` cancels the bare hopping, and the instantaneous Hamiltonian crosses a zero gap when v changes sign. Fixed bare t does not provide a uniform instantaneous gap or a uniform separation from the Fermi level in the annealed field ensemble.

Thus one cannot justify local thermalization at τ from small window increments measured against the bare band scale alone. A valid tracking scheme must retain the global imaginary-time occupation/memory, treat zero-gap crossings, or prove an annealed exception estimate rather than assuming uniform static-field separation. This example does not rule out such a scheme, and no linear-scale lattice conclusion follows from this barrier.

## 17. Testing chessboard on positive balanced sectors: ordinary spatial RP fails

**PROVED.** Positivity of the balanced path weights is not spatial reflection positivity. In fact ordinary spatial RP of the Holstein phonon-coordinate functional fails at every finite β and Ω in every nonempty, nonfull balanced canonical sector, on reflection-symmetric finite lattices divided into reflected halves A and B. This is distinct from the valid imaginary-time RP of §5 and from the negative projected imaginary-time covariance of §8.

Let `b_(0,i)=β^(−1)∫_0^β b_i(u)du` be the static Matsubara coordinate. Use the density coupling `b_i n_i` and the balanced sector `N_up=N_down=k`, with `0<k<V`, so `N=2k` is fixed. Denote the Duhamel covariance of bounded electronic observables by

```
Cov_D(X,Y)=β^(−1)∫_0^β ⟨X(s)Y⟩ds−⟨X⟩⟨Y⟩.
```

The following exact zero-mode source identity includes the full interacting phonon trace:

```
E exp(Σ_i ξ_i b_(0,i))
 = exp[(Uτ/2)Σ_i ξ_i²]
   Z[H+UτΣ_i ξ_i n_i]/Z[H].                         (38)
```

Proof: the free harmonic zero coordinates have independent Gaussian densities proportional to `exp[−βΣ_i b_(0,i)²/(2U)]`, independent of the nonzero modes. Complete the Gaussian square with the source, then translate each zero coordinate by `Uτξ_i`. The electronic history acquires exactly the constant potential `UτΣ_i ξ_i n_i`. This change of variables is valid before and after taking the balanced trace and includes its exact determinant weight. Differentiating (38) at zero gives

```
E b_(0,i)=−U⟨n_i⟩,
Cov(b_(0,i),b_(0,j))
 =Uτ δ_(ij)+U² Cov_D(n_i,n_j).                      (39)
```

The second derivative of `log Z(H+εX)` is `β² Cov_D(X,X)`; this follows directly by differentiating the exponential under the trace and using cyclicity. The cross version follows by polarization. No free-Gaussian conditional approximation is used.

Set `n_A=Σ_(i∈A)n_i`, `n_B=N−n_A`, `B_A=Σ_(i∈A)b_(0,i)`, and `B_B` analogously. Reflection symmetry gives `⟨n_A⟩=⟨n_B⟩=N/2`. Because A and B are disjoint, (39) gives

```
E[(B_A+UN/2)(B_B+UN/2)]
 =U² Cov_D(n_A,n_B)
 =−U² Var_D(n_A)<0.                                 (40)
```

For completeness strict positivity of `Var_D(n_A)` is not inferred from a classical occupation distribution. Center n_A and expand in an eigenbasis of the full Hamiltonian. Each term is its squared matrix element times the strictly positive logarithmic mean of the corresponding Gibbs weights. The Gibbs density is faithful at finite β, and n_A is not constant in this fixed sector, so at least one term is positive. Harmonic confinement supplies the heat trace; n_A is bounded, justifying the expansion and the Duhamel integral.

The observable `F=B_A+UN/2` is supported in the spatial half A, and its spatial reflection is `B_B+UN/2`. Equation (40) contradicts spatial RP. Although F is unbounded, its second moment is finite by the confining Gaussian integral; bounded truncations of the same local function converge in L² and retain a negative reflected form for sufficiently large truncation. This supplies bounded test functions as well.

As an exact check, the two-site Holstein atom at t=0 with one electron of each spin has two doublon configurations of weight `exp(2βU)` each and two separated configurations of weight `exp(βU)` each. There `Var_D(n_A)=exp(βU)/(1+exp(βU))`, so the reflected phonon form is `−U²exp(βU)/(1+exp(βU))`, for every Ω. This check is not needed for the arbitrary-t proof (40).

Consequently the proposed spatial chessboard route cannot start merely from positivity of balanced-sector weights: its ordinary spatial RP premise is false for Holstein at fixed density, including fixed t>0. This does not rule out a modified reflection, a grand-canonical spatial comparison, integrated Gaussian domination without spatial RP, or the already proved imaginary-time restriction comparison. No assertion about phonon-only spatial RP in canonical SSH is made. The surviving positive tool is (26)–(34), while a pairing-sensitive spatial estimate remains OPEN.

The new proofs use only finite-dimensional spin algebra, coordinate kinetic energy, Schatten Hölder, and the parent finite-volume path representation. No new third-party paper content or references have been introduced. The analytic semigroup approximation, continuous-time chaining, forest obstruction, spin decomposition and canonical spatial-RP counterexample are not machine-certified by the finite replay.

## 18. Executed certificate and retained mutants

The prescribed laptop-only offline command

```
UV_OFFLINE=1 UV_CACHE_DIR=D:/dev/.cache/uv uv run --no-project --with numpy --with scipy --with mpmath python -B replay.py --all-mutants
```

was run from this folder after the new replay edits. Its observed exit was **0**, and its reported wall time was **6.68 seconds**. All twenty-four child commands `python -B replay.py --mutant NAME` exited **1**, including the ten new modes checking coordinate-gauge sign, optimal centering, dyadic allocation, natural-amplitude power, mirrored unequal-time contraction, drift unitarity, fixed balanced positivity, SU(2) multiplicity, instantaneous occupation and the negative canonical spatial form. RESULT §9 lists their exact modes and rejection guards. The first fourteen retain the full SHARPNESS.md checks.

CERTIFIED covers the finite exact-rational algebra, seven earlier covariance interval points, and one 60-digit finite-matrix Schatten `S_2` spot check. It does not certify the full analytic semigroup proof, continuous-time probability argument, canonical spin decomposition, spatial-RP theorem, or any lattice pairing ceiling. Those are analytical claims with the proofs above and independent review pending. No Lean, TeX, repository tests, lint or broader suites were run. Offline laptop compute only; spend $0; no Git or outward publication.
