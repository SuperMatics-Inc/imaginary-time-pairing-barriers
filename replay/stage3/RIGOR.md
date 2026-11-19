# Stage 3: averaged fermions, slow-field convexity, and charged block transport

## 1. Scope and conventions

**PROVED analytical lemmas and method barriers; independent review pending.** No pairing-temperature ceiling is proved here. All calculations concern the pure two-spin Holstein or bond-SSH models in the parent conventions: `b=FX`, `U=F²/K>0`, `Ω=ℏsqrt(K/m)>0`, `τ=k_B T>0`, `β=1/τ`, with oscillator kinetic coefficient `D=UΩ²`. SSH hopping is `−t+b_e`; Holstein couples `b_i(n_i−1)` when the centered convention is convenient. That centering is just the oscillator translation/chemical-potential change permitted in the parent files. Energies, including τ, U and Ω, have the same units. Imaginary time has inverse-energy units. A field amplitude has energy units.

The full finite-volume Hamiltonian H is self-adjoint, bounded below, and has a finite heat trace because the scalar oscillator potential is quadratic and the electronic matrix coupling is linear and finite dimensional. The grand-canonical heat trace includes a real chemical potential in H. Fixed balanced sectors are specified separately. No theorem about an arbitrary signed canonical history measure is used.

Only the already established parent coordinate representation and elementary finite-dimensional spectral/CAR algebra are used. No third-party paper content, including Pedra–Salmhofer, is invoked. CERTIFIED below concerns the replay's finite rational and interval checks, not a machine verification of these analytic proofs.

## 2. Averaging histories does not preserve Wick's rule

### Proposition 1: an exact interacting atomic cumulant

Consider a single centered Holstein atom, without hopping, grand canonically at centered chemical potential zero:

```
H = (UΩ²/2)p² + b²/(2U) + b(n_up+n_down−1).
```
For electronic occupancy n=0,1,2, translating b by `U(n−1)` gives the same oscillator Hamiltonian minus `U(n−1)²/2`. The translation is unitary on the oscillator space. After dividing by its common oscillator partition function, the four electronic weights are

```
w_empty=w_double=exp(βU/2),   w_up=w_down=1,
Z_el=2(exp(βU/2)+1).
```
Consequently, at every Ω>0,

```
<n_up>=<n_down>=1/2,
<n_up n_down>=exp(βU/2)/[2(exp(βU/2)+1)],
<n_up n_down>−<n_up><n_down>=(1/4)tanh(βU/4)>0.       (S1)
```
The same computation holds for the normalized annealed static integral Ω=0 by completing its Gaussian square. Thus this obstruction is already present in the static limit.

The averaged one-electron density matrix is diagonal in spin, with both occupations 1/2. Off-diagonal spin contractions vanish by separate spin-number conservation. Applying number-conserving Wick's rule to that averaged matrix would give `<n_up n_down>=1/4`, or connected cumulant zero, contradicting (S1). In particular, at τ=CΩ with U fixed, the cumulant tends to 1/4, not zero.

Here is the path-integral interpretation without an interchange-of-averages assumption. At a fixed atomic history the one-period electronic evolution depends only on `b_0=β^−1∫b`: the occupation operators commute at all times. Each spin's grand-canonical conditional state is a Gaussian state with occupation f(b_0). Both spins see the same b_0. The full positive interacting history weight averages `f(b_0)²`, not `(E f(b_0))²`. Their difference is the variance of f under that weight. Averaging the propagator is not averaging its determinant. The field does not enter the interacting four-point function only through an averaged one-electron propagator.

This is a **method barrier**, not an atomic or lattice transition claim. The isolated atom has no spatial pair transport. Static mixtures can be non-Gaussian while still satisfying the parent no-ODLRO theorem.

### Canonical caveat

A fixed-number free-electron reference need not be Gaussian even before histories are averaged. For one spin on two zero-energy orbitals in the N=1 sector, the density matrix is uniform on the two singly occupied states: `<n_1>=<n_2>=1/2`, but `<n_1 n_2>=0`, whereas Wick applied to the averaged one-particle matrix gives 1/4. Balanced positive canonical weights in parent §15 therefore do not themselves license a Gaussian determinant identity. One can instead use a grand-canonical Gaussian reference followed by a legitimate source/ensemble argument, or retain the genuine canonical moments.

## 3. What can be bounded after integration: chronological moments and determinants

### Proposition 2: full grand-canonical ordered moments have norm bound one

Let A_j be bounded electron operators inserted in chronological order, with nonnegative heat intervals δ_j summing to β. Schatten Hölder gives

```
|Tr[e^(−δ_1 H) A_1 ... e^(−δ_m H) A_m]|
 ≤ (Tr e^(−βH)) Π_j ||A_j||.                        (S2)
```
Indeed assign exponent `p_j=β/δ_j` to each nonzero heat interval and operator norm to each A_j. The reciprocals sum to one, and `||e^(−δ_jH)||_(S_pj)=Z^(δ_j/β)`. Zero intervals are limits or are merged with the neighboring bounded insertions. Multiplying yields (S2). It remains valid with the unbounded oscillator space because its heat factors are in the asserted Schatten classes and the electron insertions are bounded.

For normalized single-orbital creation/annihilation operators the norms are one. Every full grand-canonical chronological electron moment is therefore at most one in modulus, after the phonons have been integrated. This does not supply spatial decay or Wick factorization. It is not asserted for charge-changing insertions in a fixed-N heat trace: those insertions leave that sector, so its partition function cannot replace Z in the same Hölder argument.

### Proposition 3: the averaged grand-canonical two-point kernel has a chronological determinant bound one

Define for 0<s<β

```
G_ab(s)=Z^−1 Tr[e^(−βH) c_a(s)c_b†],
c_a(s)=e^(sH)c_a e^(−sH).
```
Extend to the time-ordered, antiperiodic fermionic kernel, using the fixed equal-time convention `G(0+)=<c c†>`. Its opposite boundary value is `G(0−)=−<c† c>` with appropriate matrix indices. A claimed bound on a determinant of this kernel concerns that determinant, not the actual interacting higher-point function.

To prove the bound, take a full-H eigenbasis. For each pair of states m,n for which c removes one electron, put `ω_mn=E_n−E_m`, and define the positive matrix-valued atomic measure

```
dμ_ab(ω)=Z^−1 Σ_mn (e^(−βE_m)+e^(−βE_n))
          <m|c_a|n> overline(<m|c_b|n>) δ_(ω_mn)(dω).
```
Every atom is a nonnegative multiple of a rank-one Gram matrix, and CAR gives

```
∫dμ_ab=<c_a c_b†>+<c_b† c_a>=δ_ab.                 (S3)
```
The Lehmann expansion directly yields

```
G_ab(s)=∫ e^(−sω)/(1+e^(−βω)) dμ_ab(ω), 0<s<β.    (S4)
```
The thermal factor in (S4) is between zero and one for every real ω and 0≤s≤β: for ω≥0 this is immediate; for ω<0 it is `e^(−(β−s)|ω|)/(1+e^(−β|ω|))`.

Construct an auxiliary one-particle Hilbert space with one basis vector per transition (m,n), diagonal one-particle energy ω_mn, and field vectors whose components are the square roots of the weights in dμ times the corresponding c matrix elements. Equation (S3) makes the map from the original orbital vectors an isometry. On any finite subset of transitions it is a contraction, so each truncated auxiliary electron field has norm at most one. Its finite-dimensional free thermal covariance is exactly the truncated sum (S4), with the same negative-time branch. Free CAR Wick's identity says that a chronological product of auxiliary fields is the corresponding covariance determinant, up to the fixed ordering sign. Applying (S2) to this finite free auxiliary model bounds every balanced determinant by one. The finite-transition covariances converge to (S4), also at the chosen one-sided equal-time values, because the thermal factors are bounded by one and μ has finite total mass. Each finite-order determinant is continuous in its entries. Taking the limit proves

```
|det[G_(a_i b_j)(t_i−s_j)]_(i,j=1)^k|≤1.           (S5)
```
No finite heat trace of the infinite auxiliary Fock Hamiltonian is needed: only its finite truncations and the convergent covariance entries were used.

The statement survives positive-semidefinite diagonal-one vertex interpolation. Represent its factors s_ij as inner products of unit vectors ξ_i and replace each auxiliary field vector belonging to vertex i by its tensor product with ξ_i. The field norms remain at most one, and the covariance gets precisely the factor s_ij. The same chronological argument proves (S5) for the interpolated covariance.

This is a **chronological CAR determinant bound**, not a claim of a single time-independent Gram factorization with bounded Gram constant for the discontinuous chronological kernel at arbitrarily many time points. In particular it does not invoke a determinant theorem from the literature. The determinant is a moment in an auxiliary *Gaussian* state. The physical interacting moment need not equal it, by Proposition 1. Its connected higher-point vertices must be retained in any forest expansion. Ignoring them is already wrong at four fermion fields. These vertices are not Ω-small locally, by (S1). That fact by itself does not prohibit a pairing-sensitive spatial cancellation estimate.

## 4. An absolute integrated-field moment is not a small linear-scale Gram budget

### Proposition 4: an explicit exponential-moment barrier

For a free harmonic coordinate let `v=b−b_0` denote its zero-mean dynamic part. Its nonzero Fourier covariance coefficients at `ν_n=2πn/β` are

```
D_n=UΩ²/[β(Ω²+ν_n²)], n≠0.
```
The real linear functional `L=∫_0^β cos(ν_1u)v(u)du` is centered Gaussian with

```
Var L=βUΩ²/[2(Ω²+ν_1²)].                            (S6)
```
To verify the factor 1/2, only n=±1 survive in the double integral against the two cosines. Each cosine coefficient contributes β/2, giving `2 D_1(β/2)²`. For any real c>0, `∫|v|≥L`, so its Gaussian moment gives

```
log E exp[c∫_0^β |v(u)|du]
 ≥ c² Var L/2
 = c² Uβ/[4(1+(2π/(βΩ))²)].                         (S7)
```
At βΩ=1 this is greater than `c²U/(164Ω)`, using the certified inequality π²<10. More generally, for `a=βΩ≤2π`,

```
log E exp[c∫|v|] ≥ c² UΩ²/(320τ³).                 (S8)
```
Indeed `1+(2π/a)²≤2(2π/a)²` and `32π²<320`. At τ=CΩ[log(1/Ω)]^p with fixed p≥0, the right side diverges as a positive constant times `1/(Ω log(1/Ω)^(3p))`. Dimensionally c is dimensionless, and the exponent in (S7)–(S8) is dimensionless.

This free projected dynamic process is also **exactly** the dynamic process of the full interacting isolated Holstein atom: the electronic weight depends only on b_0. Thus (S7) is not just a hypothetical process allowed by an upper sub-Gaussian increment estimate.

There is a fixed-t, fixed-density finite-square version for Holstein. On V sites perform an orthogonal oscillator change of variables with uniform coordinate `q=V^−1/2Σ_i b_i`. It couples to total electron number as `qN/sqrt(V)`. N commutes with the rest of the electron Hamiltonian, also along arbitrary histories of the other modes. The zero-mean dynamic part v_q therefore drops out of the electronic one-period weight and stays an independent free projected Gaussian in the full positive grand-canonical or balanced canonical measure. Write `L_i=L_q/sqrt(V)+L_i,rest` for the cosine functional of one site's dynamic field. Time translation invariance gives `E L_i,rest=0`. Independence and Jensen imply

```
E exp[c∫|v_i|] ≥ E exp(cL_i)
 ≥ exp[c² Var(L_q)/(2V)].                           (S9)
```
All moments needed here exist by harmonic confinement; if the absolute moment were infinite the lower bound would still be sufficient. This works for any fixed finite square, any t>0 and any balanced nontrivial filling, including density 1/2. The coefficient is 1/V, so (S9) is not a nonvanishing per-site thermodynamic lower bound. It does refute an Ω-uniform polynomial majorant asserted simultaneously for all finite volumes.

The barrier is deliberately limited: a scheme which pays the full absolute moment `E exp(c∫|v|)` as its Gram/activity majorant cannot make that majorant small at a linear-times-polylog scale using the slab budget alone. It is **not** a theorem that the optimal averaged Gram or determinant bound grows: Proposition 3 gives a determinant bound one for the averaged two-point kernel. Nor does it rule out correlations between static and dynamic field averages, or cancellations before absolute values.

## 5. Exact negative nonzero-Matsubara Hessian in square SSH at fixed density

### Proposition 5: uniform slow-field convexity fails below an explicit two-thirds scale

Use the square torus of side four, V=16, with spin sectors `N_up=N_down=4`, total density n=1/2. The bare t>0 and U>0 are arbitrary and fixed. Choose a three-site path 1–2–3 in this torus. Set every static bond coordinate to t to cancel its bare hopping, except bond 1–2, where set `b_12=t−δ`. Thus the frozen one-particle Hamiltonian is

```
h_0=−δK_12,
K_ij=|i><j|+|j><i|.
```
Let the only dynamic coordinate be `b_23(u)−t=x cos(νu)`, where `ν=2π/β=2πτ`, and choose the external frozen gap parameter `δ=ν`. These static coordinates remain bounded as τ tends to zero; no large Holstein potential or removed phonon mode is being assumed.

Define S_eff(x) to be the full oscillator history action minus the log of the exact balanced canonical electron history weight along this path. The weight is strictly positive at x=0, so its log is analytic for small real x. The history is an allowed smooth periodic history with zero temporal mean. Oscillator normalization constants and the static oscillator action have zero x-Hessian.

For the single dynamic bond, direct integration of cos² and sin² gives

```
S_ph''(0)=β[1+(ν/Ω)²]/(2U).                         (S10)
```
The factor 1/2 is essential: the oscillator action is `(2U)^−1∫[b²+dot b²/Ω²]`; each square integrates to β/2. The linear static term integrates to zero.

One spin has energies −δ,+δ and fourteen zero energies. Put q=βδ and

```
A=C(14,4)=1001,   B=C(14,3)=364,   C=C(14,2)=91,
Z_4=A+2B cosh(q)+C=B[3+2cosh(q)].                   (S11)
```
Configurations with neither, exactly one, or both nonzero-energy levels occupied give the three terms. Their occupation difference is

```
<n_->−<n_+>=2B sinh(q)/Z_4.                          (S12)
```
The hopping insertion K_23 connects site 3 to each of the ± eigenvectors with squared amplitude 1/2. Its nonzero bosonic Matsubara susceptibility in this exact canonical free-electron state is

```
χ_23(iν)=δ/(δ²+ν²) (<n_->−<n_+>).
```
For completeness this formula does not assume canonical Wick factorization. For an off-diagonal one-body transition a↔b with energy difference Δ, the two ordered thermal contributions combine into `2Δ(<n_a>−<n_b>)/(Δ²+ν²)` times the squared one-body matrix element. This follows either by integrating their spectral exponentials, or by the identity `[d_a†d_b,d_b†d_a]=n_a−n_b` and the KMS relation for these number-conserving transitions inside the N=4 sector. Summing the transitions −↔3 and 3↔+ cancels `<n_3>` and gives the displayed formula.

For a source x cos(νu), the one-spin second derivative of log Z is `βχ_23(iν)/2`: stationarity reduces the double integral against two cosines to β/2 times the Matsubara susceptibility. The two independent spin traces multiply, hence the balanced log-weight second derivative is βχ_23(iν), not βχ_23(iν)/2. With δ=ν and q=2π,

```
χ_23(iν)=sinh(2π)/[ν(3+2cosh(2π))],
S_eff''(0)/β=[1+(2πτ/Ω)²]/(2U)
             −sinh(2π)/[2πτ(3+2cosh(2π))].          (S13)
```

Here is an explicit sufficient condition for negativity, with no asymptotic substitution. Since `e^(2π)>4`,

```
4sinh(2π)>3+2cosh(2π),
χ_23(iν)>1/(4ν)=1/(8πτ).                           (S14)
```
The first inequality is equivalent to `e^(2π)−3e^(−2π)>3`; its left side is increasing in e^(2π) and already exceeds three at four. The elementary certified bounds π<4 and π³<32 give

```
4πτ+16π³τ³/Ω² <16τ+512τ³/Ω² ≤ U
```
whenever

```
τ≤U/32,          τ³≤UΩ²/1024.                       (S15)
```
These inequalities make the positive term in (S13) smaller than `1/(8πτ)`, which is strictly smaller than χ_23. Thus `S_eff''(0)<0` throughout (S15), including their boundaries.

For every fixed C>0 and fixed p≥0, τ=CΩ[log(1/Ω)]^p eventually satisfies (S15): τ/U tends to zero and `τ³/(UΩ²)=(C³/U)Ω[log(1/Ω)]^(3p)` tends to zero. A uniform convexity or positive-Hessian Gaussian-domination argument for the exact fermion-integrated slow-field action therefore fails at every such scale, even at balanced regular density 1/2 and with bounded frozen fields. It already fails in a nonzero Matsubara direction, not just a static double-well direction.

An exact RG which integrates all electron modes into a phonon-only action must reproduce (S13) and cannot turn this same action into a uniformly convex one. This does **not** exclude keeping low-energy fermions as explicit variables, or proving that these frozen configurations have an adequate annealed exception cost. The small fermion energy δ=2πτ is precisely a low-energy mode, so integrating only higher scales while retaining it has not yet produced the phonon-only convex action tested here. A pointwise quadratic lower bound on the action, or convexity of a distinct Legendre effective functional, is not a positive Hessian bound on S_eff and is not refuted by this proposition. No claim of an SSH physical ordered phase or a sharp physical exponent is made.

## 6. Focused annealed-block attempt: the available spectral lemma and the missing gluing

### Proposition 6: a frozen two-electron Combes–Thomas estimate

Let h be any static, self-adjoint nearest-neighbor one-electron Hamiltonian on a finite square, with arbitrary real diagonal and hopping modulus at most J>0. On the antisymmetric two-electron Hilbert space let `h_2=h⊗I+I⊗h`, restricted to that space. Pair configurations are connected by moving either electron along one nearest-neighbor link; their row degree is at most eight and hopping modulus at most J. Let d(p,q) be their graph distance. For z with `|Im z|=τ`, set

```
a=log(1+τ/(16J)).
```
Then

```
|[(h_2−z)^−1]_(p,q)| ≤ (2/τ)exp[−a d(p,q)].        (S16)
```
Proof: weight configurations by `W_p=exp[a d(p,q)]`. A hop changes the distance by at most one, so every absolute row and column sum of `W h_2 W^−1−h_2` is at most `8J(e^a−1)=τ/2`. The diagonal is unchanged. The row/column Schur bound gives the same operator-norm estimate. Self-adjointness gives `||(h_2−z)^−1||≤1/τ`. The Neumann series for the weighted resolvent then gives norm at most 2/τ, and extracting its (p,q) entry proves (S16). The exterior-space signs change no absolute bound. In disconnected components the off-diagonal entry is zero. Its range is of order J/τ when τ/J is small. All energies and the exponent have the stated units.

This estimate needs no uniform electronic gap and includes arbitrarily low frozen electron energies. It is a useful frozen-field lemma, not a bound on the exact interacting pair susceptibility. In particular the latter is a bosonic zero-frequency response, not automatically the resolvent at iτ in (S16). A Wick or other justified insertion representation would be required to connect them. Free frozen electronic thermal pair correlations can be handled using the parent static covariance estimates; after averaging histories the missing interacting vertices from Proposition 1 reappear.

### Proposition 7: being good does not give a pointwise charged temporal gluing estimate

On the same square SSH torus cancel every static hopping by b_e=t. On one bond 1–2 choose the smooth history

```
h(u)=−r sin(2πu/β) K_12,
β=a/Ω,  r=ρsqrt(UΩ),  a,ρ>0 fixed.
```
All its one-particle evolutions commute. Its full period matrix is I because the sine integral is zero. Its half-period evolution is `exp(ηK_12)` with

```
η=rβ/π=(ρa/π)sqrt(U/Ω).                            (S17)
```
The oscillator increment across any subinterval of this history is at most 2r. A whole thermal-circle slab, starting where the sine is zero, has excursion r, so it is good for the natural cutoff 2r. This can be in the range of the parent's certified slab budget: the associated `x=(2r)²/(32Dβ)=ρ²/(8a)` exceeds 16 when ρ²≥128a. Its derivative is at most `2πr/β`, hence tends to zero like Ω^(3/2). This is not a rapid switching construction.

Grand canonically its normalized one-spin contraction from site 1 at time zero to site 2 at β/2 is `sinh(η)/2`, because M=I and the equal-time vacancy matrix is I/2. Its two-spin onsite pair-creation/annihilation transport insertion is the product of the spin contractions, hence `sinh²(η)/4` up to the fixed harmless fermion-ordering sign. The electronic history weight is `det(1+I)²=4^V`, independent of r and Ω.

There is an exact balanced canonical version, so no Gaussian canonical assumption is necessary. In the k-particle one-spin sector, M=I and its normalized trace is the identity state on the exterior k-space. Its vacancy matrix is `(1−k/V)I`. Conjugating the annihilation operator by the half-period evolution gives the contraction `(1−k/V)sinh(η)`. The two independent k-particle spin traces therefore give the charged two-time pair insertion

```
(1−k/V)²sinh²(η),
```
and the balanced electronic history weight is `C(V,k)²`, again constant. For V=16,k=4 the prefactor is 9/16 and the density is 1/2.

This insertion grows faster than any polynomial in 1/Ω even though the history is good and each instantaneous frozen h satisfies (S16). Thus instantaneous spectral bounds plus a good/bad amplitude flag do not give a uniform pointwise *charged temporal boundary/gluing* bound with polynomial prefactor. This is not a counterexample to the physical equal-time spatial pair-order theorem: the electron density at equal times here is the identity state, and the test is an intermediate charged two-time insertion. It also does not refute an annealed bound on these insertions; integrating the field before absolute values may compensate the amplification.

The focused annealed attempt therefore reaches a precise gap, rather than a false closure. One needs a **charged-block inequality** which integrates these temporally open insertions under the exact coupled weight, retains the low-energy electron boundary data, and bounds their *spatial* transport/connected attachment cost. The source-free joint excursion bound (parent (33)) controls the probability of bad flags, but not the insertion-weighted contribution within the good set. Proposition 7 shows why replacing the desired annealed inequality by its pointwise good-set counterpart is false. Proposition 4 shows why paying an unconditional absolute integrated-field exponential is insufficient. Neither proposition rules out an insertion-sensitive annealed comparison. No such comparison is established here.

## 7. Physics sanity check: what the atomic Franck–Condon factor really proves

### Proposition 8: exact atomic low-sideband pair weight

For the centered atom of Proposition 1, the oscillator ground wavefunctions in occupancies n=0 and n=2 are Gaussians of probability variance `UΩ/2`, centered at +U and −U. Explicitly their normalized wavefunctions are

```
φ_c(b)=(πUΩ)^(-1/4) exp[−(b−c)²/(2UΩ)].
```
Completing the square in their inner product gives the electronic pair ground-to-ground matrix element

```
|<0,φ_(+U)|c_down c_up|2,φ_(−U)>|=exp(−U/Ω),
```
and its spectral weight is `exp(−2U/Ω)`. The amplitude and probability are different; neither is a measured value.

Translating a Gaussian ground state by 2U is an oscillator displacement of squared coherent amplitude `Λ=2U/Ω`. Expanding the exponential displacement on the vacuum gives amplitudes `e^(−Λ/2)(±sqrt(Λ))^m/sqrt(m!)`. This expansion follows by the oscillator commutator and the elementary identity `exp(α(a†−a))=exp(−α²/2)exp(αa†)exp(−αa)`. Squaring shows that the pair transition to the m-th oscillator sideband has exact weight

```
w_m=exp(−Λ)Λ^m/m!,    Λ=2U/Ω,    Σ_m w_m=1.        (S18)
```
The sector ground energies coincide, and that sideband has energy mΩ above them. For any fixed A≥0, the zero-temperature atomic pair spectral weight up to AΩ is exactly

```
exp(−2U/Ω) Σ_(m=0)^floor(A) (2U/Ω)^m/m!.            (S19)
```
It is exponentially small for fixed A and U>0 as Ω tends to zero, despite the full pair operator having norm one. Its dominant sidebands have oscillator number of order U/Ω and therefore energy of order U, not Ω. These are exact atomic spectral statements, not a finite-density lattice tunnelling bound.

**[INFERENCE, limited physical heuristic.]** In a regime with well-separated adiabatic distortion minima and a persistent barrier, (S18) suggests exponentially small coherent bipolaron transport and possibly a pairing/coherence scale smaller than linear in Ω. A theorem would need to control motion between such minima at fixed t, include soft phonon excitations, and relate the finite-temperature stiffness or spatial electronic pair correlations to that motion. Those hypotheses have not been proved for arbitrary fixed coupling and regular density. Weakly distorted/itinerant configurations cannot be excluded by the atomic overlap, so there is no supported universal exponential claim.

The current parent evidence makes this limitation essential. LINEAR2/RIGOR §1 proves an order-t kinetic expectation at regular density 1/2 for `0<U≤t/4, τ≤t/64`, uniformly in Ω; vacuum narrowing cannot bound that kinetic term. Its §2 proves only the t=0 second-order pair-transfer coefficient

```
P_2=2t²e^(−U/Ω)∫_0^∞ exp[−Us−(U/Ω)e^(−Ωs)] ds,
P_2 ~ t²sqrt(2π/(UΩ))exp(−2U/Ω),
```
not a fixed-t many-body stiffness bound. The complement includes oscillator excitations at energy Ω, so its inverse cannot uniformly be replaced by 1/U. No extension of that coefficient to a lattice transport ceiling has been made here. A genuine polaron-tunnelling route remains an alternative research question, not a demonstrated replacement theorem.

## 8. Certification plan and limits

The retained replay uses exact rational arithmetic for the canonical binomial coefficients, the cosine kinetic/spin factors, a Machin arctangent certificate for π, and the polynomial sufficient conditions. Its interval checks independently evaluate the atomic cumulant, the finite Lehmann thermal normalization, a chronological 2×2 covariance determinant, the SSH spectral susceptibility/Hessian at an admitted parameter point, and the atomic sideband weight. All mutants have named rejection guards. RESULT.md records the actual command and exits after execution.

For the π certificate put θ=arctan(1/5), φ=arctan(1/239). Double-angle identities give `tan(2θ)=5/12` and `tan(4θ)=120/119`; the subtraction identity gives `tan(4θ−φ)=1`. Since `0<4θ−φ<1<π/2`, this angle is π/4. Hence `π=16θ−4φ`. Alternating arctangent Taylor sums give exact rational enclosing bounds for each θ and φ and thus for π. The replay checks π<22/7, π²<10 and π³<32. Its lower exponential certificate `e²>1+2+2+8/6=19/3>4`, together with π>1, supplies the strict inequality needed in (S14).

The finite replay is not a formal proof of all-order chronological determinants, infinite-dimensional heat-trace Hölder, the SSH susceptibility theorem at every parameter, an annealed charged-block estimate, or a physical pairing ceiling. The analytical proofs above state precisely what is and is not proved.
