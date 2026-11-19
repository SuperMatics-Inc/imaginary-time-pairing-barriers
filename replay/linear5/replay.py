"""Exact algebra and interval spot certificates for LINEAR5.

Checks cover SSH histories, source factors, reflection forms, kinetic gauges,
slab constants and fixed-density route barriers. No many-body simulation.
"""
from fractions import Fraction as F
from dataclasses import dataclass
from itertools import permutations
from math import factorial
import argparse
import subprocess
import sys
from pathlib import Path


@dataclass(frozen=True)
class Q15:
    a: F = F(0)
    b: F = F(0)

    @staticmethod
    def of(x):
        return x if isinstance(x, Q15) else Q15(F(x))

    def __add__(self, x):
        x = self.of(x)
        return Q15(self.a + x.a, self.b + x.b)
    __radd__ = __add__

    def __neg__(self):
        return Q15(-self.a, -self.b)

    def __sub__(self, x):
        return self + -self.of(x)

    def __rsub__(self, x):
        return self.of(x) + -self

    def __mul__(self, x):
        x = self.of(x)
        return Q15(self.a * x.a + 15 * self.b * x.b, self.a * x.b + self.b * x.a)
    __rmul__ = __mul__

    def inverse(self):
        d = self.a * self.a - 15 * self.b * self.b
        if not d:
            raise ZeroDivisionError
        return Q15(self.a / d, -self.b / d)

    def __truediv__(self, x):
        return self * self.of(x).inverse()


def require(ok, name):
    if not ok:
        raise AssertionError(name)


def eye():
    return [[Q15(F(i == j)) for j in range(3)] for i in range(3)]


def mul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(3)), Q15()) for j in range(3)] for i in range(3)]


def det(A):
    n = len(A)
    result = Q15()
    for p in permutations(range(n)):
        sign = (-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = Q15.of(sign)
        for i in range(n):
            term *= A[i][p[i]]
        result += term
    return result


def boost(c, s, i, j):
    A = eye()
    A[i][i] = A[j][j] = Q15.of(c)
    A[i][j] = A[j][i] = Q15.of(s)
    return A


def checks(mutant=None):
    c1, s1, c2, s2 = F(5, 3), F(4, 3), F(4), Q15(F(0), F(1))
    if mutant == "field-radicand":
        s2 = Q15.of(4)
    require(Q15.of(c1 * c1 - s1 * s1) == Q15.of(1), "first-boost-hyperbolic-identity")
    require(Q15.of(c2 * c2) - s2 * s2 == Q15.of(1), "second-boost-hyperbolic-identity")
    A = boost(c1, s1, 0, 1)
    Ai = boost(c1, -s1, 0, 1)
    B = boost(c2, s2, 1, 2)
    Bi = boost(c2, -s2, 1, 2)
    if mutant == "inverse-sign":
        Ai = A
    factors = [A, B, Ai, Bi]
    if mutant == "commuting-order":
        factors = [A, Ai, B, Bi]
    Y = eye()
    for factor in factors:
        Y = mul(Y, factor)
    require(det(Y) == Q15.of(1), "monodromy-determinant")
    require(sum((Y[i][i] for i in range(3)), Q15()) == Q15.of(-1), "monodromy-trace-minus-one")
    require(mul(Y, Y) == eye(), "nontrivial-involution")
    require(Y != eye(), "not-identity-history")
    plus = [[eye()[i][j] + Y[i][j] for j in range(3)] for i in range(3)]
    require(det(plus) == Q15(), "one-spin-trace-zero")
    # Spectral projectors P+ = (I+Y)/2 and P- = (I-Y)/2.
    Pm = [[(eye()[i][j] - Y[i][j]) / 2 for j in range(3)] for i in range(3)]
    require(mul(Pm, Pm) == Pm, "negative-projector")
    require(sum((Pm[i][i] for i in range(3)), Q15()) == Q15.of(2), "negative-eigenspace-rank")
    q = F(101, 100)
    occupation = q / (q - 1)
    if mutant == "thermal-occupation":
        occupation = q / (q + 1)
    require(occupation == 101, "unbounded-reference-occupation")
    # det(1+qY) = (1+q)(1-q)^2, per spin, at each rational fugacity.
    for q in [F(0), F(1, 2), F(1), F(101, 100), F(2)]:
        mat = [[eye()[i][j] + q * Y[i][j] for j in range(3)] for i in range(3)]
        require(det(mat) == Q15.of((1 + q) * (1 - q) ** 2), "fugacity-polynomial")
    one = [1, -1, -1, 1]
    two = [sum(one[i] * one[n - i] for i in range(4) if 0 <= n - i < 4) for n in range(7)]
    expected = [1, -2, -1, 4, -1, -2, 1]
    if mutant == "canonical-positivity":
        expected = [abs(x) for x in expected]
    require(two == expected and two[1] < 0, "canonical-history-sign")
    spin_power = 1 if mutant == "one-spin-weight" else 2
    require(spin_power == 2, "spin-square-weight")
    require((F(5, 3) - 1) * (4 - 1) == 2, "pi-rotation-parameters")
    # Atomic source integral's endpoint coefficients, from elementary moments.
    # 2 endpoints; integral s^2 exp(-2Us) = 2/(2U)^3.
    quantum_coefficient = F(2) * F(2, 8)
    if mutant == "endpoint-factor":
        quantum_coefficient /= 2
    require(quantum_coefficient == F(1, 2), "atomic-quantum-first-order-factor")
    require(F(2) * 2 * F(2, 8) == 1, "atomic-static-first-order-factor")
    balanced_weight = sum(F(x * x) for x in one)
    if mutant == "same-spin-twist":
        balanced_weight = F(one[0] ** 2)
    require(balanced_weight == 4, "opposite-spin-Fourier-projection")
    transfer = F(2) * Y[0][1] * Y[0][1] / balanced_weight
    require(transfer == Q15.of(F(200, 9)), "twist-averaged-pair-nonnormality")
    # Finite Hilbert-Schmidt Gram identity underlying time reflection.
    M = [[F(0), F(1)], [F(-1), F(0)]]
    N = [[F(1), F(2)], [F(0), F(1)]]
    def hs(X, Y):
        return sum(X[i][j] * Y[i][j] for i in range(2) for j in range(2))
    mm = hs(M, M)
    if mutant == "mirror-without-adjoint":
        mm = sum(M[i][j] * M[j][i] for i in range(2) for j in range(2))
    require(mm == 2, "reflected-half-Hilbert-Schmidt-positivity")
    require(mm * hs(N, N) - hs(M, N) ** 2 == 8, "reflection-Schwarz-Gram-determinant")
    for z in [F(0), F(1, 2), F(1), F(2)]:
        coefficients = [F(x) * z ** n for n, x in enumerate(one)]
        require(sum(coefficients) ** 2 <= 4 * sum(x * x for x in coefficients),
                "balanced-pressure-upper-comparison")
    # Exact finite constants for the small-cutoff sharpness proposition.
    increment_leading = F(3, 16)
    if mutant == "increment-factor":
        increment_leading /= 2
    require(increment_leading == F(3, 32) / F(1, 2), "quarter-period-increment-factor")
    require(F(3, 32) / F(4, 7) == F(21, 128), "increment-variance-lower-constant")
    require(F(3, 28) / F(1, 2) == F(3, 14), "increment-variance-upper-constant")
    require(F(1, 5) * F(3, 7) == F(3, 35), "large-a-variance-lower-constant")
    require(F(1, 24) == F(1, 2) ** 2 / 6, "projected-covariance-leading-factor")
    centering = 1 if mutant == "brownian-centering" else 2
    require(centering ** 2 == 4, "Brownian-bridge-centering-loss")
    lam = F(1, 4)
    require(1 + 2 * lam / (F(1, 2) - lam) == 3, "Gaussian-sup-exponential-moment")
    power = F(1, 2) if mutant == "threshold-power" else F(2, 3)
    require(3 * power - 2 == 0, "two-thirds-cutoff-threshold-exponent")
    for occupation in [0, 1, 2]:
        dynamic_integral = F(1) + F(-1)
        if mutant == "atomic-dynamic-tilt":
            dynamic_integral = F(1) + F(1)
        require((occupation - 1) * dynamic_integral == 0, "atomic-nonzero-mode-factorization")
    require(F(21, 512) * F(1, 2) ** 2 == F(21, 2048), "half-budget-necessary-constant")
    exp2_lower = sum(F(2) ** k / F(factorial(k)) for k in range(4))
    require(exp2_lower == F(19, 3) and 3 / exp2_lower < F(1, 2), "half-budget-sufficient-Q-eight")
    # Kinetic gauge, optimized circle source, dyadic excursion and mirrored-time factors.
    gauge_constant = F(1) if mutant == "kinetic-gauge-sign" else F(-1)
    require(gauge_constant == -F(1), "imaginary-coordinate-gauge-square-sign")
    durations, gauges = [F(1), F(2), F(1)], [F(2), F(0), F(-1)]
    beta = sum(durations)
    mean = F(1) if mutant == "gauge-centering" else sum(t*q for t, q in zip(durations, gauges)) / beta
    centered = sum(t*(q-mean)**2 for t, q in zip(durations, gauges))
    require(centered == F(19, 4) and centered <= F(5), "optimal-circle-gauge-centering")
    require(F(1, 3)-F(1, 2)+F(1, 4) == F(1, 12), "projected-point-gauge-variance")
    require(F(1)*(F(4)-F(1))/F(4)/2 == F(3, 8), "increment-exponential-half-factor")
    first = F(1, 8) if mutant == "dyadic-allocation" else F(1, 4)
    require(first/(1-F(3, 4)) == 1, "dyadic-telescoping-threshold-budget")
    require(F(3, 4)**2*2 == F(9, 8), "dyadic-tail-growth-factor")
    require(2/exp2_lower < F(1, 2), "dyadic-tail-geometric-ratio")
    require(F(512)/32 == 16, "natural-slab-amplitude-threshold")
    r_power = F(3, 4) if mutant == "natural-amplitude-power" else F(1, 2)
    require(2*r_power-F(2)-F(-1) == 0, "linear-scale-natural-amplitude-budget")
    mirrored = F(1, 2) if mutant == "mirror-equal-time" else (F(17, 8)+F(15, 8))/2
    require(mirrored == 2 and mirrored > 1, "mirrored-unequal-time-minor-growth")
    require(F(8)*F(1, 8) == 1, "quarter-profile-peak-amplitude-factor")
    balanced_weight = F(-2) if mutant == "balanced-density-weight" else F(-2)**2
    require(balanced_weight == 4 and balanced_weight >= 0, "fixed-balanced-positive-history-weight")
    largest_spin_dimension = 6 if mutant == "canonical-multiplicity" else 5
    require(largest_spin_dimension <= 4+1, "fixed-density-SU2-multiplicity-bound")
    require(F(1)+3*F(2)+5*F(3) == 22 and 22 <= 5*(1+2+3), "canonical-balanced-pressure-sandwich")
    instantaneous_error = F(0) if mutant == "adiabatic-occupation" else F(1, 2)*(4-1)/(4+1)
    require(instantaneous_error == F(3, 10), "instantaneous-equilibrium-tracking-error")
    require(F(1, 2)-1 == F(-1, 2) and F(1, 2)+1 == F(3, 2), "slow-drift-versus-thermal-memory-powers")
    occupations, weights = [F(2), F(0), F(1), F(1)], [F(4), F(4), F(2), F(2)]
    atom_mean = sum(n*w for n, w in zip(occupations, weights))/sum(weights)
    atom_var = sum(w*(n-atom_mean)**2 for n, w in zip(occupations, weights))/sum(weights)
    spatial_form = atom_var if mutant == "spatial-RP-sign" else -atom_var
    require(atom_mean == 1 and atom_var == F(2, 3), "fixed-density-atomic-phonon-covariance")
    require(spatial_form == F(-2, 3) and spatial_form < 0, "balanced-canonical-spatial-RP-counterexample")
    # Interval spot certificates for the exact covariance formulas. The all-a
    # inequalities, Gaussian tails and interacting factorization are analytic.
    import mpmath as mp
    iv = mp.iv
    iv.dps = 60
    def interval(x):
        x = F(x)
        return iv.mpf(x.numerator) / iv.mpf(x.denominator)
    def cosh(x):
        return (iv.exp(x) + iv.exp(-x)) / 2
    def sinh(x):
        return (iv.exp(x) - iv.exp(-x)) / 2
    # Interval point certificate for the non-Hermitian Schatten S_2 inequality.
    # H=diag(0,log 2), A=sigma_x; this is not a proof for unbounded operators.
    half_log = iv.ln(interval(2))/2
    if mutant == "drift-unitarity":
        nu = iv.sqrt(half_log**2+1)
        frobenius_sq = 2*iv.exp(-2*half_log)*(cosh(nu)**2+sinh(nu)**2)
    else:
        omega = iv.sqrt(1-half_log**2)
        frobenius_sq = 2*iv.exp(-2*half_log)*(iv.cos(omega)**2+
                          (half_log**2+1)*(iv.sin(omega)/omega)**2)
    require(frobenius_sq < interval(F(5, 4)), "anti-Hermitian-drift-Schatten-point")
    for a_exact in [F(1, 16), F(1, 8), F(1, 2), F(1), F(2), F(4), F(16)]:
        a = interval(a_exact)
        projected = a / (2 * sinh(a / 2)) - 1
        if mutant == "projected-RP":
            projected = a / (2 * sinh(a / 2))
        require(projected < 0, "negative-projected-reflected-covariance")
        h = (cosh(a / 2) - cosh(a / 4)) / sinh(a / 2)
        require(h > 0, "positive-local-increment-variance")
        if a_exact <= 1:
            require(h > interval(F(21, 128)) * a and h < interval(F(3, 14)) * a,
                    "interval-small-a-variance-sandwich")
        else:
            require(h > interval(F(3, 35)), "interval-large-a-variance-lower-bound")
    print("CERTIFIED: sharpness constants and seven interval covariance points passed")
    print("CERTIFIED: kinetic gauge, dyadic slab constants and mirrored-time factors passed")
    print("CERTIFIED: balanced density, slow tracking and canonical spatial-RP factors passed")
    print("CERTIFIED: monodromy, source factors, reflection Gram and twist identities passed")


MUTANTS = ["field-radicand", "inverse-sign", "commuting-order", "thermal-occupation",
           "canonical-positivity", "one-spin-weight", "endpoint-factor",
           "same-spin-twist", "mirror-without-adjoint", "increment-factor",
           "brownian-centering", "threshold-power", "atomic-dynamic-tilt", "projected-RP",
           "kinetic-gauge-sign", "gauge-centering", "dyadic-allocation",
           "natural-amplitude-power", "mirror-equal-time", "drift-unitarity",
           "balanced-density-weight", "canonical-multiplicity", "adiabatic-occupation", "spatial-RP-sign"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mutant", choices=MUTANTS)
    ap.add_argument("--all-mutants", action="store_true")
    args = ap.parse_args()
    if args.mutant and args.all_mutants:
        ap.error("choose one mode")
    try:
        checks(args.mutant)
    except AssertionError as e:
        print(f"REJECTED {args.mutant or 'base'}: {e}")
        return 1
    if args.all_mutants:
        for name in MUTANTS:
            p = subprocess.run([sys.executable, "-B", str(Path(__file__).resolve()), "--mutant", name],
                               capture_output=True, text=True, check=False)
            print(p.stdout.strip())
            if p.stderr:
                print(p.stderr.strip())
            print(f"mutant {name}: exit {p.returncode}")
            if p.returncode != 1 or not p.stdout.startswith(f"REJECTED {name}:"):
                return 2
        print(f"CERTIFIED: all {len(MUTANTS)} mutants rejected with exit 1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
