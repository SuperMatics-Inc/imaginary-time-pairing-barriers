"""Finite exact-rational/interval certificates for STAGE3; not a many-body simulation.

All mutants deliberately alter a named scientific factor or inference. No network,
files written, eigenvalue fit, or physical pairing-temperature estimate is used.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse
import subprocess
import sys
import mpmath as mp

MUTANTS = (
    "wick-average", "canonical-wick", "lehmann-weight", "thermal-normalization",
    "chronological-sign", "fourier-variance", "ssh-canonical", "ssh-spin",
    "ssh-kinetic", "hessian-sign", "good-block", "franck-condon",
)


class Rejected(Exception):
    pass


def require(condition, guard):
    if not condition:
        raise Rejected(guard)


def atan_enclosure(x, count=20):
    """Alternating exact-rational bounds, with count an even term count."""
    assert 0 < x < 1 and count % 2 == 0
    lower = sum(((-1) ** k * x ** (2 * k + 1) / (2 * k + 1)
                 for k in range(count)), F(0))
    upper = lower + x ** (2 * count + 1) / (2 * count + 1)
    return lower, upper


def iv_fraction(q):
    q = F(q)
    return mp.iv.mpf(q.numerator) / q.denominator


def enclose_zero(value):
    return value.a <= 0 and value.b >= 0


def check(mutant=None):
    mp.iv.dps = 70
    iv = mp.iv
    one, zero = iv.mpf(1), iv.mpf(0)
    pi = iv.pi

    # Rational Machin proof and constants used in the analytic sufficient condition.
    tan2 = 2 * F(1, 5) / (1 - F(1, 5) ** 2)
    tan4 = 2 * tan2 / (1 - tan2 ** 2)
    tan_difference = (tan4 - F(1, 239)) / (1 + tan4 * F(1, 239))
    require(tan2 == F(5, 12) and tan4 == F(120, 119)
            and tan_difference == 1, "machin-angle-identity")
    l5, u5 = atan_enclosure(F(1, 5))
    l239, u239 = atan_enclosure(F(1, 239))
    pi_lower, pi_upper = 16 * l5 - 4 * u239, 16 * u5 - 4 * l239
    require(3 < pi_lower < pi_upper < F(22, 7), "rational-pi-enclosure")
    require(pi_upper ** 2 < 10 and pi_upper ** 3 < 32,
            "pi-square-cube-constants")
    require(pi > iv_fraction(pi_lower) and pi < iv_fraction(pi_upper),
            "interval-pi-matches-rational-certificate")
    exp2_lower = sum((F(1), F(2), F(2), F(8, 6)), F(0))
    require(exp2_lower == F(19, 3) and exp2_lower > 4,
            "strict-matsubara-susceptibility-lower-bound")

    # Exact atomic electronic reduction at beta*U=4, all Omega.
    w = iv.exp(2)
    z_atom = 2 * (w + 1)
    occupation = (w + 1) / z_atom
    double_occupation = w / z_atom
    cumulant = double_occupation - occupation ** 2
    claimed_cumulant = (w - 1) / (4 * (w + 1))
    if mutant == "wick-average":
        claimed_cumulant = zero
    require(enclose_zero(cumulant - claimed_cumulant)
            and claimed_cumulant > 0,
            "averaging-does-not-preserve-wick")
    canonical_double = F(1, 4) if mutant == "canonical-wick" else F(0)
    require(canonical_double == 0 and F(1, 2) ** 2 - canonical_double == F(1, 4),
            "canonical-fixed-number-is-not-gaussian")

    # A one-orbital Lehmann measure must integrate to the CAR identity.
    e_minus = iv.exp(-2)
    lehmann_weight = ((1 - e_minus) if mutant == "lehmann-weight"
                      else (1 + e_minus)) / (1 + e_minus)
    require(enclose_zero(lehmann_weight - 1), "lehmann-sum-not-difference")
    denominator = one if mutant == "thermal-normalization" else (1 + iv.exp(2))
    negative_energy_covariance = iv.exp(1) / denominator
    require(negative_energy_covariance > 0 and negative_energy_covariance < 1,
            "thermal-normalization-bounds-negative-energy-branch")

    # Chronological 2x2 covariance at epsilon=2, beta=1, times 1/4,3/4;0,1/2.
    thermal_denominator = 1 + iv.exp(-2)
    c11 = c22 = iv.exp(iv_fraction(F(-1, 2))) / thermal_denominator
    c21 = iv.exp(iv_fraction(F(-3, 2))) / thermal_denominator
    c12 = -c21 if mutant != "chronological-sign" else c21
    minor = c11 * c22 - c12 * c21
    expected_minor = 1 / (iv.exp(1) + iv.exp(-1))
    require(enclose_zero(minor - expected_minor) and minor > 0 and minor < 1,
            "chronological-antiperiodic-sign")

    # Projected Fourier moment and the explicit exponential-majorant lower bound.
    omega, beta, u = F(1, 100), F(100), F(1)
    nu = 2 * pi / iv_fraction(beta)
    d1 = iv_fraction(u * omega ** 2 / beta) / (iv_fraction(omega) ** 2 + nu ** 2)
    variance_from_fourier = 2 * d1 * iv_fraction(beta / 2) ** 2
    beta_factor = F(1, 1) / beta if mutant == "fourier-variance" else beta
    variance_formula = iv_fraction(u * omega ** 2 * beta_factor / 2) / (
        iv_fraction(omega) ** 2 + nu ** 2)
    require(enclose_zero(variance_from_fourier - variance_formula),
            "two-conjugate-fourier-modes-and-beta-factor")
    log_mgf_lower = variance_from_fourier / 2
    require(log_mgf_lower > iv_fraction(F(100, 164)),
            "linear-scale-absolute-moment-is-exponential")
    require(F(1, 32) * F(1, 10) == F(1, 320),
            "small-a-cubic-scale-constant")

    # Canonical 16-site N_up=N_down=4 susceptibility, independently from occupations.
    a, b, c = comb(14, 4), comb(14, 3), comb(14, 2)
    if mutant == "ssh-canonical":
        a += 1
    require((a, b, c) == (1001, 364, 91) and a + c == 3 * b,
            "sixteen-site-canonical-binomial-weights")
    tau, omega, u = F(1, 64), F(1, 16), F(1)
    require(tau <= u / 32 and tau ** 3 <= u * omega ** 2 / 1024,
            "admitted-hessian-parameter-point")
    frequency = 2 * pi * iv_fraction(tau)
    q = 2 * pi
    eq, emq = iv.exp(q), iv.exp(-q)
    coshq, sinhq = (eq + emq) / 2, (eq - emq) / 2
    z = a + 2 * b * coshq + c
    n_minus, n_plus = (b * eq + c) / z, (b * emq + c) / z
    n_three = (comb(13, 3) + 2 * comb(13, 2) * coshq + comb(13, 1)) / z
    delta = frequency
    spectral_chi = delta / (delta ** 2 + frequency ** 2) * (
        (n_minus - n_three) + (n_three - n_plus))
    formula_chi = sinhq / (frequency * (3 + 2 * coshq))
    require(enclose_zero(spectral_chi - formula_chi)
            and formula_chi > 1 / (4 * frequency),
            "canonical-transition-spectral-susceptibility")
    cosine_square_average = F(1, 2)
    kinetic_factor = F(1) if mutant == "ssh-kinetic" else cosine_square_average
    require(kinetic_factor == F(1, 2), "single-bond-cosine-kinetic-half-factor")
    spin_factor = F(1, 2) if mutant == "ssh-spin" else F(1)
    require(spin_factor == 2 * cosine_square_average,
            "balanced-two-spin-susceptibility-factor")
    harmonic_hessian_over_beta = iv_fraction(kinetic_factor / u) * (
        1 + (frequency / iv_fraction(omega)) ** 2)
    electron_hessian_over_beta = iv_fraction(spin_factor) * formula_chi
    sign = 1 if mutant == "hessian-sign" else -1
    hessian_over_beta = harmonic_hessian_over_beta + sign * electron_hessian_over_beta
    require(hessian_over_beta < 0, "negative-nonzero-matsubara-hessian")
    require(16 * tau <= u / 2 and 512 * tau ** 3 / omega ** 2 <= u / 2,
            "explicit-all-parameter-sufficient-constants")
    # Exact CT row/column perturbation budget on the two-electron graph.
    j, ct_tau = F(2), F(1, 7)
    e_a_minus_one = ct_tau / (16 * j)
    require(8 * j * e_a_minus_one == ct_tau / 2,
            "two-electron-combes-thomas-degree-eight")

    # The all-good long slab still has a large charged temporal insertion.
    eta = iv.mpf(4)
    sinh_eta = (iv.exp(eta) - iv.exp(-eta)) / 2
    vacancy = F(3, 4)
    canonical_transport = iv_fraction(vacancy ** 2) * sinh_eta ** 2
    claimed_transport = zero if mutant == "good-block" else canonical_transport
    require(enclose_zero(canonical_transport - claimed_transport)
            and claimed_transport > 100,
            "good-amplitude-flag-does-not-control-charged-gluing")
    require(F(16) ** 2 / (8 * F(1)) >= 16,
            "smooth-good-history-lies-in-slab-budget-range")

    # Atomic Gaussian overlap: distinguish amplitude from probability.
    u, omega = F(1), F(1, 8)
    separation_squared = (2 * u) ** 2
    gaussian_probability_variance = u * omega / 2
    amplitude_exponent = -separation_squared / (8 * gaussian_probability_variance)
    claimed_exponent = -2 * u / omega if mutant == "franck-condon" else -u / omega
    require(amplitude_exponent == claimed_exponent == F(-8),
            "franck-condon-amplitude-versus-spectral-weight")
    lam = 2 * u / omega
    low_sideband_polynomial = 1 + lam + lam ** 2 / 2
    require(lam == 16 and low_sideband_polynomial == 145,
            "atomic-sideband-poisson-factor")
    low_sideband_weight = iv.exp(iv_fraction(-lam)) * iv_fraction(low_sideband_polynomial)
    require(low_sideband_weight > 0 and low_sideband_weight < iv_fraction(F(1, 10000)),
            "atomic-low-energy-weight-is-not-total-operator-weight")

    print("CERTIFIED: rational pi, Fourier constants, canonical/spin/kinetic factors passed")
    print("CERTIFIED: atomic cumulant, Lehmann normalization, chronological minor passed")
    print("CERTIFIED: canonical SSH susceptibility and admitted negative Hessian passed")
    print("CERTIFIED: CT budget, charged gluing witness and atomic sideband checks passed")
    print("INTERVAL: atomic cumulant at beta*U=4 =", claimed_cumulant)
    print("INTERVAL: SSH S_eff''/beta at U=1,tau=1/64,Omega=1/16 =", hessian_over_beta)
    print("INTERVAL: atomic sideband weight m<=2 at U/Omega=8 =", low_sideband_weight)
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutant", choices=MUTANTS)
    parser.add_argument("--all-mutants", action="store_true")
    args = parser.parse_args()
    try:
        check(args.mutant)
    except Rejected as exc:
        print(f"REJECTED {args.mutant or 'base'}: {exc}")
        return 1
    if args.mutant:
        print(f"UNREJECTED mutant: {args.mutant}")
        return 2
    if args.all_mutants:
        for name in MUTANTS:
            proc = subprocess.run(
                [sys.executable, "-B", str(Path(__file__).resolve()), "--mutant", name],
                capture_output=True, text=True, check=False,
            )
            print(proc.stdout.strip())
            if proc.stderr:
                print(proc.stderr.strip())
            print(f"mutant {name}: exit {proc.returncode}")
            if proc.returncode != 1 or not proc.stdout.startswith(f"REJECTED {name}:"):
                return 1
        print(f"CERTIFIED: all {len(MUTANTS)} mutants rejected with exit 1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
