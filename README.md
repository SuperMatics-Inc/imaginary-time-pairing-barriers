# Replay supplement for Imaginary-time estimates and method barriers beyond the two-thirds adiabatic bound

Davis Rens, SuperMatics Inc. ([ORCID](https://orcid.org/0000-0001-5839-5197)).

This repository accompanies the preprint [Imaginary-time estimates and method barriers beyond the two-thirds adiabatic bound](https://doi.org/10.5281/zenodo.24374043). The replay supplement is [10.5281/zenodo.24374045](https://doi.org/10.5281/zenodo.24374045). The source repository is [imaginary-time-pairing-barriers](https://github.com/SuperMatics-Inc/imaginary-time-pairing-barriers). Related prior notes are [Suppression of electronic pair order for slow lattice phonons](https://doi.org/10.5281/zenodo.24366508), [Pair density matrices and axial helicity for quantum bond phonons on cylinders](https://doi.org/10.5281/zenodo.24369897), and [A two-thirds-power adiabatic bound for electronic pair order](https://doi.org/10.5281/zenodo.24372238).

## Scope and evidence

This is a paper of positive estimates and method-specific barriers. It is not a linear ceiling. For the spinful square-lattice Holstein and bond Su-Schrieffer-Heeger models with independent harmonic phonons, fixed bare hopping $t>0$, coupling $U>0$, and regular density, whether $k_B T_{\rm pair}\le C\Omega$ (or $\Omega$ times a fixed power of $\log(1/\Omega)$) holds remains OPEN. The paper states no physical family that violates such a ceiling and makes no claim about any material or measured transition temperature. Calculated values here are not measurements.

The replay scripts check finite exact-rational algebra and multi-precision interval spot points. They do not formalize the unbounded semigroup approximation, continuous-time dyadic chaining, the SU(2) representation decomposition, the Lehmann Hilbert-space limit, any infinite-volume argument, or any lattice pairing ceiling.

| Result | Label | Evidence |
|---|---|---|
| Imaginary-time reflection positivity, balanced and SU(2) zero-source sandwiches, interacting increment generating function and natural-amplitude slab budget, ordered-moment and auxiliary-CAR determinant bounds, gap-free frozen two-electron Combes-Thomas-type estimate, atomic connected source response | PROVED | `paper/main.tex`, `replay/linear5/RIGOR.md`, `replay/stage3/RIGOR.md` |
| Eight method barriers (unrestricted history, mirrored Wick growth, adiabatic tracking, canonical Holstein spatial reflection positivity, non-Gaussianity of averaged fermions, slow-field Hessian below $\Omega^{2/3}$, absolute-moment majorant, good-block temporal gluing) and the atomic Franck-Condon sideband weight | PROVED | `paper/main.tex` and the two `RIGOR.md` notes |
| Exact-rational finite identities and guards, with 24 mutant checks (`replay/linear5`) and 12 mutant checks (`replay/stage3`), 36 in total | CERTIFIED | `replay/linear5/replay.py --all-mutants`, `replay/stage3/replay.py --all-mutants` |
| Multi-precision interval spot points (for example the 60-digit finite-matrix Schatten $S_2$ spot check and the 70-digit interval points at the stated parameters) | NUMERICAL | `replay/linear5/replay.py`, `replay/stage3/replay.py` |
| A linear $\Omega$ ceiling for the pure models at fixed $t,U>0$ | OPEN | not claimed |

Interval points are finite spot checks at the stated parameters. They do not prove any statement about a parameter family.

## Review status

Sections 1 through 7 and 12 through 18 of the LINEAR5 RIGOR note received PASS in `REVIEW-RIGOR-1-18`; the STAGE3 note received PASS in `REVIEW-STAGE3`. The manuscript was reviewed in `REVIEW-PAPER` and a subsequent delta review. All these reviews were carried out by separate AI agents. Any "review pending" labels retained in the bundled notes predate those reviews. This preprint has not been peer reviewed.

## Contents and replay

- `paper/` contains the manuscript source and bibliography.
- `replay/linear5/` contains the replay (`replay.py`, 24 mutants), the companion proof note `RIGOR.md`, and the handoff `RESULT.md`.
- `replay/stage3/` contains the second replay (`replay.py`, 12 mutants), the companion proof note `RIGOR.md`, and the handoff `RESULT.md`.
- `CITATION.cff` gives citation metadata. `MANIFEST.sha256` pins the released source tree.

The two `replay.py` files, the two `RIGOR.md` notes, and the two `RESULT.md` files are retained byte for byte from the reviewed working copies. Their internal status lines (for example "independent review pending"), parent-directory links, and local command paths predate the independent reviews and are intentionally not edited here. The paper's Data availability paragraph records the review status of the underlying notes and of the typeset manuscript as published. Sections 8 to 11 of the older sharpness material are cited from the preceding preprint rather than reproduced.

With Python 3.12+, install `mpmath` and run both mutant-checking drivers from their respective directories:

```sh
python -m pip install mpmath==1.4.1
(cd replay/linear5 && python -B replay.py --all-mutants)
(cd replay/stage3 && python -B replay.py --all-mutants)
sha256sum -c MANIFEST.sha256
```

The manuscript can be built from `paper/` with Tectonic:

```sh
tectonic --keep-logs main.tex
```

The GitHub Actions workflow verifies the source manifest, runs both replay scripts with `--all-mutants`, checks prose punctuation, and builds the manuscript with Tectonic. It uploads the resulting PDF as a workflow artifact.

## Licenses and provenance

The replay code and documentation in this repository are licensed under Apache-2.0. The manuscript text (`paper/`) is distributed under CC-BY-4.0 through its Zenodo record. Not peer reviewed.

AI agents generated proofs, certificates, code, and manuscript drafts and edits under D. Rens' direction and review.
