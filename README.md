# A fixed-band seam obstruction for a cubed-sphere Laplacian

**Provisional research note. Independently checked; novelty not established.**

## Result in one paragraph

For a specific second-order, conservative, Galerkin-assembled Laplace–Beltrami operator on the gnomonic cubed sphere (shared panel-edge nodes, trapezoid mass, face-flux diagonal terms plus cell-centred cross terms), the truncation error on panel-edge ("seam") rows is exactly O(h) with a closed-form coefficient, while chart-interior rows are O(h²). The seam coefficient, in Fermi coordinates along the seam, is c₁ = A∂ₜ + B∂ₜₜ + ½∂ₙₙ + D∂ₜₜₜ with explicit A, B, D. Summing the residual over a band of fixed width k around the seams and taking h → 0 gives a k-independent functional ∫ f·G dt with an explicit nonzero density G. Consequently **no modification of the stiffness matrix supported on the band — however wide tangentially, provided it keeps either conservation (zero column sums) or symmetry, with bounded diagonal mass changes — can make the band residual o(h)**, with an explicit lower bound. A second, independent derivation shows the band functional equals the integrated h² truncation density of the *interior* stencil, which is why seam-local corrections cannot touch it.

## Status labels

| Aspect | Label |
|---|---|
| Mathematical validity | **Proved** (exact symbolic algebra; see PROOF.md). Independently checked by a second derivation route and a second model session. |
| Method | **Known** — standard SBP compatibility/telescoping argument (PRIOR_ART.md §1). |
| Phenomenon | **Anticipated** — diagonal-norm "interior 2s / boundary s" barrier (PRIOR_ART.md §2). |
| Specific result (this operator's exact seam coefficient; fixed-band obstruction in this correction class) | **Not found in examined sources; Shashkin–Goyman–Tolstykh (2023) full text examined on 2026-11-01: it documents cubed-sphere boundary-order reduction and associated accuracy effects but does not derive the restricted fixed-band obstruction; search not exhaustive; novelty not established** (PRIOR_ART.md §2–4). |
| Solution-accuracy implications | **None claimed** (LIMITATIONS.md). |
| Publication status | Research note only. Not peer reviewed. Not submitted anywhere. |

## Files

- `THEOREM.md` — self-contained definition of the operator, the correction class, the witnesses, and the precise statement with its constant.
- `PROOF.md` — Lemmas 1–5, the proof, and the independent interior-stencil derivation of G.
- `PRIOR_ART.md` — what is known, what anticipates the result, what was not found, and the unresolved comparisons.
- `LIMITATIONS.md` — scope, exclusions, load-bearing hypotheses, provenance.
- `verification/` — exact symbolic checks: `proof_session/` (from the derivation session) and `independent/` (from the review session, different routes); `run_all.sh` reproduces everything; `expected_output/` holds reference outputs.
- `provenance/` — SHA-256 manifest and the provenance record.

## Reproduce

```
pip install sympy mpmath numpy
bash verification/run_all.sh
```
About 3 minutes on a laptop: the proof-session scripts take ~1 minute together, `independent_seam.py` ~20 s, and `independent_boundary.py` is the longest (~1.5–2 minutes: its Part A derives the interior h² density symbolically before the boundary reduction). Every script uses exact arithmetic; the outputs should match `verification/expected_output/` apart from timings. Each script can also be run on its own with `python3 -I`.

## How to engage

The most useful outcomes would be: a reference showing this result (or a stronger one) is already in the literature; a counterexample within the stated class; an application of the interior/quadrature criterion of PROOF.md §6 to an established multiblock operator; or a result about solution rather than residual accuracy for this operator. Please open an issue.

## Provenance and AI disclosure

The derivation and the proof-session scripts were produced by an AI coding agent (Claude Opus 5.5) from the operator's design document, under the owner's direction, on 2026-10-31; the independent check and the second derivation were produced by a different model (Claude Fable 5.1) in a separate session; a third session (Claude, same day) re-derived the density by the interior route again, proved the conditional general form and the counterexamples, and ran the numerical pilot summarized in LIMITATIONS.md. The owner directed, reviewed and is responsible for the content. Details in LIMITATIONS.md and `provenance/`.

## License

Text (all Markdown files): CC BY 4.0 — see `LICENSE`. Code (`verification/`): MIT — see `verification/LICENSE`.

## Citation

Michael Blumhardt, "A fixed-band seam obstruction for a cubed-sphere Laplacian", research note, 2026. https://github.com/blumhardt-m/cubed-sphere-seam-obstruction, tag v1.0.1.
