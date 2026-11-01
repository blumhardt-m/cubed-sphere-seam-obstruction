# Limitations and scope

## What the theorem is

A statement about **local truncation error** on a fixed-width band of rows around the seams of **one specific discretization** (THEOREM.md §1), for **one correction class** (THEOREM.md §2), with **one explicit witness**. It is proved; the lower-bound constant is explicit but weak (≈ 0.0192/((1+C)(2k+1)) for u_A).

## What it is not

- **Not a statement about solution accuracy.** O(h) truncation confined to an O(h)-measure set routinely coexists with second-order solution convergence (supraconvergence; p+2 convergence for 2p/p SBP second-derivative operators). No claim about the convergence of discrete solutions, eigenvalues or wave propagation is made or implied. A separate numerical pilot (2026-10-31, owner's records) on an independent re-implementation of this operator from the design text found, for five smooth analytic wave solutions on n = 17–129, relative W-L2 and max-norm solution errors converging at observed order ≈ 2.00, with seam residual max order ≈ 1.02, seam-region residual W-L2 order ≈ 1.52 (the h^{3/2} scaling of the corollary) and O(1) vertex residuals; in a post hoc source decomposition the seam residual accounted for a substantial share of the finite error coefficient (component norm ≈ 58 % of the total at n = 129, partly cancelling the interior contribution) without reducing the observed order. These are finite-grid observations on a reconstruction, not theorems and not statements about the original implementation.
- **Not a general theorem about multiblock interfaces.** It concerns this operator. The mechanism (PROOF.md §6: the band functional equals the integrated interior truncation density) generalizes to a conditional statement — a fixed-band obstruction follows whenever the integrated interior truncation density plus quadrature error has a nonzero limit F — but F must be computed per scheme, and a universal obstruction for conservative multiblock schemes is false (uniform centered differences with labelled interfaces have F = 0).
- **Not a statement that a different operator is needed, or that any particular remedy works.** Changing the interior stencil, the mass quadrature outside the band, the DOF structure (e.g. duplicated interface nodes with penalty coupling), or allowing band width to grow with h are all outside the class and are not addressed.
- **Not a result about the cube vertices.** Vertex rows and rows whose stencils contain corner cells are not analysed; the witnesses are supported away from them.

## Hypotheses that are load-bearing

- The witness class U₀ requires Δu = ∂ₙΔu = 0 on the seam. This is what makes the mass-perturbation term O(h); without it the argument fails. The theorem therefore refutes o(h) residuals for *some* smooth u, which suffices for the obstruction, but the lower bound is only for the stated witnesses.
- Band corrections must have bounded relative diagonal mass change (|δW| ≤ C W) and must not touch rows outside B_k.
- The target refuted is o(h) on **all** band rows (max norm) — and, by the corollary, order better than h^{3/2} in global weighted RMS. A target imposed on seam rows only is not addressed.

## Exclusions, explicitly

Nonconservative non-symmetric D; non-diagonal or unbounded relative mass changes; supports growing with h; corrections outside the band; wider-stencil or interpolating formulations that change the global DOF structure; k = k(h) (the bound decays like 1/(2k+1) and gives no rate statement).

## Provenance and review limitations

- The derivation and the proof-session scripts were produced by an AI coding agent (Claude Opus 5.5) working from the operator's design text, in a single session on 2026-10-31, under the repository owner's direction. The independent check and the second derivation of G were produced by a different model (Claude Fable 5.1) in a separate session the same day, by routes that share no intermediate step with the proof session.
- Both the proof and its checks were produced for the same owner. No external numerical analyst has reviewed the result. That is the purpose of publishing it.
- The exact algebra is reproducible from `verification/run_all.sh` in a clean Python environment (sympy, mpmath, numpy); nothing depends on code outside this repository.
- Novelty is not established (PRIOR_ART.md §4).
- Numerical cross-checks against a direct implementation of the operator (seam residual regression coefficients at 15 seam points agreeing with the exact A, B, ½, D to ~1e-6; band functionals at n = 513 within 0.4 % of the exact limits) exist in the owner's records and were consistent, but are not part of this release and are not proof steps.
- **2026-11-01 analytic audit (AI-assisted).** The four non-symbolic steps of the proof (mass-perturbation estimate, pigeonhole constants and pair quantification, weighted band-norm corollary, cutoff geometry) were rederived by an AI agent in a fresh session from a blind packet containing the operator, correction class, witness data and the §4 identities, with the theorem's bounds withheld and `THEOREM.md`/`PROOF.md` not accessible until the derivation was hashed and frozen (Stage A). A second session then compared the frozen derivation against the published text (Stage B). Result: the fixed-witness bounds (i), (ii) and the corollary agree; the published band-mass factor 36 is valid but looser than the 24 obtained independently; adversarial examples demonstrate the failure of the corresponding unrestricted inverse-mass, signed-weight-norm, and single-witness symmetric assertions, which the published argument does not make. The comparison found that the general witness construction in §3 was stated for |t| < 0.6155, which is not uniformly admissible with the fixed cutoff; v1.0.1 restricts it to |t| < 0.45, a verified sufficient condition. Independence limits: the blind packet was prepared from theorem material, so both stages share its assumptions; Stage A's access controls are self-attested; the §4 symbolic identities were assumed by both stages and remain machine-checked (`verification/`) but not rederived; no human verification has occurred. Status: fixed-witness analytic bounds corroborated by blind AI-assisted reconstruction and post-freeze comparison; symbolic identities assumed; human verification outstanding.

## Relationship to the discretization's development

The operator analysed here was a candidate component of a small wave-simulation instrument. This repository is a mathematical analysis of that candidate; it is not a validated implementation, it is not the instrument, and it carries no qualification or performance claim about either.
