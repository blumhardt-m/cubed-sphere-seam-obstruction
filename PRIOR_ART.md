# Prior art and novelty status

**Novelty status: NOT ESTABLISHED.** The method of proof is known. The phenomenon is anticipated qualitatively. The specific result — the exact seam coefficient for this operator and the fixed-band obstruction in the class 𝒞_B(k) — was not found in any examined source, but one directly relevant text (Shashkin, Goyman & Tolstykh 2023, Part I) could not be examined and the search was not exhaustive. Readers who know of prior work covering this result are asked to open an issue.

## 1. What is known and is used here

**Compatibility / telescoping identity (Lemma 1).** Pairing an algebraic property of the operator (zero column sums, or symmetry) with the residual summed against a test vector is the standard summation-by-parts (SBP) compatibility argument. The conservation form 1ᵀ(·) is the classical telescoping proof of discrete conservation (Nordström & Carpenter, ICASE Report 98-19, Theorem 5; Carpenter, Nordström & Gottlieb, J. Comput. Phys. 148, 1999). A modern general statement of the exactness/compatibility duality is Glaubitz, Nordström & Öffner, "Summation-by-parts operators for general function spaces", SIAM J. Numer. Anal. 61(2), 2023, Lemma 4.3 / Theorem 4.4 (arXiv:2203.05479). None of these concern second-derivative operators on multiblock spheres or seam truncation coefficients; the theorem here neither uses nor yields their quadrature-existence statements.

**Interior second-order consistency (Lemma 2)** of a conservative face-flux scheme with arithmetically averaged coefficients, plus a cell-centred cross term with averaged corner differences, is routine.

## 2. What anticipates the phenomenon

**The diagonal-norm order barrier.** For SBP first-derivative operators with a diagonal norm, interior order 2s forces boundary-closure order ≤ s: Kreiss & Scherer (1974); Strand (1994); Linders, Lundquist & Nordström, "On the order of accuracy of finite difference operators on diagonal norm based summation-by-parts form", SIAM J. Numer. Anal. 56(2), 2018, doi:10.1137/17M1139333 — Theorems 4–5 (uniform grid) and 9–10 (finitely many arbitrarily placed boundary nodes): s ≥ τ, and for s = τ the diagonal norm is a quadrature of order exactly 2τ. The proof runs through the norm-as-quadrature identity, a close cousin of the interior/quadrature identity in PROOF.md §6. Related: Lundquist & Nordström (2015 preprint, Proposition 1) bound the exactness of adjoint-compatible interpolation between two SBP norms by an explicit moment mismatch — the role a nonzero F plays here. For second-derivative SBP operators the s-order boundary closure is presented as a construction (Mattsson & Nordström, J. Comput. Phys. 199, 2004; Wang & Kreiss, arXiv:1503.08926; Wang & Petersson, arXiv:1809.04310), with global convergence of order s+2 recovered by normal-mode or energy analysis. Interface rows of SBP-SAT couplings carry the same O(h^s) truncation as boundary rows (Wang & Kreiss, §3).

The (S) variant of the theorem is the interface analogue of this barrier for a particular second-derivative operator. The (C) variant — conservation only, no symmetry — is a strictly weaker hypothesis not covered by these results. Neither variant follows from them: the SBP theorems are one-dimensional, for operators of the form P⁻¹Q, at physical boundaries.

**Cubed-sphere edge accuracy.** Ullrich, Jablonowski & Lauritzen, J. Comput. Phys. 229 (2010) observe empirically that the panel edges control the global convergence order and give one-sided stencil error terms at edges (first order in the edge-normal direction) in their Appendix A; Katta, Nair & Kumar, Mon. Wea. Rev. 143 (2015) and Thuburn, Cotter & Dubos, Geosci. Model Dev. 7 (2014) report edge/corner-dominated error without truncation analysis. Shashkin, Goyman & Tretyak (arXiv:2501.07323, 2025) state for SBP shallow-water models on cubed-sphere grids that "2s-order accurate operators at the inner points can simultaneously have only s-order accuracy near the edges", citing the diagonal-norm barrier, and perform no truncation analysis of interface rows.

**O(h) truncation with O(h²) convergence (supraconvergence).** Kreiss, Manteuffel, Swartz, Wendroff & White, Math. Comp. 47 (1986); Manteuffel & White, Math. Comp. 47 (1986): conservative schemes on nonuniform grids yield second-order solutions despite first-order truncation. This is the reason the theorem here does not bear on solution accuracy (see LIMITATIONS.md).

**Discrete-Laplacian no-go results.** Wardetzky, Mathur, Kälberer & Grinspun, "Discrete Laplace operators: no free lunch", SGP 2007, is qualitative (no rates, no interfaces).

## 3. What was not found

No examined source contains:
- a closed-form seam truncation coefficient (c₁, or its Fermi form A, B, ½, D) for a Galerkin-assembled Laplace–Beltrami operator on the gnomonic cubed sphere;
- a fixed-band non-improvability theorem for a second-derivative interface in a class permitting arbitrary tangential coupling under conservation alone or symmetry alone;
- the observation that the band functional is determined by the interior truncation density (PROOF.md §6) stated as such, together with the conditional general form, the 1-D jump-spacing example and the uniform-grid counterexample that bound its reach.

## 4. Unresolved comparisons

- **Shashkin, Goyman & Tolstykh, "Summation-by-parts finite-difference shallow water model on the cubed-sphere grid. Part I: Non-staggered grid", J. Comput. Phys. 474 (2023) 111797.** Full text not accessed (publisher paywall; an SSRN preprint, abstract id 4176834, exists but was not retrieved; no arXiv version). Its interface treatment and any interface accuracy analysis are **not assessed**.
- The older multiblock-interface SBP literature and the Russian-language INM RAS literature were not surveyed.
- Search effort: three bounded searches over public sources only (the proof session; the review session; a subsequent diagnostic campaign that retrieved the Linders–Lundquist–Nordström and Lundquist–Nordström full texts from DiVA).
