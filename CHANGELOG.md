# Changelog

## v1.0.1 — 2026-11-01

v1.0.1 corrects an overbroad auxiliary witness-construction statement. The theorem's two fixed-witness obstruction conclusions, their supporting proof arguments, and their quantitative constants remain unchanged.

- THEOREM.md §3: the general construction u_f was stated for f ∈ C_c^∞(|t| < 0.6155). With the fixed normal cutoff (ρ = 0 for |n| ≥ 0.20) that range is not uniformly admissible: for f supported near the ends of the seam segment the constructed function leaves Ω. Restricted to |t| < 0.45, a verified sufficient condition for containment. The two witnesses actually used (u_A, support |t| < 0.42; v, support |t| < 0.45) were and remain admissible (PROOF.md Lemma 5).
- PROOF.md §7 step 1: O(kh) → O((k+1)h) so the k = 0 case is not vacuous.
- LIMITATIONS.md: records the 2026-11-01 AI-assisted blind rederivation and publication comparison and their independence limits.
- PRIOR_ART.md: Shashkin, Goyman & Tolstykh (2023) Part I examined in full; moved from unresolved (§4) to related prior work (§2); novelty-status paragraph updated accordingly.
- CITATION.cff: version 1.0.1.

Found by: Stage A blind analytic rederivation (fresh-context AI session, theorem bounds withheld) followed by Stage B comparison of the frozen derivation against the published text. The construction-scope defect was detectable only by the comparison step.

## v1.0.0 — 2026-10-31

Initial release.
