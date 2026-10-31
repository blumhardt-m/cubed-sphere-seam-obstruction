# Proof

Notation and the statement are in `THEOREM.md`. Every step below is exact algebra, elementary analysis, or exact symbolic algebra in sympy (rational and surd arithmetic; no floating point enters any proof step). Scripts are in `verification/`; §8 maps each lemma to its script. Numerical values appear only as remarks and are not used.

## 1. Lemma 1 — exact band identities

From −Ku = W(r + Δu) and −K′u = W′(r′ + Δu): Du = W r − δW Δu − W′ r′. Since every nonzero row of D lies in B_k:

- **(C)** Σ_{B_k} (Du)_j = 1ᵀDu = 0, hence Σ_{B_k} W r(u) = Σ_{B_k} δW Δu + Σ_{B_k} W′ r′(u).
- **(S)** vᵀDu = Σ_{B_k} v_j (Du)_j and uᵀDv = Σ_{B_k} u_j (Dv)_j; symmetry gives equality, hence
  Σ_{B_k} W[v r(u) − u r(v)] = Σ_{B_k} δW (vΔu − uΔv) + Σ_{B_k} W′[v r′(u) − u r′(v)].

No locality of D is used: the sums run over all of B_k and the columns of D are unrestricted. ∎

## 2. Lemma 2 — chart-interior consistency (band rows 1…k)

Take a node whose full stencil lies in one chart (four full-weight faces, four cells, W = √g h²) and any smooth U on a neighbourhood of the closed chart. The residual R_int(ξ₀, η₀, h) at that node is a smooth function of (ξ₀, η₀, h) including at h = 0, and its h⁰ and h¹ Taylor coefficients vanish identically (`d3_interior_series.py`; independently `independent_boundary.py`, Part A). A generic 4-jet suffices because the h^j coefficient depends only on the j-jet of U (node displacements are O(h)). Hence |r_j| ≤ c_u h² uniformly on compact sets.

PX rows n−1−i and PY rows i, 1 ≤ i ≤ k, use exactly this stencil: their stencils reach the seam row of their own chart, where U = u∘chart is smooth, and all weights are full. ∎

The h² coefficient r₂[U] of the interior residual is a differential expression of order ≤ 4 in U with coefficients rational in (ξ, η); its closed form is printed by `independent_boundary.py` and is used in §6.

## 3. Lemma 3 — exact seam-row expansion

At a seam node PX (1, e) ≡ PY (−1, e) with interior η index, differentiating the assembled energy with respect to the seam dof gives

  (Ku)_b = Σ_{j=±1} C_j (u_b − u(1, e+jh))
         + Σ_{side∈{PX,PY}} [ c (u_b − u_side(1,0)) − (κ₊/2)(u_b − u_side(1,+1)) + (κ₋/2)(u_b − u_side(1,−1)) ],

with C_j = ¼(√g(1,e)+√g(1,e+jh))(g²²(1,e)+g²²(1,e+jh)), c = ¼(√g(1,e)+√g(1−h,e))(g¹¹(1,e)+g¹¹(1−h,e)), κ_± = κ(1−h/2, e±h/2), W_b = √g(1,e) h². Here u_side(1, j) is the row-1 node of that chart at η index b+j; PY's coefficients equal PX's because √g, g¹¹, g²² are even in ξ and κ is odd in ξ, with the cell-corner roles mirrored.

**Mirror symmetry.** The PX node (1−h, e+jh) and the PY node (−1+h, e+jh) are the points (1, 1−h, e+jh)/ρ and (1−h, 1, e+jh)/ρ: mirror images under x↔y with identical coefficients; seam nodes lie on the mirror plane. For mirror-odd u, Δu = 0 on the plane and r_b = 0 exactly for every h.

**Series.** Let U = u∘PX, extended smoothly past ξ = 1 (valid for x > 0). A PY node (−1+h, e+jh) has PX coordinates (1/(1−h), (e+jh)/(1−h)). Expanding for a generic 4-jet of U at (1, e):

- (Ku)_b has zero h⁰ and h¹ coefficients;
- r_b has zero h⁰ coefficient (consistency: −(Ku)_b/W_b → ΔU);
- the h¹ coefficient is exactly
  c₁ = (e²+2)U_ξ + (e²+2)U_ξξ + e(e²+2)U_ξη + e(4e⁴+21e²+41)/(4(e²+2)) U_η − (e²/2+1)U_ηη − e(e²+2)²/6 U_ηηη.

**Fermi form.** Substituting U = F(n(ξ,η), t(ξ,η)) and applying the chain rule exactly, c₁ = A F_t + B F_tt + ½F_nn + D F_ttt with A, B, D as in THEOREM.md §4 and no F_n, F_nt, F_ntt, F_nnn or F_nnt terms. As a check, Δ at the seam comes out as F_nn + F_tt.

**Uniformity.** R_seam(e, h) := −(Ku)_b/W_b − ΔU is smooth in (e, h) for e in a compact subset of (−1, 1) and h ∈ [0, h₁]: (Ku)_b is h² times a smooth function (Hadamard's lemma, since its h⁰ and h¹ coefficients vanish) and W_b/h² = √g(1,e) > 0. With R_seam(e, 0) = 0 and ∂_h R_seam(e, 0) = c₁(e):

  |r_b − h c₁(e_b)| ≤ c_u h²  uniformly on supp u. ∎

(`d1_seam_series.py`, `d2_fermi_reduction.py`; independently `independent_seam.py`, which builds the row by symbolic differentiation of the assembled energy with PY's own metric and confirms the Fermi form additionally on explicit polynomial test functions.)

## 4. Lemma 4 — limits of the band functionals

For u, v ∈ U₀ define E_k(u) = h⁻² Σ_{B_k} W r(u) and S_k(u, v) = h⁻² Σ_{B_k} W[v r(u) − u r(v)]. For small h only nodes whose stencils meet supp u contribute: the PX|PY seam row and PX/PY rows 1…k over a compact η range (Lemma 5 places the other seams and the vertices at positive distance).

- Rows 1…k: at most 2k·O(1/h) nodes, each with W ≤ h² and |r| ≤ c h² (Lemma 2): contribution O(kh) → 0.
- Seam row: h⁻² Σ_b √g(1,e_b) h² (h c₁(e_b) + O(h²)) = Σ_b h √g(1,e_b) c₁(e_b) + O(h), a Riemann sum of a smooth compactly supported function.

Therefore, for every fixed k,

  E_k(u) → E_∞(u) = ∫ √g(1,e) c₁[u](e) de = ∫ ω(t) c₁[u] dt,  ω = √g · de/dt = (2(2+e²))^{−1/2},
  S_k(u, v) → S_∞(u, v) = ∫ ω (v c₁[u] − u c₁[v]) dt,

independent of k. On U₀ the conditions Δu = ∂ₙΔu = 0 at n = 0 (and their t-derivatives) give F_nn = −f″, F_nnt = −f‴, F_nnn = F_n − F_ntt with f = u|_seam; since c₁ contains no odd-in-n terms only the first is needed, and c₁[u] = L₀f with

  L₀ = A ∂ₜ + (B − ½) ∂ₜ² + D ∂ₜ³,  B − ½ = (e² − 2)/(e² + 2).

Integrating by parts (f compactly supported):

  E_∞(u) = ∫ f(t) G(t) dt,  G = −(ωA)′ + (ω(B−½))″ − (ωD)‴ = √2 (3e⁴ + 31e² − 25) / (4 (e²+2)^{5/2}),

exact, with G(0) = −25/16. ∎

The symmetric-class statement follows. For L = a∂ + b∂² + d∂³ with smooth coefficients, L* = −∂(a·) + ∂²(b·) − ∂³(d·), so L − L* has third-order coefficient 2d and first-order coefficient 2a − 2b′ + 3d″. Here d = ωD = −e/(3(e²+2)^{3/2}) ≠ 0 for e ≠ 0, which already settles non-self-adjointness; the first-order defect is ωA − (ω(B−½))′ + (3/2)(ωD)″ = e(29e²+103)/(12(e²+2)^{5/2}), also nonzero for e ≠ 0. (The quantity ωA − (ω(B−½))′ = e(6e⁴−e²+19)/(12(e²+2)^{5/2}) that appears in the original derivation is the first-order defect only when d = 0; it is not the relevant quantity here.)

## 5. Lemma 5 — explicit nonzero witnesses

- **Support.** For |t| ≤ t_m and |n| ≤ 0.20: min(x, y) = (cos n cos t − |sin n|)/√2 ≥ (cos 0.2 cos t_m − sin 0.2)/√2 and |z| = cos n |sin t| ≤ sin t_m. For t_m = 0.42 the bounds are 0.4923 > 0.4078; for t_m = 0.45 they are 0.4835 > 0.4350 (interval arithmetic in `d5_witness_bounds.py`). So supp u_A and supp v lie in Ω at positive distance from the other eleven seams and the eight vertices.
- **Conservation witness.** On |t| ≤ 0.42, tan t ≤ t/(1 − t²/2), so e² = 2tan²t ≤ 0.42436, while G < 0 iff e² < (−31 + √1261)/6 ≈ 0.7518. The polynomial 3e⁴ + 31e² − 25 is increasing in e² ≥ 0 and equals −11.30 at the bound, so G < 0 on supp ψ_A, and E_∞(u_A) = ∫ ψ_A G dt < 0 strictly (ψ_A ≥ 0, ψ_A ≢ 0).
- **Symmetric witness.** v = u_g with g ≡ 1 on |t| ≤ 0.42 ⊇ supp ψ_A. Where ψ_A ≠ 0, v ≡ 1 near the seam and c₁[v] = L₀g = 0 (L₀ has no zeroth-order term); also g = 1 on supp ψ_A. Hence S_∞(u_A, v) = ∫ ω g L₀ψ_A dt − ∫ ω ψ_A L₀g dt = E_∞(u_A) < 0. ∎

## 6. Second derivation of G from the interior stencil (independent route)

This section is not needed for the theorem; it is an independent confirmation of Lemmas 3–4 that never uses the seam row, and it explains the mechanism.

Because K·1 = 0 exactly, Σ_all W r(u) = −1ᵀKu − Σ_all WΔu = −Q_W(Δu), where Q_W is the six-chart product-trapezoid quadrature. For u ∈ U₀, √g Δu∘chart vanishes to second order at the seam (ξ = 1 on PX), so its ξ-derivative vanishes there and the Euler–Maclaurin h² term of the trapezoid rule is zero: Q_W(Δu) = O(h⁴). The complement of B_k within supp u consists of chart-interior nodes, so by Lemma 2

  E_∞(u) = −∫_{PX} √g r₂[U] dξ dη − ∫_{PY} √g r₂[U] dξ dη.

`independent_boundary.py` carries this out exactly:

- the Euler–Lagrange expression Σ(−1)^{p+q}∂^{pq}a_pq of √g·r₂ = Σ a_pq ∂^{pq}U vanishes identically — √g·r₂ is an exact divergence, as conservation of the interior scheme requires;
- integrating by parts leaves a seam boundary functional Σ_{m≤3} ∫ b_m(η) ∂_ξ^m U(1,η) dη with explicit b_m;
- converting to Fermi jets, imposing the U₀ conditions and integrating by parts in t gives, for the PX chart, the density −√2(3e⁴+31e²−25)/(8(e²+2)^{5/2}) against f, together with odd-in-n jet terms;
- the PY chart is the x↔y mirror image of PX and the design is mirror-symmetric, so its density has the odd-in-n terms with opposite sign and the same even part; the sum is exactly G.

Interpretation: the band functional E_∞ is determined entirely by the interior discretization and the mass quadrature. This is why no correction supported in a fixed band can change it, and it identifies the result as an instance of the conservation/quadrature order barrier rather than as a defect of any particular seam closure.

**Conditional general form and its limits.** The argument of §7 does not depend on the cubed sphere. Abstractly: for any sampled discretization L_h = −W_h⁻¹K_h of a differential operator A on a closed manifold, a node set B_h with Σ_{B_h} W ≤ M_B h^β, corrections with rows in B_h and |δW| ≤ C W: (conservative version) if 1ᵀD = 0 and a fixed witness u has h^{−p} Σ_{B_h} W r_h(u) → F ≠ 0 and h^{−p} Σ_{B_h} W |Au| → 0, then liminf h^{−(p−β)} ‖r′_h(u)‖_{∞,B_h} ≥ |F|/((1+C)M_B); (symmetric version) if D = Dᵀ and a fixed pair u, v has h^{−p} Σ_{B_h} W (v r_h(u) − u r_h(v)) → F_S ≠ 0 and h^{−p} Σ_{B_h} W |v Au − u Av| → 0, then liminf h^{−(p−β)} [V ‖r′_h(u)‖_{∞,B_h} + U ‖r′_h(v)‖_{∞,B_h}] ≥ |F_S|/((1+C)M_B), with U = sup|u|, V = sup|v|. The symmetric version is irreducibly a pair statement (remark after §7). When 1ᵀK_h = 0 and ∫Au = 0, F = −q_p(Au) − ∫T_p u, where q_p is the h^p quadrature-error coefficient of W_h and T_p the interior truncation density — the identity used above. This conditional statement, with explicit examples, was established in a separate diagnostic campaign on 2026-10-31 (owner's records), which also showed:

- a 1-D conservative two-block scheme with a spacing jump (spacings h and 2h, shared interface node, W_j = (d₋+d₊)/2) has F = 3 for the witness (x−1)³ near the interface, so the obstruction is not cubed-sphere-specific; its tensor lift gives a 2-D example;
- uniform centered differences with arbitrary labelled "interfaces" have F = 0 and are uniformly second order, so **a universal obstruction for conservative multiblock schemes is false**; nonzero F must be computed per scheme;
- the mass-negligibility hypothesis cannot be dropped: in the 1-D example with k ≥ 1, a diagonal mass change of relative size 1/8 at one node, combined with a conservative rank-one D, makes the band residual of that witness vanish identically. In the cubed-sphere theorem this hypothesis is supplied by the witness condition Δu = O(n²).

Vanishing of F is necessary for a fixed-band correction to reach o(h); it is not sufficient for the existence of an accurate one.

## 7. Proof of the theorem

1. **Mass term.** Band nodes in supp u lie within normal distance n ≤ kh/√2 of the seam, inside the ρ-plateau for small h, where |Δu| ≤ c n² ≤ c′(k+1)²h²; Δu = 0 outside supp u. With |δW| ≤ C W and Σ_{B_k ∩ supp} W = O(kh), the term h⁻² Σ δW Δu is O(h).
2. **Band mass bound.** √g ≤ 1 and alias weights sum to at most 1, so W_j ≤ h². B_k has at most 12(2k+1)n nodes, and n = 2/h + 1 ≤ 3/h for h ≤ 1. Hence Σ_{B_k} |W′| ≤ (1+C)·12(2k+1)·(3/h)·h² = 36(1+C)(2k+1) h.
3. **Variant (C).** By Lemma 1 (C) and step 1, |E_k(u_A) − O(h)| = h⁻²|Σ_{B_k} W′ r′(u_A)| ≤ h⁻² max_{B_k}|r′(u_A)| · 36(1+C)(2k+1) h. Lemma 4 gives E_k(u_A) → E_∞(u_A) ≠ 0 (Lemma 5), which is (i).
4. **Variant (S).** The same with Lemma 1 (S): |Σ W′[v r′(u_A) − u_A r′(v)]| ≤ M (max|r′(u_A)| + max|r′(v)|) Σ|W′|; the mass term is O(h) as in step 1 with |u_A|, |v| ≤ M; and S_k(u_A, v) → E_∞(u_A) ≠ 0, which is (ii).
5. **Corollary (C).** From step 3, |Σ_{B_k} W′ r′(u_A)| ≥ (|E_∞|/2) h² for small h; by Cauchy–Schwarz, (Σ_{B_k}|W′ r′|)² ≤ (Σ_{B_k}|W′|)(Σ_{B_k}|W′| r′²), and step 2 bounds the first factor by 36(1+C)(2k+1)h.
6. **Corollary (S).** From step 4, |Σ_{B_k} W′[v r′(u_A) − u_A r′(v)]| ≥ (|E_∞|/2) h² for small h. By Cauchy–Schwarz and the triangle inequality in the |W′|-weighted norm, the left side is at most (Σ_{B_k}|W′|)^{1/2} · (V ‖r′(u_A)‖_{|W′|,B} + U ‖r′(v)‖_{|W′|,B}), with U = sup|u_A|, V = sup|v|; step 2 bounds the first factor. The bound is on the pair and cannot be assigned to u_A alone, by the next remark. ∎

**Remark (single-witness symmetric cancellation).** Let u be any witness, δW = 0, and b := (W r(u))|_{B_k}, extended by zero. If bᵀu ≠ 0, set D = b bᵀ/(bᵀu). Then D = Dᵀ, every nonzero row and column of D lies in B_k, and Du = b. The identity Du = W r(u) − W′ r′(u) gives W r′(u) = 0 on B_k, i.e. r′(u) ≡ 0 on the band. Hence no single-witness lower bound can hold in variant (S): the theorem's variant (ii) and Corollary (S) are necessarily statements about a pair. Such a D violates 1ᵀD = 0, since 1ᵀb = Σ_{B_k} W r(u) = h²E_k(u) ≠ 0 for u ∈ U₀ with E_∞(u) ≠ 0; this is why variant (i) needs only one witness. (For the non-generic case bᵀu = 0 a rank-two symmetric D does the same.) This remark is a statement about one fixed witness; it does not produce a correction that is accurate for all smooth functions.

## 8. Scripts

| Script | Establishes | Origin |
|---|---|---|
| `proof_session/d1_seam_series.py` | Lemma 3 series, PX basis | proof session |
| `proof_session/d2_fermi_reduction.py` | Fermi form, U₀ reduction, ω, G, symmetric densities | proof session |
| `proof_session/d3_interior_series.py` | Lemma 2 | proof session |
| `proof_session/d5_witness_bounds.py` | Lemma 5 bounds (interval arithmetic) | proof session |
| `independent/independent_seam.py` | Lemma 3 by energy differentiation; Fermi form by chain rule and by explicit test functions; G; symmetric densities | review session |
| `independent/independent_boundary.py` | Lemma 2 with the closed-form density r₂; §6 (Euler–Lagrange check, boundary reduction, G) | review session |

All scripts run with `python3 -I` and import nothing beyond sympy, mpmath and numpy. `verification/run_all.sh` runs them all; `verification/expected_output/` holds the outputs they are expected to reproduce.
