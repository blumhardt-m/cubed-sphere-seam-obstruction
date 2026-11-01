# The restricted seam obstruction theorem

This file states the result precisely and self-containedly. The proof is in `PROOF.md`; what the result does not say is in `LIMITATIONS.md`.

## 1. The operator

**Surface and charts.** The unit sphere S², covered by six gnomonic charts. Two of them are used here:

- PX: (ξ, η) ↦ (1, ξ, η)/ρ,  PY: (ξ, η) ↦ (−ξ, 1, η)/ρ,  with ρ = √δ, δ = 1 + ξ² + η², (ξ, η) ∈ [−1, 1]².

Induced metric (identical on every chart): √g = δ^{−3/2}, g¹¹ = (1+ξ²)δ, g²² = (1+η²)δ, g¹² = ξηδ. Write κ = √g·g¹² = ξη δ^{−1/2}.

**Mesh.** n nodes per chart edge, spacing h = 2/(n−1), nodes at ξ_a = −1 + ah, η_b = −1 + bh. Nodes on a shared chart edge are a single global degree of freedom (one entry of the state vector u ∈ ℝ^N). The PX|PY seam identifies PX (a = n−1, b) with PY (a = 0, b), same b, no reversal.

**Mass.** W is diagonal: W_g = Σ_{aliases (F,a,b) of g} w_F(a,b) √g_F(a,b) h², with the product trapezoid weight w_F(a,b) = wx(a)·wx(b), wx(0) = wx(n−1) = ½, otherwise 1. Hence W_g = √g h² at an interior node and at an ordinary seam node (½ + ½), and ¾√g h² at a cube vertex.

**Stiffness.** K = Σ_F S_Fᵀ K_F S_F (Galerkin assembly over charts, S_F the 0/1 restriction), where K_F is the Hessian of the chart energy E_F = Σ(faces) + Σ(cells):

- x-face between (a,b) and (a+1,b): E = ½ · wx(b) · mean(√g) · mean(g¹¹) · (u_{a+1,b} − u_{a,b})²;
- y-face between (a,b) and (a,b+1): E = ½ · wx(a) · mean(√g) · mean(g²²) · (u_{a,b+1} − u_{a,b})²;
  (mean = arithmetic mean of the nodal values at the face's two nodes; faces along a chart boundary carry that chart's half weight, and the neighbouring chart supplies the other half);
- cell with corners u₀₀, u₁₀, u₀₁, u₁₁: E = κ_c · P_c · Q_c · h², with κ_c evaluated analytically at the cell centre, P_c = ½[(u₁₀−u₀₀)+(u₁₁−u₀₁)]/h, Q_c = ½[(u₀₁−u₀₀)+(u₁₁−u₁₀)]/h.

Every face and cell is entered exactly once. K is symmetric with K·1 = 0.

**Operator and residual.** L = −W⁻¹K approximates the Laplace–Beltrami operator Δ. For smooth u sampled at the nodes, r(u) = Lu − Δu = −Ku/W − Δu.

## 2. The correction class 𝒞_B(k)

Fix k ≥ 0, independent of h. The band B_k is the set of global dofs having an alias (F, a, b) with min(a, b, n−1−a, n−1−b) ≤ k: the seam rows and the k rows on either side, around all twelve seams, including the vertices.

A correction is K′ = K + D, W′ = W + δW such that

- every nonzero row of D lies in B_k (columns unrestricted — arbitrary tangential extent, any sign, no positivity, interpolation weights allowed, arbitrary dependence on h);
- **either (C)** 1ᵀD = 0 (the corrected operator stays conservative: 1ᵀK′ = 0), **or (S)** D = Dᵀ (W′L′ stays symmetric);
- δW is diagonal, supported in B_k, with |δW_j| ≤ C·W_j for a constant C independent of h, and W′_j ≠ 0.

The state remains the same single global vector u, and r′(u) = −K′u/W′ − Δu.

## 3. Seam coordinates and the witnesses

The PX|PY seam is the great circle x = y. Fermi coordinates about it: ν = (1,−1,0)/√2, n = arcsin(p·ν) (signed geodesic distance from the seam), t = atan2(p·ẑ, p·(1,1,0)/√2) (arclength along it). On the seam η = √2 tan t, and the open seam segment is |t| < atan(1/√2) ≈ 0.6155. In these coordinates Δ = ∂ₙₙ − tan n ∂ₙ + sec²n ∂ₜₜ.

**Witness class U₀:** smooth u on S² with compact support inside Ω = {x > 0, y > 0, |z| < min(x,y)} (the PX and PY charts across their common open seam, away from all other seams and vertices) such that Δu = 0 and ∂ₙΔu = 0 on the seam.

**Construction:** for f ∈ C_c^∞(|t| < 0.45) and a C^∞ cutoff ρ(n) equal to 1 for |n| ≤ 0.12 and 0 for |n| ≥ 0.20,

  u_f = [f(t) − ½ n² f″(t)] ρ(n),   so that Δu_f = n²(2f″ − ½f⁗) + O(n⁴) on the plateau |n| ≤ 0.12.

With this fixed normal cutoff, the tangential restriction |t| < 0.45 is a verified sufficient condition for supp u_f ⊂ Ω (PROOF.md, Lemma 5): Ω narrows to zero normal width at the two ends of the seam segment (the cube vertices), and the strip |n| ≤ 0.20 is no longer contained in Ω by |t| = 0.5. For f supported closer to an endpoint the cutoff width would have to shrink accordingly; no such f is used here. (v1.0.0 stated the construction for |t| < 0.6155, which is not uniformly admissible with this cutoff; see CHANGELOG.md.)

**Specific witnesses:**
- u_A = u_f with f = ψ_A(t) = β(t/0.42), β(s) = exp(1 − 1/(1−s²)) for |s| < 1 (zero otherwise);
- v = u_g with g ∈ C_c^∞(|t| < 0.45), g ≡ 1 on |t| ≤ 0.42.

Both supports lie in Ω at positive distance from its boundary (PROOF.md, Lemma 5).

## 4. The limiting seam quantities

Define (PROOF.md, Lemmas 3–4):

- Seam truncation coefficient. For any smooth u and seam node at parameter e = η, the seam-row residual is r = h·c₁(e) + O(h²) with, in Fermi form,
    c₁ = A·F_t + B·F_tt + ½·F_nn + D·F_ttt,
    A = √2 e(47e²+139)/(12(e²+2)²),  B = (3e²−2)/(2(e²+2)),  D = −√2 e/(3(e²+2)),  e = √2 tan t,
  where F(n,t) = u in Fermi coordinates. No odd-in-n derivatives occur.
- Conservation density. For u ∈ U₀ with f = u|_seam,
    E_∞(u) := lim_{h→0} h⁻² Σ_{j∈B_k} W_j r_j(u) = ∫ f(t) G(t) dt,   G = √2(3e⁴ + 31e² − 25)/(4(e²+2)^{5/2}).
  The limit does not depend on k. G(0) = −25/16, and G < 0 for |t| < 0.54999.

## 5. Theorem

> Assume the operator of §1, h → 0 along any sequence, and fix k ≥ 0.
>
> **(i) Conservative corrections.** For the witness u_A, every sequence of corrections in 𝒞_B(k) satisfying 1ᵀD = 0 obeys
>
>   liminf_{h→0} h⁻¹ · max_{j∈B_k} |r′_j(u_A)|  ≥  |E_∞(u_A)| / (36(1+C)(2k+1))  > 0.
>
> **(ii) Symmetric corrections.** For the pair (u_A, v), every sequence of corrections in 𝒞_B(k) satisfying D = Dᵀ obeys
>
>   liminf_{h→0} h⁻¹ · ( max_{B_k}|r′(u_A)| + max_{B_k}|r′(v)| )  ≥  |E_∞(u_A)| / (36(1+C)(2k+1)·M)  > 0,
>
> with M = max(sup|u_A|, sup|v|).
>
> In particular no correction in either class achieves band residuals that are o(h), for any fixed k.

E_∞(u_A) < 0 strictly, by a sign argument (G < 0 on supp ψ_A, ψ_A ≥ 0, ψ_A ≠ 0). Numerically E_∞(u_A) ≈ −0.690415, so the constant in (i) is ≈ 0.0192/((1+C)(2k+1)).

**Corollary (residual norms).** Write ‖r‖_{|W′|,B} = (Σ_{j∈B_k} |W′_j| r_j²)^{1/2}; the absolute values are needed because the class only requires W′_j ≠ 0 (when W′ > 0 this is the ordinary W′-weighted L2 norm on the band). From the same inequalities and Cauchy–Schwarz, for small h:

- **(C) conservative corrections:** ‖r′(u_A)‖_{|W′|,B} ≥ (|E_∞(u_A)|/2) h^{3/2} / (36(1+C)(2k+1))^{1/2}. So for the witness u_A no conservative correction achieves a band- or global-weighted-RMS residual of order better than h^{3/2}, nor a band-relative weighted RMS (dividing by Σ_{B_k}|W′| = O(h)) of order better than h.
- **(S) symmetric corrections:** V·‖r′(u_A)‖_{|W′|,B} + U·‖r′(v)‖_{|W′|,B} ≥ (|E_∞(u_A)|/2) h^{3/2} / (36(1+C)(2k+1))^{1/2}, with U = sup|u_A|, V = sup|v|. The bound is joint: it does not attach to u_A alone.

The pair in (S) is essential, not a proof artifact. For any single witness u with b := (W r(u))|_{B_k} and bᵀu ≠ 0 (generic), the symmetric band-supported correction D = b bᵀ/(bᵀu), δW = 0, gives Du = b and hence r′(u) ≡ 0 on B_k (PROOF.md §7, remark). That D is not conservative (1ᵀD = (1ᵀb) bᵀ/(bᵀu), and 1ᵀb = h²E_k(u) ≠ 0), so (C) excludes it; (S) does not.

**Supplementary statement (symmetric class).** The weighted limiting seam operator ω·L₀, with ω = (2(2+e²))^{−1/2} and L₀ = A∂ₜ + (B−½)∂ₜ² + D∂ₜ³ acting on f, fails formal self-adjointness at every seam point except the equator e = 0. For L = a∂ + b∂² + d∂³ the adjoint defect L − L* has third-order part 2d and first-order part 2a − 2b′ + 3d″; here ωD = −e/(3(e²+2)^{3/2}) ≠ 0 and ωA − (ω(B−½))′ + (3/2)(ωD)″ = e(29e²+103)/(12(e²+2)^{5/2}) ≠ 0 for e ≠ 0. (An earlier statement of the first-order defect omitted the (3/2)d″ term; the conclusion is unchanged.)

## 6. What is excluded (see LIMITATIONS.md)

Nonconservative non-symmetric D; non-diagonal or unbounded relative mass changes; supports growing with h; corrections outside the band; targets imposed on seam rows only; formulations that change the global dof structure; any statement about solution (as opposed to residual) accuracy; the vertex rows.
