"""Independent derivation of G from the chart-INTERIOR truncation density (route 2).

Part A (this block, from indep_interior.py): Lemma 2 and the exact h^2 density r2 of the interior stencil.
Part B (below, from indep_boundary.py): Euler-Lagrange check, seam boundary reduction, Fermi conversion, G.
Self-contained: no pickle, no imports from other project files. Exact sympy arithmetic throughout.

Original docstring of Part A follows.

Independent check of Lemma 2 (interior rows O(h^2)) and of the limit E_inf (Lemma 4)
by a route that never uses the seam row:

   sum_all W r(u) = -1^T K u - Q_W(Delta u) = -Q_W(Delta u)        (K 1 = 0 exactly)
   =>  E_k(u) = h^-2 sum_{B_k} W r  =  -h^-2 Q_W(Delta u) - h^-2 sum_{B_k^c} W r
   For u in U_0: Q_W(Delta u) = O(h^4) (Euler-Maclaurin, d_n(sqrt g Delta u)=0 on the seam),
   and the complement is chart-interior:  h^-2 sum W r -> int sqrt(g) r2[U] dxi deta.
   Hence   E_inf(u) = - int_{PX} sqrt(g) r2[U] - int_{PY} sqrt(g) r2[U],
   with r2 the h^2 truncation density of the chart-INTERIOR stencil.

Part 1: derive r2 exactly (sympy) from the C.1/C.2 energies at a generic interior node.
Part 2: evaluate -int sqrt(g) r2[u_A] over PX and PY numerically (Gauss-Legendre) for the
        package's witness u_A, and compare with int psi_A G dt from the published G.
"""
import sys, time
import sympy as sp
from sympy import Rational as R

h = sp.symbols("h", positive=True)
xi0, eta0 = sp.symbols("xi0 eta0", real=True)
X, Y = sp.symbols("X Y", real=True)
ORD = 4
J = {(p, q): sp.Symbol(f"U{p}{q}") for p in range(ORD + 1) for q in range(ORD + 1 - p)}
Upoly = sum(J[p, q] * X**p * Y**q / (sp.factorial(p) * sp.factorial(q)) for (p, q) in J)


def metric(xi, eta):
    d = 1 + xi**2 + eta**2
    return dict(sqrtg=d ** R(-3, 2), g11=(1 + xi**2) * d, g22=(1 + eta**2) * d, g12=xi * eta * d)


def Uat(i, j):
    return Upoly.subs({X: i * h, Y: j * h}, simultaneous=True)


def fx(nA, nB, uA, uB):
    mA, mB = metric(*nA), metric(*nB)
    return R(1, 2) * R(1, 4) * (mA["sqrtg"] + mB["sqrtg"]) * (mA["g11"] + mB["g11"]) * (uB - uA) ** 2


def fy(nA, nB, uA, uB):
    mA, mB = metric(*nA), metric(*nB)
    return R(1, 2) * R(1, 4) * (mA["sqrtg"] + mB["sqrtg"]) * (mA["g22"] + mB["g22"]) * (uB - uA) ** 2


def cell(centre, u00, u10, u01, u11):
    m = metric(*centre)
    kappa = m["sqrtg"] * m["g12"]
    Pc = R(1, 2) * ((u10 - u00) + (u11 - u01)) / h
    Qc = R(1, 2) * ((u01 - u00) + (u11 - u10)) / h
    return kappa * Pc * Qc * h**2


def interior_r2():
    u = {(i, j): sp.Symbol(f"u_{i+1}_{j+1}") for i in (-1, 0, 1) for j in (-1, 0, 1)}
    node = lambda i, j: (xi0 + i * h, eta0 + j * h)
    E = 0
    for j in (-1, 0, 1):
        for i in (-1, 0):
            E += fx(node(i, j), node(i + 1, j), u[i, j], u[i + 1, j])
    for i in (-1, 0, 1):
        for j in (-1, 0):
            E += fy(node(i, j), node(i, j + 1), u[i, j], u[i, j + 1])
    for i in (-1, 0):
        for j in (-1, 0):
            E += cell((xi0 + (i + R(1, 2)) * h, eta0 + (j + R(1, 2)) * h), u[i, j], u[i + 1, j], u[i, j + 1], u[i + 1, j + 1])
    Ku = sp.diff(E, u[0, 0]).subs({u[i, j]: Uat(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)})
    Ku = sp.series(sp.expand(Ku), h, 0, 5).removeO()
    W = metric(xi0, eta0)["sqrtg"] * h**2
    # exact Laplace-Beltrami at the node
    xi, eta = sp.symbols("xi eta", real=True)
    Ug = Upoly.subs({X: xi - xi0, Y: eta - eta0}, simultaneous=True)
    m = metric(xi, eta)
    F1 = m["sqrtg"] * (m["g11"] * sp.diff(Ug, xi) + m["g12"] * sp.diff(Ug, eta))
    F2 = m["sqrtg"] * (m["g12"] * sp.diff(Ug, xi) + m["g22"] * sp.diff(Ug, eta))
    lap = sp.expand(sp.simplify(((sp.diff(F1, xi) + sp.diff(F2, eta)) / m["sqrtg"]).subs({xi: xi0, eta: eta0})))
    k = {j: sp.simplify(Ku.coeff(h, j)) for j in range(5)}
    sg = metric(xi0, eta0)["sqrtg"]
    r0 = sp.simplify(-k[2] / sg - lap)
    r1 = sp.simplify(-k[3] / sg)
    r2 = sp.expand(sp.simplify(-k[4] / sg))
    return k[0], k[1], r0, r1, r2, lap



# ===================== Part B: boundary reduction and G =====================

xi, eta = sp.symbols("xi eta", real=True)
t_start = time.time()
k0, k1, r0, r1, r2, _lap = interior_r2()
print("Part A (Lemma 2): interior Ku h^0 =", k0, " Ku h^1 =", k1, " r h^0 =", r0, " r h^1 =", r1)
coef = {key: sp.simplify(r2.coeff(J[key])) for key in J}
print("Part A: interior h^2 density r2[U] (coefficient of each jet d_xi^p d_eta^q U at (xi0,eta0)):")
for key in sorted(coef):
    if coef[key] != 0:
        print(f"   ({key[0]},{key[1]}): {sp.factor(coef[key])}")
print(f"[Part A done in {time.time()-t_start:.0f}s]")
coef = {k: sp.simplify(v.subs({xi0: xi, eta0: eta})) for k, v in coef.items() if v != 0}
sg = metric(xi, eta)["sqrtg"]
a = {k: sp.simplify(sg * v) for k, v in coef.items()}          # L[U] = sum a_pq d^pq U

# ---- Step A: Euler-Lagrange expression ----
EL = sum((-1) ** (p + q) * sp.diff(a[(p, q)], xi, p, eta, q) for (p, q) in a)
print("Step A: Euler-Lagrange expression of sqrt(g) r2 =", sp.simplify(EL))

# ---- Step B: seam boundary functional at xi = 1 ----
# B[U] = int_{eta} sum_m b_m(eta) * (d_xi^m U)(1, eta) deta
b = {m: 0 for m in range(4)}
for (p, q), apq in a.items():
    for j in range(p):
        b[p - 1 - j] += (-1) ** (q + j) * sp.diff(apq, xi, j, eta, q).subs(xi, 1)
b = {m: sp.simplify(v) for m, v in b.items()}
print("Step B: boundary coefficients b_m(eta) multiplying d_xi^m U at xi=1:")
for m in range(4):
    print(f"   m={m}: {sp.factor(b[m])}")

# ---- Step C: Fermi conversion on the seam ----
e = sp.symbols("e", real=True)                      # seam parameter eta = e
N, T = sp.symbols("N T", real=True)
ORD = 3
Fc = {(p, q): sp.Symbol(f"F{p}{q}") for p in range(ORD + 1) for q in range(ORD + 1 - p)}
Fpoly = sum(Fc[p, q] * N**p * T**q / (sp.factorial(p) * sp.factorial(q)) for (p, q) in Fc)
rho = sp.sqrt(1 + xi**2 + eta**2)
n_e = sp.asin((1 - xi) / (sp.sqrt(2) * rho))
t_e = sp.atan(sp.sqrt(2) * eta / (1 + xi))
t0 = sp.atan(e / sp.sqrt(2))
Fcomp = Fpoly.subs({N: n_e, T: t_e - t0}, simultaneous=True)
dxiU = {m: sp.simplify((sp.diff(Fcomp, xi, m) if m else Fcomp).subs({xi: 1, eta: e})) for m in range(4)}
integrand_e = sp.expand(sum(b[m].subs(eta, e) * dxiU[m] for m in range(4)))   # d(eta) density against F-jets
# U_0 reductions: F20 = -F02, F21 = -F03, F30 = F10 - F12
red = sp.expand(integrand_e.subs({Fc[2, 0]: -Fc[0, 2], Fc[2, 1]: -Fc[0, 3], Fc[3, 0]: Fc[1, 0] - Fc[1, 2]}))
cof = {k: sp.simplify(red.coeff(Fc[k])) for k in Fc}
print("Step C: seam density (per d eta) after U_0 reduction, coefficients of F-jets:")
for k in sorted(cof):
    if cof[k] != 0:
        print(f"   F_{'n'*k[0]}{'t'*k[1]:<5}: {sp.factor(cof[k])}")
odd = [k for k in cof if cof[k] != 0 and k[0] % 2 == 1]
print("   odd-in-n jets surviving on PX (they change sign on PY by the x<->y mirror and cancel in the sum):", odd)

# change variable eta -> t : d eta = (2+e^2)/sqrt2 dt ; density per dt = cof * (2+e^2)/sqrt2
Dt = lambda g: sp.simplify(sp.diff(g, e) * (2 + e**2) / sp.sqrt(2))
jac = (2 + e**2) / sp.sqrt(2)
c0, c1_, c2_, c3_ = [sp.simplify(cof.get(k, 0) * jac) for k in ((0, 0), (0, 1), (0, 2), (0, 3))]
# int (c0 f + c1 f' + c2 f'' + c3 f''') dt = int f (c0 - c1' + c2'' - c3''') dt
Gtilde = sp.simplify(c0 - Dt(c1_) + Dt(Dt(c2_)) - Dt(Dt(Dt(c3_))))
G_pub = sp.sqrt(2) * (3*e**4 + 31*e**2 - 25) / (4*(e**2 + 2)**R(5, 2))
print("\nPX boundary density against f (per dt):", sp.factor(Gtilde))
print("E_inf density = -2 * that =", sp.factor(-2 * Gtilde))
print("published G                =", G_pub)
print("difference (-2*Gtilde - G_pub) =", sp.simplify(-2 * Gtilde - G_pub))
