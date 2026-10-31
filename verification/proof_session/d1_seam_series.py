"""D1 (step 2): exact series of the frozen-design seam-row residual at a PX|PY seam node.

Independent implementation from the design text (C.1 face energies, C.2 cell energies,
B mass, F convention) and the gnomonic chart metric. Imports no external project
code. All algebra exact (sympy); no floating point in the derivation.

U(xi, eta) = u o (PX gnomonic chart), extended smoothly past xi = 1 (x > 0 hemisphere).
Seam node: PX (xi=1, eta=e) == PY (xi'=-1, eta=e).
PY node (xi'=-1+h, eta=e+jh) is the point (1-h, 1, e+jh)/rho, i.e. PX coords
(1/(1-h), (e+jh)/(1-h)).
Generic U: Taylor polynomial at (1, e) with symbolic coefficients U[p,q] = d_xi^p d_eta^q U.
"""
import sympy as sp

h, e = sp.symbols("h e", real=True)
X, Y = sp.symbols("X Y", real=True)          # displacements xi-1, eta-e
ORD = 4
Uc = {(p, q): sp.Symbol(f"U{p}{q}") for p in range(ORD + 1) for q in range(ORD + 1 - p)}
Upoly = sum(Uc[p, q] * X ** p * Y ** q / (sp.factorial(p) * sp.factorial(q)) for (p, q) in Uc)

def U_at(xi, eta, order):
    """U at (xi, eta) expressed as a series in h up to h^order."""
    dx = sp.series(xi - 1, h, 0, order + 1).removeO()
    dy = sp.series(eta - e, h, 0, order + 1).removeO()
    expr = sp.expand(Upoly.subs({X: dx, Y: dy}, simultaneous=True))
    return sum(expr.coeff(h, j) * h ** j for j in range(order + 1))

def delta(xi, eta): return 1 + xi ** 2 + eta ** 2
def s(xi, eta): return delta(xi, eta) ** sp.Rational(-3, 2)            # sqrt(g)
def G11(xi, eta): return (1 + xi ** 2) * delta(xi, eta)                # g^11
def G22(xi, eta): return (1 + eta ** 2) * delta(xi, eta)               # g^22
def kap(xi, eta): return s(xi, eta) * xi * eta * delta(xi, eta)        # sqrt(g) g^12

def ser(expr, order):
    return sp.series(expr, h, 0, order + 1).removeO()

def seam_row(order=3):
    """(K u)_b and W_b at the seam node, as series in h. Coefficients from C.1/C.2:
    x-face (a,a+1) weight wy(b)=1 for an interior b; y-face along the chart edge has wx=1/2
    on each chart (both charts contribute: total full weight); cell energy
    kappa_c*Pc*Qc*h^2 = kappa_c/4*((u4-u1)^2-(u2-u3)^2), corners 1=(a,b) 2=(a+1,b) 3=(a,b+1) 4=(a+1,b+1)."""
    ub = U_at(sp.Integer(1), e, order)
    # PX nodes, row n-2 (xi = 1-h); PY nodes, row 1 (xi' = -1+h) in PX coords
    PX1 = {j: U_at(1 - h, e + j * h, order) for j in (-1, 0, 1)}
    PY1 = {j: U_at(1 / (1 - h), (e + j * h) / (1 - h), order) for j in (-1, 0, 1)}
    S0 = {j: U_at(sp.Integer(1), e + j * h, order) for j in (-1, 1)}
    Ku = 0
    # y-faces along the seam: PX (wx=1/2) + PY (wx=1/2, same metric values since even in xi)
    for j in (-1, 1):
        cPX = sp.Rational(1, 2) * sp.Rational(1, 4) * (s(1, e) + s(1, e + j * h)) * (G22(1, e) + G22(1, e + j * h))
        cPY = sp.Rational(1, 2) * sp.Rational(1, 4) * (s(-1, e) + s(-1, e + j * h)) * (G22(-1, e) + G22(-1, e + j * h))
        Ku += ser(cPX + cPY, order) * (ub - S0[j])
    # PX: x-face between (n-2,b) and (n-1,b)
    cx = sp.Rational(1, 4) * (s(1, e) + s(1 - h, e)) * (G11(1, e) + G11(1 - h, e))
    Ku += ser(cx, order) * (ub - PX1[0])
    # PX cells: node b is corner 2 of upper cell (a=n-2; cols b,b+1) and corner 4 of lower cell
    kup = kap(1 - h / 2, e + h / 2); klo = kap(1 - h / 2, e - h / 2)
    Ku += ser(-kup / 2, order) * (ub - PX1[1])     # dE/du2 = -kappa/2 (u2-u3)
    Ku += ser(klo / 2, order) * (ub - PX1[-1])     # dE/du4 =  kappa/2 (u4-u1)
    # PY: x-face between (0,b) and (1,b)
    cxp = sp.Rational(1, 4) * (s(-1, e) + s(-1 + h, e)) * (G11(-1, e) + G11(-1 + h, e))
    Ku += ser(cxp, order) * (ub - PY1[0])
    # PY cells (a=0): node b is corner 1 of upper cell, corner 3 of lower cell
    kupP = kap(-1 + h / 2, e + h / 2); kloP = kap(-1 + h / 2, e - h / 2)
    Ku += ser(-kupP / 2, order) * (PY1[1] - ub)    # dE/du1 = -kappa/2 (u4-u1)
    Ku += ser(kloP / 2, order) * (PY1[-1] - ub)    # dE/du3 =  kappa/2 (u2-u3)
    W = (sp.Rational(1, 2) * s(1, e) + sp.Rational(1, 2) * s(-1, e)) * h ** 2
    return sp.expand(Ku), W

def lap_exact():
    """Delta U at (1,e) in PX coordinates, from (1/sqrt g) d_i(sqrt g g^ij d_j U)."""
    xi, eta = sp.symbols("xi eta", real=True)
    Ug = Upoly.subs({X: xi - 1, Y: eta - e}, simultaneous=True)
    sg = s(xi, eta); g11 = G11(xi, eta); g22 = G22(xi, eta); g12 = xi * eta * delta(xi, eta)
    F1 = sg * (g11 * sp.diff(Ug, xi) + g12 * sp.diff(Ug, eta))
    F2 = sg * (g12 * sp.diff(Ug, xi) + g22 * sp.diff(Ug, eta))
    L = (sp.diff(F1, xi) + sp.diff(F2, eta)) / sg
    return sp.expand(sp.simplify(L.subs({xi: 1, eta: e})))

if __name__ == "__main__":
    import json, sys, time
    t0 = time.time()
    Ku, W = seam_row(3)
    lap = lap_exact()
    s1 = s(1, e)
    coeff = {j: sp.expand(Ku.coeff(h, j)) for j in range(4)}
    out = {}
    print("Ku h^0:", sp.simplify(coeff[0]))
    print("Ku h^1:", sp.simplify(coeff[1]))
    # residual r = -Ku/W - lap = -(Ku/h^2)/s1 - lap ; Ku/h^2 = coeff2 + h coeff3 + ...
    r0 = sp.simplify(-coeff[2] / s1 - lap)
    r1 = sp.expand(sp.simplify(-coeff[3] / s1))
    print("r h^0 (must be 0):", r0)
    c1 = {k: sp.simplify(r1.coeff(Uc[k])) for k in Uc}
    rest = sp.simplify(r1 - sum(c1[k] * Uc[k] for k in Uc))
    print("r h^1 remainder outside U-basis (must be 0):", rest)
    for k in sorted(c1):
        if c1[k] != 0:
            print(f"  c1[U{k[0]}{k[1]}] =", c1[k])
    out = {"Ku_h0": str(sp.simplify(coeff[0])), "Ku_h1": str(sp.simplify(coeff[1])),
           "r_h0": str(r0), "r_h1_rest": str(rest),
           "c1_PX_basis": {f"U{k[0]}{k[1]}": str(v) for k, v in c1.items() if v != 0}}
    json.dump(out, open(sys.argv[1], "w"), indent=1)
    print(f"elapsed {time.time()-t0:.1f}s")
