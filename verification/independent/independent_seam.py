"""Independent re-derivation of the PX|PY seam-row residual series for the frozen
cubed-sphere operator (design sections B, C.1, C.2, F).

Differences from the package's d1/d2:
  * The seam row (Ku)_b is obtained by symbolic differentiation of the assembled
    quadratic ENERGY (face energies C.1 + cell energies C.2) with respect to u_b,
    not from hand-derived per-term gradient formulas.
  * PY contributions use PY's own chart metric evaluated at PY coordinates
    (no appeal to evenness/oddness in xi).
  * The Fermi form is checked two ways: generic-jet chain rule AND explicit test
    functions u(x,y,z) (polynomials) evaluated through both routes.
Exact sympy arithmetic throughout.
"""
import sympy as sp
from sympy import Rational as R

h, e = sp.symbols("h e", real=True)
X, Y = sp.symbols("X Y", real=True)
ORD = 4
J = {(p, q): sp.Symbol(f"U{p}{q}") for p in range(ORD + 1) for q in range(ORD + 1 - p)}
Upoly = sum(J[p, q] * X**p * Y**q / (sp.factorial(p) * sp.factorial(q)) for (p, q) in J)


def U_at_PX(xi, eta, order):
    """U = u o PX-chart, generic jet at (1,e), evaluated at PX coords (xi,eta), series in h."""
    dx = sp.series(xi - 1, h, 0, order + 1).removeO()
    dy = sp.series(eta - e, h, 0, order + 1).removeO()
    ex = sp.expand(Upoly.subs({X: dx, Y: dy}, simultaneous=True))
    return sum(ex.coeff(h, j) * h**j for j in range(order + 1))


def PY_to_PX(xip, eta):
    """PY chart point (xi',eta) -> sphere (-xi',1,eta)/rho -> PX coords (y/x, z/x)."""
    x, y, z = -xip, 1, eta
    return (y / x, z / x)


def metric(xi, eta):
    d = 1 + xi**2 + eta**2
    return dict(sqrtg=d ** R(-3, 2), g11=(1 + xi**2) * d, g22=(1 + eta**2) * d,
                g12=xi * eta * d)


def face_energy_x(nodeA, nodeB, uA, uB, wy):
    """x-face between (a,b) and (a+1,b): 0.5*wy(b)*(deta/dxi=1)*mean(sqrtg)*mean(g11)*(du)^2"""
    mA, mB = metric(*nodeA), metric(*nodeB)
    return R(1, 2) * wy * R(1, 2) * (mA["sqrtg"] + mB["sqrtg"]) * R(1, 2) * (mA["g11"] + mB["g11"]) * (uB - uA) ** 2


def face_energy_y(nodeA, nodeB, uA, uB, wx):
    mA, mB = metric(*nodeA), metric(*nodeB)
    return R(1, 2) * wx * R(1, 2) * (mA["sqrtg"] + mB["sqrtg"]) * R(1, 2) * (mA["g22"] + mB["g22"]) * (uB - uA) ** 2


def cell_energy(centre, u00, u10, u01, u11):
    """C.2: kappa_c * Pc * Qc * h^2, Pc = 0.5*((u10-u00)+(u11-u01))/h, Qc = 0.5*((u01-u00)+(u11-u10))/h."""
    m = metric(*centre)
    kappa = m["sqrtg"] * m["g12"]
    Pc = R(1, 2) * ((u10 - u00) + (u11 - u01)) / h
    Qc = R(1, 2) * ((u01 - u00) + (u11 - u10)) / h
    return kappa * Pc * Qc * h**2


def seam_row(order=3):
    # node symbols: PX rows a = n-1 (seam, 'S') and n-2 ('P'), PY rows 0 (seam) and 1 ('Q'); eta offsets j
    uS = {j: sp.Symbol(f"uS{j+1}") for j in (-1, 0, 1)}   # shared seam dofs
    uP = {j: sp.Symbol(f"uP{j+1}") for j in (-1, 0, 1)}   # PX xi = 1-h
    uQ = {j: sp.Symbol(f"uQ{j+1}") for j in (-1, 0, 1)}   # PY xi' = -1+h
    E = 0
    # --- PX chart (xi = 1 seam side) ---
    for j in (-1, 0, 1):
        E += face_energy_x((1 - h, e + j*h), (1, e + j*h), uP[j], uS[j], wy=1)   # interior b => wy = 1
    for j in (-1, 0):
        E += face_energy_y((1, e + j*h), (1, e + (j+1)*h), uS[j], uS[j+1], wx=R(1, 2))  # chart boundary a=n-1
        E += face_energy_y((1 - h, e + j*h), (1 - h, e + (j+1)*h), uP[j], uP[j+1], wx=1)
        # cell [n-2,n-1] x [b+j, b+j+1]: corners 00=(n-2,b+j) 10=(n-1,b+j) 01=(n-2,b+j+1) 11=(n-1,b+j+1)
        E += cell_energy((1 - h/2, e + (j + R(1, 2))*h), uP[j], uS[j], uP[j+1], uS[j+1])
    # --- PY chart (xi' = -1 seam side), its OWN metric at its OWN coordinates ---
    for j in (-1, 0, 1):
        E += face_energy_x((-1, e + j*h), (-1 + h, e + j*h), uS[j], uQ[j], wy=1)
    for j in (-1, 0):
        E += face_energy_y((-1, e + j*h), (-1, e + (j+1)*h), uS[j], uS[j+1], wx=R(1, 2))
        E += face_energy_y((-1 + h, e + j*h), (-1 + h, e + (j+1)*h), uQ[j], uQ[j+1], wx=1)
        # cell [0,1] x [b+j,b+j+1]: corners 00=(0,b+j) 10=(1,b+j) 01=(0,b+j+1) 11=(1,b+j+1)
        E += cell_energy((-1 + h/2, e + (j + R(1, 2))*h), uS[j], uQ[j], uS[j+1], uQ[j+1])
    Ku_b = sp.diff(E, uS[0])          # (K u)_b = dE/du_b since E = 1/2 u^T K u
    # substitute the smooth field
    sub = {}
    for j in (-1, 0, 1):
        sub[uS[j]] = U_at_PX(sp.Integer(1), e + j*h, order)
        sub[uP[j]] = U_at_PX(1 - h, e + j*h, order)
        xi_px, eta_px = PY_to_PX(-1 + h, e + j*h)
        sub[uQ[j]] = U_at_PX(xi_px, eta_px, order)
    Ku_b = Ku_b.subs(sub)
    Ku_b = sp.series(sp.expand(Ku_b), h, 0, order + 1).removeO()
    # mass: W_b = sum over aliases w*sqrtg*h^2, w = 1/2 at each chart's boundary index
    W_b = (R(1, 2) * metric(1, e)["sqrtg"] + R(1, 2) * metric(-1, e)["sqrtg"]) * h**2
    return sp.expand(Ku_b), sp.simplify(W_b)


def laplace_beltrami_PX():
    xi, eta = sp.symbols("xi eta", real=True)
    Ug = Upoly.subs({X: xi - 1, Y: eta - e}, simultaneous=True)
    m = metric(xi, eta)
    sg = m["sqrtg"]
    F1 = sg * (m["g11"] * sp.diff(Ug, xi) + m["g12"] * sp.diff(Ug, eta))
    F2 = sg * (m["g12"] * sp.diff(Ug, xi) + m["g22"] * sp.diff(Ug, eta))
    L = (sp.diff(F1, xi) + sp.diff(F2, eta)) / sg
    return sp.expand(sp.simplify(L.subs({xi: 1, eta: e})))


if __name__ == "__main__":
    Ku, W = seam_row(3)
    lap = laplace_beltrami_PX()
    sg1 = metric(1, e)["sqrtg"]
    print("W_b / h^2 =", sp.simplify(W / h**2), "  (expect sqrt(g)(1,e) = (2+e^2)^(-3/2))")
    k = {j: sp.simplify(Ku.coeff(h, j)) for j in range(4)}
    print("Ku h^0 =", k[0]); print("Ku h^1 =", k[1])
    r0 = sp.simplify(-k[2] / sg1 - lap)
    print("r h^0 (consistency, must be 0) =", r0)
    c1 = sp.expand(sp.simplify(-k[3] / sg1))
    print("\nc1 in PX jet basis:")
    for key in sorted(J):
        c = sp.simplify(c1.coeff(J[key]))
        if c != 0:
            print(f"  U_{'xi'*key[0]}{'eta'*key[1]:<12}: {sp.factor(c)}")
    # compare with the package's published closed form
    pub = ((e**2+2)*J[1,0] + (e**2+2)*J[2,0] + e*(e**2+2)*J[1,1]
           + e*(4*e**4+21*e**2+41)/(4*(e**2+2))*J[0,1] - (e**2/2+1)*J[0,2] - e*(e**2+2)**2/6*J[0,3])
    print("\nc1 - published c1 =", sp.simplify(c1 - pub))

    # ---------------- Fermi form, route A: generic jet chain rule ----------------
    xi, eta = sp.symbols("xi eta", real=True)
    N, T = sp.symbols("N T", real=True)
    Fc = {(a, b): sp.Symbol(f"F{a}{b}") for a in range(4) for b in range(4 - a)}
    Fpoly = sum(Fc[a, b] * N**a * T**b / (sp.factorial(a) * sp.factorial(b)) for (a, b) in Fc)
    rho = sp.sqrt(1 + xi**2 + eta**2)
    n_of = sp.asin((1 - xi) / (sp.sqrt(2) * rho))          # p.nu, nu=(1,-1,0)/sqrt2, p=(1,xi,eta)/rho
    t_of = sp.atan(sp.sqrt(2) * eta / (1 + xi))            # atan2(z, (x+y)/sqrt2)
    t0 = sp.atan(e / sp.sqrt(2))
    Fcomp = Fpoly.subs({N: n_of, T: t_of - t0}, simultaneous=True)
    def Ujet(p, q):
        d = Fcomp if (p == q == 0) else sp.diff(Fcomp, xi, p, eta, q)
        return sp.simplify(d.subs({xi: 1, eta: e}))
    subF = {J[key]: Ujet(*key) for key in J if c1.coeff(J[key]) != 0}
    c1F = sp.expand(c1.subs(subF))
    print("\nc1 in Fermi basis (route A):")
    for key in sorted(Fc):
        c = sp.simplify(c1F.coeff(Fc[key]))
        if c != 0:
            print(f"  F_{'n'*key[0]}{'t'*key[1]:<6}: {sp.factor(c)}")
    A_pub = sp.sqrt(2)*e*(47*e**2+139)/(12*(e**2+2)**2)
    B_pub = (3*e**2-2)/(2*(e**2+2))
    D_pub = -sp.sqrt(2)*e/(3*(e**2+2))
    print("A - A_pub =", sp.simplify(c1F.coeff(Fc[0,1]) - A_pub))
    print("B - B_pub =", sp.simplify(c1F.coeff(Fc[0,2]) - B_pub))
    print("C_nn - 1/2 =", sp.simplify(c1F.coeff(Fc[2,0]) - R(1,2)))
    print("D - D_pub =", sp.simplify(c1F.coeff(Fc[0,3]) - D_pub))
    # Laplacian check in Fermi basis
    lapF = sp.simplify(sp.expand(lap.subs({J[key]: Ujet(*key) for key in J if key[0]+key[1] <= 2})))
    print("Delta at seam in Fermi basis (expect F20+F02):", lapF)

    # ---------------- Fermi form, route B: explicit test functions ----------------
    # u(x,y,z) polynomial; compute c1[u] via PX jets of u o PX, and via Fermi derivatives of u o p(n,t).
    x, y, z = sp.symbols("x y z", real=True)
    nn, tt = sp.symbols("n t", real=True)
    tests = [x**2*y - z**3 + x*y*z, x**4 - 3*y**2*z**2 + y*z**3 + x*z, (x+2*y-z)**3 + x*y**2*z]
    ev = {e: R(3, 7)}   # a generic rational seam point
    for u in tests:
        Ux = u.subs({x: 1/rho, y: xi/rho, z: eta/rho})
        jets = {J[key]: sp.diff(Ux, xi, key[0], eta, key[1]).subs({xi: 1, eta: e}) if key != (0,0) else Ux.subs({xi:1, eta:e}) for key in J}
        lhs = sp.simplify(c1.subs(jets).subs(ev))
        # Fermi: p(n,t) = cos n cos t (1,1,0)/sqrt2 + cos n sin t zhat + sin n (1,-1,0)/sqrt2
        px = sp.cos(nn)*sp.cos(tt)/sp.sqrt(2) + sp.sin(nn)/sp.sqrt(2)
        py = sp.cos(nn)*sp.cos(tt)/sp.sqrt(2) - sp.sin(nn)/sp.sqrt(2)
        pz = sp.cos(nn)*sp.sin(tt)
        F = u.subs({x: px, y: py, z: pz})
        tval = sp.atan(e/sp.sqrt(2)).subs(ev)
        at = {nn: 0, tt: tval}
        Ft = sp.diff(F, tt).subs(at); Ftt = sp.diff(F, tt, 2).subs(at); Fnn = sp.diff(F, nn, 2).subs(at); Fttt = sp.diff(F, tt, 3).subs(at)
        rhs = (A_pub*Ft + B_pub*Ftt + R(1,2)*Fnn + D_pub*Fttt).subs(ev)
        print("route B test: c1[u] (PX jets) - (A F_t + B F_tt + F_nn/2 + D F_ttt) =", sp.nsimplify(sp.simplify(lhs - rhs)), "  [value", sp.N(lhs, 15), "]")

    # ---------------- G via integration by parts ----------------
    omega = 1 / sp.sqrt(2 * (2 + e**2))
    Dt = lambda f: sp.simplify(sp.diff(f, e) * (2 + e**2) / sp.sqrt(2))   # d/dt = (de/dt) d/de
    a_, b_, d_ = A_pub, B_pub - R(1, 2), D_pub
    G = sp.simplify(-Dt(omega*a_) + Dt(Dt(omega*b_)) - Dt(Dt(Dt(omega*d_))))
    G_pub = sp.sqrt(2)*(3*e**4 + 31*e**2 - 25)/(4*(e**2+2)**R(5,2))
    print("\nG =", G); print("G - G_pub =", sp.simplify(G - G_pub), "; G(0) =", G.subs(e, 0))
    M1 = sp.simplify(omega*a_ - Dt(omega*b_))
    print("omega*A - (omega*(B-1/2))' =", sp.factor(M1), " [pub: e(6e^4-e^2+19)/(12(e^2+2)^(5/2))]")
    print("omega*D =", sp.simplify(omega*d_))
