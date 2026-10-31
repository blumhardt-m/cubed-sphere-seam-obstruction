"""D2 (step 2): convert the exact seam coefficient c1 (D1) to Fermi coordinates of the seam
great circle, reduce on the witness class U0, and derive the limiting functionals exactly.

Fermi coords (as in the original numerical check): nu=(1,-1,0)/sqrt2, n = arcsin(p.nu),
t = atan2(p.z, p.(1,1,0)/sqrt2).  On the seam e = sqrt2 tan t.
U(xi,eta) = F(n(xi,eta), t(xi,eta)); F_ab = d_n^a d_t^b F at (0, t0).
Witness class U0: smooth, compact support near the open seam segment, with
Delta u = 0 and d_n Delta u = 0 on the seam  =>  F_nn = -F_tt, F_nnt = -F_ttt,
F_nnn = F_n - F_ntt   (from Delta = d_nn - tan n d_n + sec^2 n d_tt).
Pure exact algebra (sympy)."""
import json, sys
import sympy as sp
sys.path.insert(0, sys.argv[2])            # this package's checks/ directory only
from d1_seam_series import seam_row, lap_exact, Uc, s, h, e

xi, eta = sp.symbols("xi eta", real=True)
N, T = sp.symbols("N T", real=True)
ORD = 3
Fc = {(a, b): sp.Symbol(f"F{a}{b}") for a in range(ORD + 1) for b in range(ORD + 1 - a)}
Fpoly = sum(Fc[a, b] * N ** a * T ** b / (sp.factorial(a) * sp.factorial(b)) for (a, b) in Fc)

rho = sp.sqrt(1 + xi ** 2 + eta ** 2)
n_expr = sp.asin((1 - xi) / (sp.sqrt(2) * rho))
t0 = sp.atan(e / sp.sqrt(2))
# t = atan2(eta/rho, (1+xi)/(sqrt2 rho)) = atan(sqrt2*eta/(1+xi)) since 1+xi>0
t_expr = sp.atan(sp.sqrt(2) * eta / (1 + xi))

# U_pq = d_xi^p d_eta^q F(n,t) at (1,e), via multivariate chain rule (exact)
def Upq_in_F(p, q):
    expr = Fpoly.subs({N: n_expr, T: t_expr - t0}, simultaneous=True)
    d = sp.diff(expr, xi, p, eta, q) if (p or q) else expr
    return sp.simplify(d.subs({xi: 1, eta: e}))

Ku, W = seam_row(3)
lap = lap_exact()
s1 = s(1, e)
r1 = sp.expand(-Ku.coeff(h, 3) / s1)
c1U = {k: sp.simplify(r1.coeff(Uc[k])) for k in Uc}
needed = [k for k in c1U if c1U[k] != 0]
sub = {Uc[k]: Upq_in_F(*k) for k in needed}
c1F = sp.expand(sum(c1U[k] * sub[Uc[k]] for k in needed))
coefF = {k: sp.simplify(c1F.coeff(Fc[k])) for k in Fc}
print("c1 in Fermi basis (exact):")
for k in sorted(coefF):
    if coefF[k] != 0:
        print(f"  F{k[0]}{k[1]}: ", coefF[k])
# also express laplacian in F basis as a consistency check: should be F20 + F02 at n=0
lapF = sp.expand(lap.subs({Uc[k]: Upq_in_F(*k) for k in Uc if k[0] + k[1] <= 2}))
print("Delta at seam in F basis (should be F20+F02):", sp.simplify(lapF))

# convert coefficients to functions of t: e = sqrt2 tan t
tt = sp.symbols("t", real=True)
def in_t(expr): return sp.simplify(expr.subs(e, sp.sqrt(2) * sp.tan(tt)))
A_ = coefF[(0, 1)]; B_ = coefF[(0, 2)]; Cnn = coefF[(2, 0)]; D_ = coefF[(0, 3)]; Etnn = coefF[(2, 1)]
odd = {k: coefF[k] for k in coefF if k[0] % 2 == 1}
print("odd-in-n coefficients (mirror-odd data):", {f"F{k[0]}{k[1]}": str(v) for k, v in odd.items()})
# reduction on U0: F20 -> -F02, F21 -> -F03, F30 -> F10 - F12
red = sp.expand(c1F.subs({Fc[2, 0]: -Fc[0, 2], Fc[2, 1]: -Fc[0, 3], Fc[3, 0]: Fc[1, 0] - Fc[1, 2]}))
redc = {k: sp.simplify(red.coeff(Fc[k])) for k in Fc}
print("c1 reduced on U0:")
for k in sorted(redc):
    if redc[k] != 0:
        print(f"  F{k[0]}{k[1]}: ", redc[k])
omega_e = 1 / sp.sqrt(2 * (2 + e ** 2))           # = s(1,e) * d eta/dt
assert sp.simplify(omega_e - s1 * sp.diff(sp.sqrt(2) * sp.tan(tt), tt).subs(tt, t0)) == 0
# L0 = a d + b d^2 + d3 d^3 acting on f(t) = F|_{n=0}
a = redc[(0, 1)]; b = redc[(0, 2)]; d3 = redc[(0, 3)]
wa, wb, wd = [sp.simplify(omega_e * c) for c in (a, b, d3)]
# derivatives in t: d/dt = (de/dt) d/de, de/dt = sqrt2 sec^2 t = (2+e^2)/sqrt2
Dt = lambda f: sp.simplify(sp.diff(f, e) * (2 + e ** 2) / sp.sqrt(2))
G = sp.simplify(-Dt(wa) + Dt(Dt(wb)) - Dt(Dt(Dt(wd))))      # E_inf(f) = int f G dt
M1 = sp.simplify(wa - Dt(wb))                                # symmetric-class first-order condition
print("omega*a =", wa); print("omega*b =", wb); print("omega*d =", wd)
print("G (conservation adjoint density) =", G)
print("omega*a - (omega*b)' =", M1)
vals = {}
for ev in (sp.Integer(0), sp.Rational(1, 2), sp.Integer(1), sp.Rational(1, 4)):
    vals[str(ev)] = {"G": str(sp.nsimplify(sp.simplify(G.subs(e, ev)))), "G_float": float(G.subs(e, ev)),
                     "omega_d": str(sp.simplify(wd.subs(e, ev))), "omega_d_float": float(wd.subs(e, ev)),
                     "M1": str(sp.simplify(M1.subs(e, ev))), "M1_float": float(M1.subs(e, ev))}
    print(ev, vals[str(ev)])
json.dump({"c1_fermi": {f"F{k[0]}{k[1]}": str(v) for k, v in coefF.items() if v != 0},
           "c1_reduced_U0": {f"F{k[0]}{k[1]}": str(v) for k, v in redc.items() if v != 0},
           "omega": str(omega_e), "omega_a": str(wa), "omega_b": str(wb), "omega_d": str(wd),
           "G": str(G), "M1": str(M1), "values": vals,
           "coeffs_as_function_of_e": {"A": str(A_), "B": str(B_), "Cnn": str(Cnn), "D": str(D_), "Etnn": str(Etnn)}},
          open(sys.argv[1], "w"), indent=1)
