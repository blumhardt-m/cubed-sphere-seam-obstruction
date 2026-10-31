"""D3 (step 2): exact series of the frozen-design residual at a generic chart-interior node
(xi0, eta0) -- the stencil used by band rows 1..k (their whole stencil lies inside one chart,
with full trapezoid weights). Claim: r = O(h^2) for every smooth U (h^0, h^1 coefficients 0)."""
import sympy as sp
h, x0, y0 = sp.symbols("h x0 y0", real=True)
X, Y = sp.symbols("X Y", real=True)
ORD = 4
Uc = {(p, q): sp.Symbol(f"U{p}{q}") for p in range(ORD + 1) for q in range(ORD + 1 - p)}
Upoly = sum(Uc[p, q] * X ** p * Y ** q / (sp.factorial(p) * sp.factorial(q)) for (p, q) in Uc)
def U(i, j): return sp.expand(Upoly.subs({X: i * h, Y: j * h}, simultaneous=True))
def d(x, y): return 1 + x ** 2 + y ** 2
def s(x, y): return d(x, y) ** sp.Rational(-3, 2)
def G11(x, y): return (1 + x ** 2) * d(x, y)
def G22(x, y): return (1 + y ** 2) * d(x, y)
def kap(x, y): return x * y * d(x, y) ** sp.Rational(-1, 2)
def ser(e, o=3): return sp.series(e, h, 0, o + 1).removeO()
P = lambda i, j: (x0 + i * h, y0 + j * h)
ub = U(0, 0); Ku = 0
for i in (-1, 1):   # x-faces
    c = sp.Rational(1, 4) * (s(*P(0, 0)) + s(*P(i, 0))) * (G11(*P(0, 0)) + G11(*P(i, 0)))
    Ku += ser(c) * (ub - U(i, 0))
for j in (-1, 1):   # y-faces
    c = sp.Rational(1, 4) * (s(*P(0, 0)) + s(*P(0, j))) * (G22(*P(0, 0)) + G22(*P(0, j)))
    Ku += ser(c) * (ub - U(0, j))
kc = lambda i, j: kap(x0 + i * h / 2, y0 + j * h / 2)        # cell centre offsets (+-1/2)
Ku += ser(kc(-1, -1) / 2) * (ub - U(-1, -1))   # node = corner 4 of cell (a-1,b-1)
Ku += ser(-kc(1, -1) / 2) * (ub - U(1, -1))    # node = corner 3 of cell (a,b-1): -k/2 (u3-u2), u2=(a+1,b-1)
Ku += ser(-kc(-1, 1) / 2) * (ub - U(-1, 1))    # node = corner 2 of cell (a-1,b):  -k/2 (u2-u3), u3=(a-1,b+1)
Ku += ser(kc(1, 1) / 2) * (ub - U(1, 1))       # node = corner 1 of cell (a,b):    k/2 (u1-u4)
Ku = sp.expand(Ku)
xi, eta = sp.symbols("xi eta", real=True)
Ug = Upoly.subs({X: xi - x0, Y: eta - y0}, simultaneous=True)
sg = s(xi, eta); g12 = xi * eta * d(xi, eta)
F1 = sg * (G11(xi, eta) * sp.diff(Ug, xi) + g12 * sp.diff(Ug, eta))
F2 = sg * (g12 * sp.diff(Ug, xi) + G22(xi, eta) * sp.diff(Ug, eta))
lap = sp.expand(((sp.diff(F1, xi) + sp.diff(F2, eta)) / sg).subs({xi: x0, eta: y0}))
s0 = s(x0, y0)
print("Ku h^0:", sp.simplify(Ku.coeff(h, 0)), " Ku h^1:", sp.simplify(Ku.coeff(h, 1)))
print("r h^0 = -Ku_h2/s - lap:", sp.simplify(-Ku.coeff(h, 2) / s0 - lap))
print("r h^1 = -Ku_h3/s:", sp.simplify(-Ku.coeff(h, 3) / s0))
