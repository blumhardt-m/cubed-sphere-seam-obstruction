"""D5 (step 2): rigorous elementary bounds for the witnesses (interval arithmetic, mpmath.iv,
outward-rounded; plus exact rational checks). No external project code.
(i)  sign of G on supp psi_A: G<0 iff 3e^4+31e^2-25<0; on |t|<=0.42, e^2=2tan^2 t <= 2*(0.42/(1-0.42^2/2))^2.
(ii) support of u (|t|<=0.42,|n|<=0.20) and v (|t|<=0.45,|n|<=0.20) lies in {x>0,y>0,|z|<min(x,y)}
     at positive margin: min(x,y) >= (cos n_max cos t_max - sin n_max)/sqrt2 > sin t_max >= |z|.
(iii) omega*D and omega*a-(omega*b)' are nonzero for e != 0 (6e^4-e^2+19 > 0)."""
from fractions import Fraction as Fr
from mpmath import iv, mpf
iv.dps = 30
# (i) exact rational bound: tan t <= t/(1-t^2/2) for 0<=t<1 (sin t<=t, cos t>=1-t^2/2)
t = Fr(42, 100); tb = t / (1 - t * t / 2); e2max = 2 * tb * tb
poly = 3 * e2max ** 2 + 31 * e2max - 25
print("(i) e^2 bound on supp psi_A:", float(e2max), " 3e^4+31e^2-25 at bound:", float(poly), " negative:", poly < 0)
# monotone in e^2 >= 0, so negative on whole support
# (ii) interval check
for name, tm, nm in (("u", "0.42", "0.20"), ("v", "0.45", "0.20")):
    T, Nn = iv.mpf(tm), iv.mpf(nm)
    lhs = (iv.cos(Nn) * iv.cos(T) - iv.sin(Nn)) / iv.sqrt(2)
    rhs = iv.sin(T)
    print(f"(ii) {name}: min(x,y) lower bound {lhs}  > |z| upper bound {rhs}:", lhs.a > rhs.b)
# monotonicity: cos n cos t - sin n decreasing in |n|,|t| on [0,pi/2]; sin t increasing -> extreme corner suffices
# (iii)
print("(iii) 6e^4-e^2+19 = 6(e^2-1/12)^2 + 19 - 1/24 > 0 for all real e:", Fr(19) - Fr(1, 24) > 0)
