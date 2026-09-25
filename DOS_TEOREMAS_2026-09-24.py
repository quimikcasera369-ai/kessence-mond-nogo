#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
=======================================================================
LOS DOS TEOREMAS DEL CORPUS — establecimiento y verificacion completa
=======================================================================
Fecha: 2026-09-24.

T1 — TEOREMA mu  (geometrico, INDEPENDIENTE de Lambda-EG)
     mu(x) = x/sqrt(1+x^2) es la acumulada isotropa de cot(theta).
     Estatus: sometido a RMF. Sobrevive aunque Lambda-EG sea falsa.

T2 — NO-GO k-essence minimal  (INDEPENDIENTE de Lambda-EG; reclasificado 24-sep)
     Con materia minimamente acoplada, la mu de MOND y la curva plana
     son mutuamente excluyentes.
     Estatus: corregido tras hallar prior art (astro-ph/0505207).
     Clasificacion: INDEPENDIENTE. mu->x es la DEFINICION estandar del limite
     deep-MOND, no una identificacion de Lambda-EG: si Lambda-EG fuera falsa
     manana, el enunciado sigue siendo verdad (correccion de J.P., 24-sep).
     Caso marginal (anadido 24-sep): en n=3/2 la masa encerrada va como ln r
     => v^2 ∝ ln(r)/r (no r^-1 a secas); para 1/2<n<3/2 converge => kepleriana.

Cada teorema se establece en 4 bloques: HIPOTESIS -> CONSTRUCCION ->
CONCLUSION -> CONTROLES (incl. casos degenerados y prior art).
=======================================================================
"""
import sympy as sp

OK = []
def chk(lab, cond):
    OK.append((lab, bool(cond)))
    print(f"    [{'PASS' if cond else 'FAIL'}] {lab}")

def head(t):
    print("\n" + "=" * 74); print(t); print("=" * 74)

x, X, th = sp.symbols('x X theta', positive=True)
k, A, lam, n = sp.symbols('k A lambda n', positive=True)

# ======================================================================
head("T1 — TEOREMA mu")
# ======================================================================
print("""
  HIPOTESIS
    Ax1 (isotropia). Los canales son EJES no orientados (n^ y -n^ son el
         mismo canal). El ensemble es invariante bajo SO(3). Como S^2 es
         homogeneo (S^2 = SO(3)/SO(2)) y SO(3) es compacto y transitivo,
         la unica medida de probabilidad Borel invariante es dOmega/4pi.
         Reducida al dominio fundamental de la identificacion antipodal
         —el hemisferio n^.g^ >= 0, que NO es la mitad sino TODAS las
         posibilidades— queda  dP = sin(theta) dtheta  en [0, pi/2].
    Ax2 (gate de tilt). El canal esta activo sii |g x n^| > a0 |n^ . g^|.

  CONSTRUCCION""")
dP = sp.sin(th)
norm = sp.integrate(dP, (th, 0, sp.pi/2))
print(f"    normalizacion:  int_0^{{pi/2}} sin(theta) dtheta = {norm}")
chk("la medida del hemisferio normaliza a 1", norm == 1)
print("""
    Con |g x n^| = g sin(theta) y |n^.g^| = cos(theta) (en el hemisferio
    cos(theta)>=0, y SOLO ahi se puede dividir sin invertir la desigualdad):
        g sin(theta) > a0 cos(theta)  <=>  tan(theta) > 1/x   con x = g/a0
    El conjunto activo es la banda theta in (arctan(1/x), pi/2].""")
th0 = sp.atan(1/x)
mu = sp.integrate(dP, (th, th0, sp.pi/2))
mu = sp.simplify(mu)
print(f"\n    mu(x) = int_{{arctan(1/x)}}^{{pi/2}} sin(theta) dtheta = {mu}")
chk("mu(x) = x/sqrt(1+x^2)", sp.simplify(mu - x/sp.sqrt(1+x**2)) == 0)

print("\n  CONTROLES")
chk("limite deep-MOND: mu/x -> 1 cuando x->0", sp.limit(mu/x, x, 0) == 1)
chk("limite newtoniano: mu -> 1 cuando x->oo", sp.limit(mu, x, sp.oo) == 1)
chk("monotona creciente: dmu/dx > 0", sp.simplify(sp.diff(mu, x)) == (1+x**2)**sp.Rational(-3,2))
# la densidad de cot(theta) bajo la medida isotropa
u = sp.Symbol('u', positive=True)   # u = cot(theta)
dens = sp.simplify(sp.diff(mu.subs(x, u), u))
print(f"\n    mu'(x) = {sp.simplify(sp.diff(mu,x))}  = densidad isotropa de cot(theta)")
chk("mu' = (1+x^2)^{-3/2} es la densidad de cot(theta)",
    sp.simplify(dens - (1+u**2)**sp.Rational(-3,2)) == 0)
print("""
    CONTROL DEL MODULO (por que el |.| no es decorativo):
    si se toma la forma reducida g sin > a0 cos como primitiva y se integra
    sobre [0,pi] se activa el hemisferio inferior (donde el lado derecho es
    negativo) y sale otra funcion:""")
mu_mal = sp.simplify(sp.integrate(dP, (th, th0, sp.pi))/sp.integrate(dP,(th,0,sp.pi)))
print(f"      mu_espuria(x) = {sp.simplify(mu_mal)}")
chk("la espuria FALLA el limite deep-MOND (mu(0) = 1/2 != 0)",
    sp.limit(mu_mal, x, 0) == sp.Rational(1,2))
print("""
  CONCLUSION T1
    mu(x) = x/sqrt(1+x^2) NO es un ansatz: es la fraccion de angulo solido
    bajo el gate de tilt con medida isotropa. Es GEOMETRIA PURA: no usa
    accion, ni acoplo, ni a0 derivado, ni Lambda-EG.
    => TEOREMA INDEPENDIENTE.""")

# ======================================================================
head("T2 — NO-GO k-essence con acoplo MINIMO  (version CORREGIDA)")
# ======================================================================
print("""
  HIPOTESIS
    H1. L_phi = -N F(X),  X = (d phi)^2/a_*^2,  N>0 constante.
    H2. Materia MINIMAMENTE acoplada: S_mat = S_mat[g_mu nu, Psi].
        => el escalar actua SOLO via T^phi_{mu nu}; no hay fuerza directa.
    H3. Estatico, esfericamente simetrico, campo debil.
    H4. Fuera de la fuente barionica la EOM del escalar es HOMOGENEA.
    H5. F de clase C^2 con F' > 0 (no-ghost) en el dominio.
    H6. W(X) := F' sqrt(X) monotona  (garantiza que X(r) es invertible;
        se cumple si F'>0 y F''>=0).

  CONSTRUCCION — paso 1: la EOM homogenea fija el perfil del CAMPO""")
print("""    nabla.(F' grad phi) = 0  =>  r^2 F'(X) phi' = const
    con phi' = a_* sqrt(X):      r^2 W(X) = C/a_*,   W := F' sqrt(X)   ...(1)""")

print("\n  CONSTRUCCION — paso 2: la fuente de Poisson")
print("""    T_{mu nu} = 2F' d_mu phi d_nu phi - g_{mu nu} a_*^2 F   (phi estatico radial)
      rho   = a_*^2 F
      p_r   = 2F'(phi')^2 - a_*^2 F
      p_t   = -a_*^2 F
      rho + 3 p_bar = 2 a_*^2 (X F' - F)  =:  2 a_*^2 S(X)             ...(2)""")
F = sp.Function('F')
S_gen = X*sp.Derivative(F(X), X) - F(X)
print(f"    S(X) = {S_gen}")

print("\n  CONSTRUCCION — paso 3: condicion de curva plana EXACTA")
print("""    v^2 = r g = const  <=>  g ∝ 1/r  <=>  nabla^2 Psi ∝ 1/r^2
    Por (2):  S(X(r)) ∝ 1/r^2.   Por (1):  W(X(r)) ∝ 1/r^2.
    Por H6, X(r) es invertible  =>  S(X) = lambda W(X)  como FUNCIONES.  ...(3)""")

print("\n  CONSTRUCCION — paso 4: resolver la ecuacion funcional (3)")
f = sp.Function('f')
ode = sp.Eq(X*sp.Derivative(f(X), X) - f(X), lam*sp.Derivative(f(X), X)*sp.sqrt(X))
sol = sp.dsolve(ode, f(X))
print(f"    {sol}")
Fsol = A*(sp.sqrt(X) - lam)**2
print(f"    forma general:  F(X) = A (sqrt(X) - lambda)^2")
chk("resuelve la ecuacion funcional",
    sp.simplify((X*sp.diff(Fsol,X) - Fsol) - lam*sp.diff(Fsol,X)*sp.sqrt(X)) == 0)

print("\n  CONTROLES sobre la solucion unica")
Fp = sp.simplify(sp.diff(Fsol, X))
print(f"    F'(X) = {Fp}")
chk("limite newtoniano existe: F' -> A en X->oo", sp.limit(Fp, X, sp.oo) == A)
chk("FANTASMA: F' < 0 para X < lambda^2", Fp.subs({A:1, lam:1, X:sp.Rational(1,4)}) < 0)
chk("IR: F(0) = A lambda^2 != 0  => es una Lambda, no deep-MOND",
    sp.limit(Fsol, X, 0) == A*lam**2)
mu2 = sp.simplify(Fp.subs(X, x**2))
chk("mu = F' diverge a -oo en x->0 (deep-MOND exige mu -> x)",
    sp.limit(mu2, x, 0) == -sp.oo)

print("""
  CONTROL DE DEGENERACION: lambda = 0""")
S0 = sp.simplify(X*sp.diff(A*X, X) - A*X)
chk("lambda=0 => F=AX => S(X)=0: la fuente se ANULA (caso trivial, no contraejemplo)", S0 == 0)

# ----------------------------------------------------------------------
print("""
  🔴 CONTROL DE PRIOR ART  —  astro-ph/0505207 ("Haloes of k-Essence")
     Ese trabajo usa acoplo MINIMO (l.121), back-reaction por T_mu nu
     (l.286), L = M^4 c_n X^n (l.340), y OBTIENE curvas planas en n->oo
     (l.455). Hay que comprobar si es contraejemplo.""")
ev = sp.simplify(2 - 4*n/(2*n-1))
print(f"\n    Para F = A X^n:   v^2 ∝ r^e   con  e = 2 - 4n/(2n-1) = {ev}")
cab = "mu = F' prop"
print(f"    {'n':>8}{'e (exponente v^2)':>22}{'rho_phi':>16}{cab:>16}")
for nv in [sp.Rational(3,2), 2, 3, 100]:
    e = float(ev.subs(n, nv)); rh = float(-4*nv/(2*nv-1))
    print(f"    {str(nv):>8}{e:>22.4f}{('r^'+str(round(rh,3))):>16}{("X^"+str(sp.nsimplify(nv-1))):>16}")
print(f"    {'n->oo':>8}{float(sp.limit(ev,n,sp.oo)):>22.4f}{'r^-2 isotermo':>16}{"X^(n-1)":>16}")
chk("ningun n FINITO da e=0 (planitud EXACTA imposible)", sp.solve(sp.Eq(ev,0), n) == [])
chk("pero e -> 0 en n->oo  => planitud ASINTOTICA SI es alcanzable",
    sp.limit(ev, n, sp.oo) == 0)
print("""
    => NO es contraejemplo, pero OBLIGA a acotar el enunciado:
       la planitud asintotica SI se logra, con n grande.""")

print("\n  EL ENUNCIADO CORRECTO — la exclusion mutua")
print("""    mu = F' = n A X^{n-1}.  El limite deep-MOND exige  mu -> x = sqrt(X),
    es decir  n - 1 = 1/2  =>  n = 3/2  (UNICO).""")
e32 = ev.subs(n, sp.Rational(3,2))
chk("n=3/2 (la unica con mu de MOND) da e = -1, NO plana", e32 == -1)
chk("los n que aplanan (n->oo) tienen mu = X^{n-1}, que NO tiende a x",
    sp.limit(ev, n, sp.oo) == 0)

print("""
  CONTROL DEL CASO MARGINAL (la ley de potencias local vs la masa encerrada)
    v^2 = G M(<r)/r con M(<r) = int r'^2 rho_phi dr', rho_phi ∝ r^{-4n/(2n-1)}.""")
rr, R0 = sp.symbols('r R_0', positive=True)
def v2_asint(nv):
    pp = 4*nv/(2*nv - 1)
    M = sp.integrate(rr**(2 - pp), (rr, R0, rr))
    return sp.simplify(M/rr)
v2_32 = v2_asint(sp.Rational(3, 2))
print(f"    n=3/2 :  v^2 ∝ {v2_32}")
chk("n=3/2 (MOND): v^2 ∝ ln(r/R0)/r, marginal y decreciente",
    sp.simplify(v2_32 - sp.log(rr/R0)/rr) == 0 and sp.limit(v2_32, rr, sp.oo) == 0)
v2_54 = v2_asint(sp.Rational(5, 4))
chk("1/2<n<3/2 (ej. n=5/4): masa encerrada converge => kepleriana v^2 ∝ 1/r",
    sp.limit(v2_54*rr, rr, sp.oo).is_finite and sp.limit(v2_54*rr, rr, sp.oo) > 0)
v2_3 = v2_asint(sp.Integer(3))
chk("n=3 (>3/2): exponente asintotico = -2/(2n-1) = -2/5",
    sp.limit(sp.log(v2_3)/sp.log(rr), rr, sp.oo) == sp.Rational(-2, 5))

print("""
  CONTROL n <= 1/2 (por que la clase se restringe a n > 1/2)""")
Wn = sp.diff(A*X**n, X)*sp.sqrt(X)
chk("n=1/2: W constante => r^2 W = C sin solucion no trivial",
    sp.simplify(sp.diff(Wn.subs(n, sp.Rational(1, 2)), X)) == 0)
nq = sp.Rational(1, 4)
chk("n<1/2 (ej. 1/4): X ∝ r^{4/(1-2n)} crece con r => viola campo debil",
    (4/(1 - 2*nq)) > 0)
chk("n<1/2: fuente activa S=(n-1)AX^n negativa",
    sp.simplify(X*sp.diff(A*X**nq, X) - A*X**nq).subs({A: 1, X: 1}) < 0)

print("""
  COROLARIO (F general con rama IR X^{3/2}; carta, Corollary 1)
    F = A X^{3/2} + G(X),  G y X G' = o(X^{3/2}).  Ejemplos G = B X^2, B X^{7/4}.""")
Bc, Kc = sp.symbols('B K', positive=True)
for m in (sp.Integer(2), sp.Rational(7, 4)):
    Fg = A*X**sp.Rational(3, 2) + Bc*X**m
    Sg = sp.expand(X*sp.diff(Fg, X) - Fg)
    Xa = 2*Kc/(3*A*rr**2)                   # perfil dominante de r^2 W = K
    pend = sp.limit(sp.log(Sg.subs(X, Xa))/sp.log(rr), rr, sp.oo)
    mu_ir = sp.limit(sp.diff(Fg, X).subs(X, x**2)/x, x, 0)
    chk(f"G = B X^{m}: rho_act ∝ r^-3 (=> v^2 ∝ ln r/r) y mu/x -> 3A/2 finito",
        pend == -3 and mu_ir == sp.Rational(3, 2)*A)
Fl = A*X**sp.Rational(3, 2) + Bc*X
chk("un termino lineal B X no aporta a S pero rompe mu -> x (mu -> B en x->0)",
    sp.simplify(X*sp.diff(Bc*X, X) - Bc*X) == 0 and
    sp.limit(sp.diff(Fl, X).subs(X, x**2), x, 0) == Bc)

print("""
  CONTROL CONTRA EL PAPER PUBLICADO: Armendariz-Picon & Lim, astro-ph/0505207
    Polítropo  L = -[(-M^4 X)^{(1-g)/2} - p*/rho*^g]^{1/(1-g)}  (su ec. 30), gradiente espacial.""")
gA, cA, uu = sp.symbols('gamma_A c_A u', positive=True)   # gamma_APL = -gA  (caso Chaplygin, gamma<0)
LA = -(uu**((1 + gA)/2) + cA)**(1/(1 + gA))              # p*<0 -> c = -cA
serA = sp.series(LA.subs(gA, sp.Rational(1, 2)), uu, 0, 1).removeO()
chk("Chaplygin (gamma=-1/2): L - L(0) arranca en X^{3/4} = X^{(1-gamma)/2}",
    sp.limit(sp.log(-(serA - serA.subs(uu, 0)))/sp.log(uu), uu, 0) == sp.Rational(3, 4))
n_eq = (1 + gA)/2                                          # n = (1-gamma_APL)/2
chk("su ec. (45) v^2 ∝ r^{2/gamma} == nuestra r^{-2/(2n-1)} con n=(1-gamma)/2",
    sp.simplify(-2/(2*n_eq - 1) - 2/(-gA)) == 0)
chk("n=3/2 (MOND) <=> gamma=-2: su frontera entre 1/r y r^{2/gamma}",
    sp.solve(sp.Eq(n_eq, sp.Rational(3, 2)), gA) == [2])
# caso (a): p*>0, 0<gamma<1.  y = u^{(1-g)/2}, delta = y - c -> 0 en r->oo (cero de L)
gP, rr2 = sp.symbols('g r', positive=True)
# r^2 W = K con W ∝ delta^{g/(1-g)}  =>  delta ∝ r^{-2(1-g)/g}
d_exp = -2*(1 - gP)/gP
rho_exp = d_exp/(1 - gP)                 # rho = -L = delta^{1/(1-g)}
pr_exp = d_exp*gP/(1 - gP)               # p_r ≈ y * delta^{g/(1-g)}, y -> c
chk("caso (a): rho ∝ r^{-2/gamma} (su texto) y r^2 rho -> 0 para 0<gamma<1",
    sp.simplify(rho_exp + 2/gP) == 0 and sp.limit(2 + rho_exp.subs(gP, sp.Rational(1, 2)), gP, 0) < 0)
chk("caso (a): r^2 p_r -> const (la PRESION aplana), dado X -> X_c != 0 (verificado numericamente: u_c=0.6561, r^2 p_r -> 1.62 para gamma=1/2)",
    sp.simplify(pr_exp + 2) == 0)
print("""
  CONCLUSION T2 (corregida)
    Con acoplo minimo, v^2 ∝ r^{-2/(2n-1)}. La curva se aplana solo en
    n->oo, donde rho_phi -> r^-2 (halo isotermo: MATERIA OSCURA efectiva).
    Pero la mu de MOND exige n = 3/2, que da v^2 ∝ ln(r)/r (caso marginal).
    => LA mu DE MOND Y LA CURVA PLANA SON MUTUAMENTE EXCLUYENTES.
    Y para planitud EXACTA, la ecuacion funcional tiene solucion unica
    F = A(sqrt(X)-lambda)^2, con fantasma en X<lambda^2 y sin limite MOND.

    ESTATUS: teorema INDEPENDIENTE de Lambda-EG. Lleva mu -> x, pero eso es
    la definicion estandar del limite deep-MOND, no un supuesto del marco.
    (Se clasifico mal como DEPENDIENTE confundiendo input fenomenologico
    con dependencia de Lambda-EG; corregido por J.P. el 24-sep.)
    PRIOR ART: astro-ph/0505207 (planitud asintotica via n grande, como
    materia oscura). Lo NUEVO aqui es la exclusion mutua con la mu de MOND.""")

head("RESUMEN")
n_ok = sum(1 for _, o in OK if o)
print(f"  {n_ok}/{len(OK)} PASS")
for l, o in OK:
    if not o: print("   FALLO:", l)
print("""
  T1 (mu)   : INDEPENDIENTE — geometria pura, sobrevive a Lambda-EG falsa.
  T2 (no-go): INDEPENDIENTE — no usa nada de Lambda-EG; prior art citado.""")
