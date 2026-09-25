#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figuras de la carta no-go k-essence (2026-09-24).
Genera fig1_exponente.pdf, fig2_curvas.pdf, fig3_exacta.pdf en esta carpeta.

Todo sale de las fórmulas del texto:
  - rho_phi ∝ r^{-p},  p = 4n/(2n-1)            (Teorema 2)
  - v^2 = G M(<r)/r, con M(<r) integrada de verdad (no la ley de potencias local):
      n > 3/2 : v^2 ∝ r^{-2/(2n-1)}
      n = 3/2 : v^2 ∝ ln(r/R0)/r
      n < 3/2 : masa encerrada converge -> v^2 ∝ 1/r (kepleriana)
  - solución exacta F = A(sqrt X - lambda)^2, mu = F' = A(1 - lambda/sqrt X)  (Teorema 1)
Paleta validada (dataviz validate_palette.js, modo light: 5/5 PASS):
  #0072B2  #D55E00  #009E73  #8B4A9C
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
C = ["#0072B2", "#D55E00", "#009E73", "#8B4A9C"]
INK, INK2, GRID = "#1a1a1a", "#555555", "#e4e4e4"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelcolor": INK,
    "ytick.labelcolor": INK, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
})
W = 6.2  # ancho de \linewidth en pulgadas (A4, márgenes 2.6 cm)


def guardar(fig, nombre):
    fig.savefig(os.path.join(AQUI, nombre))
    fig.savefig(os.path.join(AQUI, nombre.replace(".pdf", ".png")), dpi=150)
    plt.close(fig)
    print("  ->", nombre)


# ---------------------------------------------------------------------
# Fig. 1 — exponente asintótico de v^2 contra n
# ---------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(W, 2.6))
n = np.logspace(np.log10(0.62), np.log10(200), 600)
local = -2 / (2 * n - 1)
real = np.where(n > 1.5, local, -1.0)

ax.axhline(0, color=INK2, lw=0.8, ls=(0, (4, 3)))
ax.text(150, 0.04, r"flat: $v^2=$ const", ha="right", va="bottom", color=INK2)
ax.plot(n[n < 1.5], local[n < 1.5], color=C[0], lw=1.2, ls=(0, (2, 2)),
        label=r"local power law $-2/(2n-1)$ (subdominant)")
ax.plot(n, real, color=C[0], lw=2, label=r"actual exponent of $v^2(r)$, $r\to\infty$")

ax.plot([1.5], [-1], "o", ms=8, color=C[1], mec="white", mew=1.5, zorder=5)
ax.annotate(r"$n=3/2$: deep-MOND $\mu\to x$" "\n" r"$v^2\propto \ln r/r$",
            (1.5, -1), (2.4, -1.55), color=INK,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
ax.plot([1], [-1], "o", ms=8, mfc="white", mec=C[3], mew=1.5, zorder=5)
ax.annotate(r"$n=1$: canonical (Fabri)" "\n" r"$S\equiv0$, no source",
            (1, -1), (0.66, -0.5), color=INK,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
ax.annotate(r"$n\to\infty$: $\rho_\phi\propto r^{-2}$, isothermal" "\n"
            "dark-matter branch (Armendariz-Picon & Lim)",
            (100, -2 / 199), (14, -0.62), color=INK,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))

ax.set_xscale("log")
ax.set_xlim(0.6, 200)
ax.set_ylim(-2.1, 0.25)
ax.set_xticks([1, 1.5, 3, 10, 30, 100])
ax.set_xticklabels(["1", "3/2", "3", "10", "30", "100"])
ax.set_xlabel(r"kinetic exponent $n$ in $F=AX^n$")
ax.set_ylabel(r"exponent $e$ in $v^2\propto r^{\,e}$")
ax.legend(loc="lower right")
guardar(fig, "fig1_exponente.pdf")


# ---------------------------------------------------------------------
# Fig. 2 — curvas de rotación del halo escalar (masa encerrada integrada)
# ---------------------------------------------------------------------
def v2(r, nn, R0=0.1):
    """v^2 = M(<r)/r con rho = r^{-p}, integrada desde R0 (normalización arbitraria)."""
    if np.isinf(nn):
        return np.ones_like(r)
    p = 4 * nn / (2 * nn - 1)
    if np.isclose(p, 3):
        M = np.log(r / R0)
    else:
        M = (r ** (3 - p) - R0 ** (3 - p)) / (3 - p)
    return M / r


fig, ax = plt.subplots(figsize=(W, 2.7))
r = np.logspace(0, 2, 400)
casos = [(1.5, r"$n=3/2$ (deep-MOND $\mu$)"), (3, r"$n=3$"),
         (10, r"$n=10$"), (np.inf, r"$n\to\infty$ (isothermal)")]
for (nn, lab), col in zip(casos, C):
    y = np.sqrt(v2(r, nn) / v2(np.array([1.0]), nn)[0])
    ax.plot(r, y, color=col, lw=2, label=lab)
    ax.text(r[-1] * 1.04, y[-1], lab.split(" ")[0], color=INK, va="center")
ax.set_xscale("log")
ax.set_xlim(1, 100)
ax.set_ylim(0, 1.12)
ax.set_xlabel(r"radius $r/r_0$")
ax.set_ylabel(r"$v(r)/v(r_0)$")
ax.legend(loc="lower left", ncol=2)
guardar(fig, "fig2_curvas.pdf")


# ---------------------------------------------------------------------
# Fig. 3 — la solución exacta F = A(sqrt X - lambda)^2 (A = lambda = 1)
# ---------------------------------------------------------------------
lam = 1.0
x = np.linspace(0.02, 3, 500)          # x = sqrt(X)
F = (x - lam) ** 2
mu = 1 - lam / x
mond = x / np.sqrt(1 + x ** 2)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(W, 2.5))
for a in (a1, a2):
    a.axvspan(0, lam, color=C[1], alpha=0.10, lw=0)
    a.set_xlim(0, 3)
    a.set_xlabel(r"$\sqrt{X}=x$")

a1.plot(x, F, color=C[0], lw=2)
a1.plot([0], [lam ** 2], "o", ms=8, color=C[0], mec="white", mew=1.5, clip_on=False, zorder=5)
a1.annotate(r"$F(0)=A\lambda^2\neq0$" "\n(a cosmological constant)",
            (0, 1), (0.9, 2.3), color=INK,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
a1.text(0.5, 3.6, "ghost\n" r"$F'<0$", ha="center", va="top", color=C[1])
a1.set_ylim(0, 4.1)
a1.set_ylabel(r"$F(X)=A(\sqrt{X}-\lambda)^2$")
a1.set_title("(a) kinetic function", loc="left", fontsize=9, color=INK)

a2.axhline(0, color=INK2, lw=0.6)
a2.plot(x, mu, color=C[0], lw=2, label=r"$\mu=F'=A(1-\lambda/x)$")
a2.plot(x, mond, color=INK2, lw=1.4, ls=(0, (4, 3)), label=r"MOND: $x/\sqrt{1+x^2}$")
a2.text(0.5, 0.75, "ghost", ha="center", color=C[1])
a2.annotate(r"$\mu\to-\infty$", (0.18, -3.4), (0.55, -2.6), color=INK,
            arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
a2.set_ylim(-3.5, 1.2)
a2.set_ylabel(r"$\mu(x)$")
a2.set_title(r"(b) interpolation function", loc="left", fontsize=9, color=INK)
a2.legend(loc="lower right")
fig.tight_layout(w_pad=2)
guardar(fig, "fig3_exacta.pdf")
