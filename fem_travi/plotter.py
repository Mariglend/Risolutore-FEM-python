"""
fem_travi/plotter.py
====================
Grafici matplotlib per script e notebook.

    plot_struttura(struttura, risultato=None)   struttura, carichi, deformata, reazioni
    plot_diagrammi(struttura, risultato)        N, V, M disegnati sulla struttura
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

from .assembler import _lunghezza_theta
from .core import CaricoDistribuito, CoppiaNodale, ForzaInCampata, ForzaNodale
from .solver import Risultato, Struttura

C_TRAVE, C_DEF, C_VINC, C_CAR, C_REAZ = "#1f2937", "#d9480f", "#2b8a3e", "#c2410c", "#7c3aed"


def _check_mpl():
    if not HAS_MPL:
        raise ImportError("matplotlib non trovato: pip install matplotlib")


def _span(s: Struttura) -> float:
    xs = [n.x for n in s.nodi]; ys = [n.y for n in s.nodi]
    return max(max(xs) - min(xs), max(ys) - min(ys), 1.0)


def _base(ax, s: Struttura, labels=True, colore=C_TRAVE, lw=2.5):
    for el in s.travi:
        ni, nj = s.nodi[el.nodo_i], s.nodi[el.nodo_j]
        ax.plot([ni.x, nj.x], [ni.y, nj.y], color=colore, lw=1.2 if el.asta else lw, zorder=2)
        if labels:
            ax.text((ni.x + nj.x) / 2, (ni.y + nj.y) / 2, f" T{el.id + 1}", fontsize=7, color="#6b7280")
    for n in s.nodi:
        if n.cerniera or any(el.asta for el in s.travi if n.id in (el.nodo_i, el.nodo_j)):
            ax.plot(n.x, n.y, "o", mfc="white", mec=colore, ms=6, zorder=4)
        else:
            ax.plot(n.x, n.y, "o", color=colore, ms=3, zorder=4)
        if labels:
            ax.text(n.x, n.y, f"  N{n.id + 1}", fontsize=7, color="#2563eb", va="bottom")


def _vincoli(ax, s: Struttura, d: float):
    for v in s.vincoli:
        n = s.nodi[v.nodo]
        (cx, sx), (cy, sy) = v.assi
        tipo = v.tipo or ("incastro" if v.ux_fisso and v.uy_fisso and v.phi_fisso else "")
        if tipo == "incastro" or (v.phi_fisso and v.ux_fisso and v.uy_fisso):
            ax.plot([n.x - d * cx, n.x + d * cx], [n.y - d * sx, n.y + d * sx], color=C_VINC, lw=4)
        elif v.phi_fisso:   # doppio pendolo
            for k in (-0.5, 0.5):
                ox, oy = n.x + k * d * cx, n.y + k * d * sx
                ax.plot([ox, ox - d * cy], [oy, oy - d * sy], color=C_VINC, lw=1.5)
            ax.plot([n.x - d * cx - d * cy, n.x + d * cx - d * cy],
                    [n.y - d * sx - d * sy, n.y + d * sx - d * sy], color=C_VINC, lw=2)
        elif v.ux_fisso or v.uy_fisso:
            p1 = (n.x - 0.6 * d * cx - d * cy, n.y - 0.6 * d * sx - d * sy)
            p2 = (n.x + 0.6 * d * cx - d * cy, n.y + 0.6 * d * sx - d * sy)
            ax.add_patch(plt.Polygon([(n.x, n.y), p1, p2], closed=True, fc="none", ec=C_VINC, lw=1.5))
            base = 1.0 if (v.ux_fisso and v.uy_fisso) else 1.35
            ax.plot([n.x - d * cx - base * d * cy, n.x + d * cx - base * d * cy],
                    [n.y - d * sx - base * d * sy, n.y + d * sx - base * d * sy], color=C_VINC, lw=2)
        else:   # molla
            ax.plot(n.x, n.y - d, marker="$∿$", color=C_VINC, ms=14)


def _carichi(ax, s: Struttura, d: float):
    for c in s.carichi:
        if isinstance(c, ForzaNodale):
            n = s.nodi[c.nodo]; F = math.hypot(c.Fx, c.Fy)
            if F == 0:
                continue
            ux, uy = c.Fx / F, c.Fy / F
            ax.annotate("", xy=(n.x, n.y), xytext=(n.x - 1.6 * d * ux, n.y - 1.6 * d * uy),
                        arrowprops=dict(arrowstyle="-|>", color=C_CAR, lw=1.6))
            ax.text(n.x - 1.7 * d * ux, n.y - 1.7 * d * uy, f"{F:.3g} kN", color=C_CAR, fontsize=7)
        elif isinstance(c, CoppiaNodale):
            n = s.nodi[c.nodo]
            ax.text(n.x, n.y + 0.6 * d, ("↺ " if c.M > 0 else "↻ ") + f"{abs(c.M):.3g} kN·m",
                    color=C_CAR, fontsize=8, ha="center")
        elif isinstance(c, (CaricoDistribuito, ForzaInCampata)):
            el = s.travi[c.trave]
            ni, nj = s.nodi[el.nodo_i], s.nodi[el.nodo_j]
            if isinstance(c, ForzaInCampata):
                L, th = _lunghezza_theta(el, s.nodi)
                px, py = ni.x + c.a * math.cos(th), ni.y + c.a * math.sin(th)
                ax.plot(px, py, "v", color=C_CAR)
                ax.text(px, py + 0.4 * d, f"{math.hypot(c.Fx, c.Fy):.3g} kN", color=C_CAR, fontsize=7)
            else:
                ax.text((ni.x + nj.x) / 2, (ni.y + nj.y) / 2 + 0.8 * d,
                        f"q = {c.qi:.3g}…{c.qj:.3g} kN/m", color=C_CAR, fontsize=7, ha="center")


def plot_struttura(struttura: Struttura, risultato: Optional[Risultato] = None,
                   scala_deformata: Optional[float] = None, mostra_labels: bool = True,
                   mostra_carichi: bool = True, titolo: str = "Struttura FEM",
                   figsize: tuple = (11, 6), salva: Optional[str] = None, mostra: bool = True):
    """Struttura con vincoli e carichi; con ``risultato`` anche deformata e reazioni."""
    _check_mpl()
    s = struttura
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal"); ax.grid(True, alpha=0.25, ls="--"); ax.set_title(titolo)
    d = 0.04 * _span(s)
    _base(ax, s, mostra_labels)
    _vincoli(ax, s, d)
    if mostra_carichi:
        _carichi(ax, s, d)
    if risultato is not None:
        umax = max(np.abs(risultato.U[0::3]).max(), np.abs(risultato.U[1::3]).max(), 1e-300)
        if scala_deformata is None:
            scala_deformata = 0.1 * _span(s) / umax
        for t in range(len(s.travi)):
            df = risultato.deformata(t, scala=scala_deformata)
            ax.plot(df["X"], df["Y"], color=C_DEF, lw=1.5, ls="--",
                    label=f"deformata (×{scala_deformata:.3g})" if t == 0 else None)
        for n, rv in risultato.reazioni().items():
            nd = s.nodi[n]
            for comp, vx, vy in (("Rx", 1, 0), ("Ry", 0, 1)):
                val = rv[comp]
                if abs(val) > 1e-9:
                    sg = np.sign(val)
                    ax.annotate("", xy=(nd.x, nd.y), xytext=(nd.x - 2 * d * vx * sg, nd.y - 2 * d * vy * sg),
                                arrowprops=dict(arrowstyle="-|>", color=C_REAZ, lw=1.8))
                    ax.text(nd.x - 2.2 * d * vx * sg, nd.y - 2.2 * d * vy * sg, f"{abs(val):.3g}",
                            color=C_REAZ, fontsize=7, fontweight="bold")
            if abs(rv["M"]) > 1e-9:
                ax.text(nd.x + d, nd.y - d, ("↺" if rv["M"] > 0 else "↻") + f"{abs(rv['M']):.3g}",
                        color=C_REAZ, fontsize=8, fontweight="bold")
        ax.legend(fontsize=8, loc="best")
    plt.tight_layout()
    if salva:
        plt.savefig(salva, dpi=150, bbox_inches="tight")
    if mostra:
        plt.show()
    return fig


def plot_diagrammi(struttura: Struttura, risultato: Risultato, figsize: tuple = (15, 5),
                   salva: Optional[str] = None, mostra: bool = True):
    """Diagrammi N, V, M disegnati sulla struttura (M dal lato delle fibre tese)."""
    _check_mpl()
    s, r = struttura, risultato
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    nomi = {"N": ("Sforzo normale N [kN]", "#2563eb"),
            "V": ("Taglio V [kN]", "#d9480f"),
            "M": ("Momento flettente M [kN·m]", "#7c3aed")}
    dati = [r.sollecitazioni(t, 81) for t in range(len(s.travi))]
    for ax, k in zip(axes, ("N", "V", "M")):
        titolo, col = nomi[k]
        ax.set_title(titolo, fontsize=10); ax.set_aspect("equal"); ax.axis("off")
        _base(ax, s, labels=False, colore="#9ca3af", lw=1.5)
        vmax = max(max(np.abs(dd[k]).max() for dd in dati), 1e-12)
        h = 0.12 * _span(s) / vmax
        segno = -1.0 if k == "M" else 1.0      # M verso le fibre tese (−y locale)
        for t, el in enumerate(s.travi):
            L, th = _lunghezza_theta(el, s.nodi)
            ni = s.nodi[el.nodo_i]
            c, sn = math.cos(th), math.sin(th)
            x, val = dati[t]["x"], dati[t][k]
            bx, by = ni.x + c * x, ni.y + sn * x
            ox, oy = bx + segno * h * val * (-sn), by + segno * h * val * c
            ax.fill(np.r_[bx, ox[::-1]], np.r_[by, oy[::-1]], color=col, alpha=0.18, lw=0)
            ax.plot(ox, oy, color=col, lw=1.2)
            for idx in {0, len(x) - 1, int(np.argmax(np.abs(val)))}:
                if abs(val[idx]) > 1e-6 * vmax:
                    ax.text(ox[idx], oy[idx], f"{val[idx]:.3g}", fontsize=7, color=col)
    plt.tight_layout()
    if salva:
        plt.savefig(salva, dpi=150, bbox_inches="tight")
    if mostra:
        plt.show()
    return fig
