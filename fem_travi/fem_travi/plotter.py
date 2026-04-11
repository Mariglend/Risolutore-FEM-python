"""
fem_travi/plotter.py
====================
Visualizzazione della struttura e dei risultati FEM con matplotlib.

Funzioni esposte
----------------
    plot_struttura(struttura, risultato=None, ...)
    plot_diagrammi(struttura, risultato, ...)
"""

from __future__ import annotations

import math
from typing import List, Optional

import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyArrowPatch
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

from .core import Carico, Nodo, Trave, Vincolo
from .solver import Risultato, Struttura
from .assembler import _lunghezza_theta, matrice_rotazione


def _check_mpl():
    if not HAS_MPL:
        raise ImportError(
            "matplotlib non trovato. Installalo con: pip install matplotlib"
        )


# ---------------------------------------------------------------------------
# Plot struttura + deformata
# ---------------------------------------------------------------------------

def plot_struttura(
    struttura: Struttura,
    risultato: Optional[Risultato] = None,
    scala_deformata: Optional[float] = None,
    mostra_labels: bool = True,
    mostra_carichi: bool = True,
    titolo: str = "Struttura FEM",
    figsize: tuple = (12, 6),
    salva: Optional[str] = None,
) -> None:
    """
    Disegna la struttura con nodi, travi, vincoli e (opzionalmente)
    la deformata amplificata.

    Parametri
    ----------
    struttura : Struttura
    risultato : Risultato | None
        Se fornito, disegna anche la deformata.
    scala_deformata : float | None
        Fattore di amplificazione. Se None viene calcolato automaticamente.
    mostra_labels : bool
        Mostra etichette nodi e travi.
    mostra_carichi : bool
        Mostra frecce dei carichi.
    salva : str | None
        Percorso file per salvare la figura (es. 'struttura.png').
    """
    _check_mpl()

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(titolo)

    nodi = struttura.nodi
    travi = struttura.travi

    # -- travi non deformate --------------------------------------------------
    for el in travi:
        ni, nj = nodi[el.nodo_i], nodi[el.nodo_j]
        ax.plot([ni.x, nj.x], [ni.y, nj.y],
                color="#534AB7", linewidth=2.5, zorder=2,
                label="Trave" if el.id == 0 else "")
        if mostra_labels:
            mx, my = (ni.x + nj.x) / 2, (ni.y + nj.y) / 2
            L, theta = _lunghezza_theta(el, nodi)
            perp = (-math.sin(theta), math.cos(theta))
            offset = 0.05 * max(
                max(n.x for n in nodi) - min(n.x for n in nodi), 0.5
            )
            ax.text(mx + perp[0]*offset, my + perp[1]*offset,
                    f"T{el.id+1}\nL={L:.2f}m",
                    fontsize=7, ha="center", color="#3C3489", zorder=5)

    # -- deformata ------------------------------------------------------------
    if risultato is not None:
        U = risultato.U
        abs_U = np.abs(U)
        max_u = float(abs_U[abs_U > 0].max()) if (abs_U > 0).any() else 1.0
        if scala_deformata is None:
            span = max(
                max(n.x for n in nodi) - min(n.x for n in nodi),
                max(n.y for n in nodi) - min(n.y for n in nodi),
                0.1,
            )
            scala_deformata = span * 0.1 / max_u

        for el in travi:
            ni, nj = nodi[el.nodo_i], nodi[el.nodo_j]
            L, theta = _lunghezza_theta(el, nodi)
            c, s = math.cos(theta), math.sin(theta)
            i, j = el.nodo_i, el.nodo_j
            u = [U[i*3]*scala_deformata, U[i*3+1]*scala_deformata,
                 U[i*3+2]*scala_deformata,
                 U[j*3]*scala_deformata, U[j*3+1]*scala_deformata,
                 U[j*3+2]*scala_deformata]

            # rotazione in coord locali
            ul_i =  c*u[0] + s*u[1]
            vl_i = -s*u[0] + c*u[1]
            phi_i = u[2]
            ul_j =  c*u[3] + s*u[4]
            vl_j = -s*u[3] + c*u[4]
            phi_j = u[5]

            steps = 30
            xs, ys = [], []
            for k in range(steps + 1):
                xi = k / steps
                N1 = 1 - 3*xi**2 + 2*xi**3
                N2 = xi * (1-xi)**2 * L
                N3 = 3*xi**2 - 2*xi**3
                N4 = xi**2 * (xi-1) * L
                ul = ul_i*(1-xi) + ul_j*xi
                vl = N1*vl_i + N2*phi_i + N3*vl_j + N4*phi_j
                xg = ni.x + (nj.x - ni.x)*xi + c*ul - s*vl
                yg = ni.y + (nj.y - ni.y)*xi + s*ul + c*vl
                xs.append(xg)
                ys.append(yg)

            ax.plot(xs, ys, color="#D85A30", linewidth=1.5,
                    linestyle="--", alpha=0.8,
                    label=f"Deformata (×{scala_deformata:.1f})" if el.id == 0 else "")

    # -- nodi -----------------------------------------------------------------
    for n in nodi:
        ax.plot(n.x, n.y, "o", color="white", markeredgecolor="#185FA5",
                markeredgewidth=2, markersize=8, zorder=4)
        if mostra_labels:
            ax.text(n.x + 0.05, n.y + 0.05, f"N{n.id+1}",
                    fontsize=8, color="#185FA5", zorder=5, fontweight="bold")

    # -- vincoli --------------------------------------------------------------
    _disegna_vincoli(ax, struttura.vincoli, nodi)

    # -- carichi --------------------------------------------------------------
    if mostra_carichi and struttura.carichi:
        _disegna_carichi(ax, struttura.carichi, struttura.travi, nodi)

    # -- reazioni (se disponibili) --------------------------------------------
    if risultato is not None:
        _disegna_reazioni(ax, risultato, nodi)

    handles, labels_leg = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=8, loc="best")

    plt.tight_layout()
    if salva:
        plt.savefig(salva, dpi=150, bbox_inches="tight")
        print(f"Figura salvata in: {salva}")
    plt.show()


# ---------------------------------------------------------------------------
# Diagrammi interni M, V, N
# ---------------------------------------------------------------------------

def plot_diagrammi(
    struttura: Struttura,
    risultato: Risultato,
    figsize: tuple = (14, 10),
    salva: Optional[str] = None,
) -> None:
    """
    Diagrammi del momento flettente M, taglio V e sforzo normale N
    per ogni trave.
    """
    _check_mpl()

    nodi = struttura.nodi
    travi = struttura.travi
    U = risultato.U
    n_travi = len(travi)

    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=False)
    titles = ["Momento flettente M [kN·m]",
              "Taglio V [kN]",
              "Sforzo normale N [kN]"]
    colors = ["#185FA5", "#D85A30", "#3B6D11"]

    for ax, title, color in zip(axes, titles, colors):
        ax.set_title(title, fontsize=11)
        ax.axhline(0, color="black", linewidth=0.8)
        ax.grid(True, alpha=0.3, linestyle="--")
        ax.set_ylabel(title.split(" ")[-2] + " " + title.split(" ")[-1])

    x_offset = 0.0
    xticks, xticklabels = [0.0], ["0"]

    for el in travi:
        ni_node, nj_node = nodi[el.nodo_i], nodi[el.nodo_j]
        L, theta = _lunghezza_theta(el, nodi)
        i, j = el.nodo_i, el.nodo_j

        # spostamenti locali
        T = matrice_rotazione(theta)
        u_glob = np.array([
            U[i*3], U[i*3+1], U[i*3+2],
            U[j*3], U[j*3+1], U[j*3+2],
        ])
        u_loc = T @ u_glob   # u_loc: [ul_i, vl_i, phi_i, ul_j, vl_j, phi_j]

        EI, EA = el.EI, el.EA
        steps = 50
        xs = np.linspace(0, L, steps + 1)

        M_vals = []
        V_vals = []
        N_vals = []

        for xi in xs:
            s = xi / L
            # Momento flettente da spostamenti locali (Euler-Bernoulli)
            vl_i, phi_i = u_loc[1], u_loc[2]
            vl_j, phi_j = u_loc[4], u_loc[5]
            ul_i, ul_j  = u_loc[0], u_loc[3]

            # derivata seconda della forma deformata → curvatura → M = EI * κ
            # Funzioni di Hermite per la flessione
            d2N1 = (-6 + 12*s) / L**2
            d2N2 = (-4 + 6*s) / L
            d2N3 = ( 6 - 12*s) / L**2
            d2N4 = (-2 + 6*s) / L
            kappa = d2N1*vl_i + d2N2*phi_i + d2N3*vl_j + d2N4*phi_j
            M = EI * kappa

            # Taglio V = -EI * d³v/dx³  (costante per trave senza carichi)
            d3N1 = 12 / L**3
            d3N3 = -12 / L**3
            d3 = d3N1*vl_i + d3N3*vl_j + (6/L**2)*(phi_i - phi_j)
            V = -EI * d3

            # Sforzo normale N = EA * du/dx
            N = EA * (ul_j - ul_i) / L

            M_vals.append(M)
            V_vals.append(V)
            N_vals.append(N)

        xs_plot = xs + x_offset
        axes[0].fill_between(xs_plot, M_vals, 0,
                             alpha=0.25, color=colors[0])
        axes[0].plot(xs_plot, M_vals, color=colors[0], linewidth=1.5)

        axes[1].fill_between(xs_plot, V_vals, 0,
                             alpha=0.25, color=colors[1])
        axes[1].plot(xs_plot, V_vals, color=colors[1], linewidth=1.5)

        axes[2].fill_between(xs_plot, N_vals, 0,
                             alpha=0.25, color=colors[2])
        axes[2].plot(xs_plot, N_vals, color=colors[2], linewidth=1.5)

        x_offset += L
        xticks.append(x_offset)
        xticklabels.append(f"{x_offset:.2f}")

        for ax in axes:
            ax.axvline(x_offset, color="gray", linewidth=0.5, linestyle=":")

    for ax in axes:
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=8)
        ax.set_xlabel("Posizione lungo le travi [m]")

    plt.suptitle("Diagrammi delle azioni interne", fontsize=13, fontweight="bold")
    plt.tight_layout()
    if salva:
        plt.savefig(salva, dpi=150, bbox_inches="tight")
        print(f"Diagrammi salvati in: {salva}")
    plt.show()


# ---------------------------------------------------------------------------
# Helper disegno
# ---------------------------------------------------------------------------

def _disegna_vincoli(ax, vincoli: list, nodi: List[Nodo]) -> None:
    for v in vincoli:
        n = nodi[v.nodo]
        fissi = (not v.ux_fisso, not v.uy_fisso, not v.phi_fisso)

        if not v.ux_fisso and not v.uy_fisso and not v.phi_fisso:
            # incastro
            ax.plot(n.x, n.y, "s", color="#3B6D11",
                    markersize=14, zorder=1, alpha=0.5)
        elif not v.ux_fisso and not v.uy_fisso and v.phi_fisso:
            # cerniera
            triangle = plt.Polygon(
                [[n.x, n.y], [n.x - 0.15, n.y - 0.25], [n.x + 0.15, n.y - 0.25]],
                closed=True, color="#3B6D11", alpha=0.5, zorder=1
            )
            ax.add_patch(triangle)
            ax.plot(n.x, n.y, "o", color="#3B6D11", markersize=5, zorder=3)
        elif v.ux_fisso and not v.uy_fisso and v.phi_fisso:
            # carrello
            triangle = plt.Polygon(
                [[n.x, n.y], [n.x - 0.15, n.y - 0.25], [n.x + 0.15, n.y - 0.25]],
                closed=True, color="#3B6D11", alpha=0.5, zorder=1
            )
            ax.add_patch(triangle)
            ax.plot(n.x, n.y, "o", color="white",
                    markeredgecolor="#3B6D11", markersize=6, zorder=3)
        else:
            ax.plot(n.x, n.y, "s", color="#3B6D11",
                    markersize=14, zorder=1, alpha=0.5)


def _disegna_carichi(
    ax, carichi: List[Carico], travi: List[Trave], nodi: List[Nodo]
) -> None:
    span = max(
        max(n.x for n in nodi) - min(n.x for n in nodi),
        max(n.y for n in nodi) - min(n.y for n in nodi),
        0.5,
    )
    arrow_len = span * 0.12
    q_scale   = span * 0.08

    for c in carichi:
        if c.type == "nodo_fy":
            n = nodi[c.obj]
            ax.annotate("",
                xy=(n.x, n.y),
                xytext=(n.x, n.y + arrow_len),
                arrowprops=dict(arrowstyle="->", color="#D85A30", lw=1.5)
            )
            ax.text(n.x + 0.05, n.y + arrow_len * 0.5,
                    f"{c.val1:.1f} kN", fontsize=7, color="#D85A30")

        elif c.type == "nodo_fx":
            n = nodi[c.obj]
            ax.annotate("",
                xy=(n.x, n.y),
                xytext=(n.x - arrow_len, n.y),
                arrowprops=dict(arrowstyle="->", color="#D85A30", lw=1.5)
            )
            ax.text(n.x - arrow_len * 0.5, n.y + 0.05,
                    f"{c.val1:.1f} kN", fontsize=7, color="#D85A30")

        elif c.type in ("uniforme", "triangolare_sx", "triangolare_dx"):
            el = travi[c.obj]
            ni, nj = nodi[el.nodo_i], nodi[el.nodo_j]
            L, theta = _lunghezza_theta(el, nodi)
            perp = (-math.sin(theta), math.cos(theta))
            n_steps = 8
            for k in range(n_steps + 1):
                t = k / n_steps
                if c.type == "uniforme":
                    q = c.val1
                elif c.type == "triangolare_sx":
                    q = c.val1 * (1 - t)
                else:
                    q = c.val1 * t
                bx = ni.x + t * (nj.x - ni.x)
                by = ni.y + t * (nj.y - ni.y)
                ax.annotate("",
                    xy=(bx, by),
                    xytext=(bx + perp[0]*q*q_scale, by + perp[1]*q*q_scale),
                    arrowprops=dict(arrowstyle="->", color="#BA7517", lw=1.0)
                )
            ax.text(
                (ni.x+nj.x)/2 + perp[0]*c.val1*q_scale*1.3,
                (ni.y+nj.y)/2 + perp[1]*c.val1*q_scale*1.3,
                f"q={c.val1:.1f} kN/m", fontsize=7, color="#BA7517", ha="center"
            )


def _disegna_reazioni(ax, risultato: Risultato, nodi: List[Nodo]) -> None:
    span = max(
        max(n.x for n in nodi) - min(n.x for n in nodi),
        max(n.y for n in nodi) - min(n.y for n in nodi),
        0.5,
    )
    arrow_len = span * 0.15

    for label, (val, unit) in risultato.reazioni_vincolari().items():
        if abs(val) < 1e-6:
            continue
        nodo_idx = int(label.split("_N")[1]) - 1
        comp = label.split("_")[0]
        n = nodi[nodo_idx]

        if comp == "Ry":
            dy = np.sign(val) * arrow_len
            ax.annotate("",
                xy=(n.x, n.y + dy),
                xytext=(n.x, n.y),
                arrowprops=dict(arrowstyle="->", color="#993556", lw=2.0)
            )
            ax.text(n.x + 0.05, n.y + dy * 0.6,
                    f"{val:.2f} kN", fontsize=7, color="#993556", fontweight="bold")
        elif comp == "Rx":
            dx = np.sign(val) * arrow_len
            ax.annotate("",
                xy=(n.x + dx, n.y),
                xytext=(n.x, n.y),
                arrowprops=dict(arrowstyle="->", color="#993556", lw=2.0)
            )
            ax.text(n.x + dx * 0.5, n.y + 0.05,
                    f"{val:.2f} kN", fontsize=7, color="#993556", fontweight="bold")
