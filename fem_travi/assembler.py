"""
fem_travi/assembler.py
======================
Calcoli a livello di elemento (sistema locale):

- matrice di rigidezza e di rotazione
- forze nodali equivalenti dei carichi in campata (esatte per Euler-Bernoulli)
- condensazione statica degli svincoli (cerniere interne, aste)
- sollecitazioni N, V, M e deformata lungo l'asse, per equilibrio a partire
  dalle forze di estremità: il diagramma del momento sotto carico distribuito
  è quindi la parabola esatta, non l'interpolazione lineare delle funzioni
  di forma.

GDL locali: [u_i, v_i, φ_i, u_j, v_j, φ_j]
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np

from .core import Nodo, Trave

# Gauss-Legendre a 5 punti su [0, 1]: esatto fino al grado 9
_GP, _GW = np.polynomial.legendre.leggauss(5)
_GP = 0.5 * (_GP + 1.0)
_GW = 0.5 * _GW


# ---------------------------------------------------------------------------
# Geometria
# ---------------------------------------------------------------------------

def _lunghezza_theta(el: Trave, nodi: List[Nodo]) -> Tuple[float, float]:
    ni, nj = nodi[el.nodo_i], nodi[el.nodo_j]
    dx, dy = nj.x - ni.x, nj.y - ni.y
    L = math.hypot(dx, dy)
    if L < 1e-12:
        raise ValueError(f"Trave {el.id + 1}: i nodi {el.nodo_i + 1} e {el.nodo_j + 1} coincidono.")
    return L, math.atan2(dy, dx)


# ---------------------------------------------------------------------------
# Matrici elemento
# ---------------------------------------------------------------------------

def k_locale(L: float, EI: float, EA: float) -> np.ndarray:
    """Matrice di rigidezza locale 6×6 (Euler-Bernoulli)."""
    a = EA / L
    b = 12 * EI / L**3
    c = 6 * EI / L**2
    d = 4 * EI / L
    e = 2 * EI / L
    return np.array([
        [ a,  0,  0, -a,  0,  0],
        [ 0,  b,  c,  0, -b,  c],
        [ 0,  c,  d,  0, -c,  e],
        [-a,  0,  0,  a,  0,  0],
        [ 0, -b, -c,  0,  b, -c],
        [ 0,  c,  e,  0, -c,  d],
    ], dtype=float)


def matrice_rotazione(theta: float) -> np.ndarray:
    """T tale che u_loc = T · u_glob."""
    c, s = math.cos(theta), math.sin(theta)
    T = np.zeros((6, 6))
    T[0, 0] = c;  T[0, 1] = s
    T[1, 0] = -s; T[1, 1] = c
    T[2, 2] = 1.0
    T[3, 3] = c;  T[3, 4] = s
    T[4, 3] = -s; T[4, 4] = c
    T[5, 5] = 1.0
    return T


def _hermite(xi: float, L: float):
    """Funzioni di forma flessionali e loro derivate prime (rispetto a x)."""
    N = np.array([
        1 - 3 * xi**2 + 2 * xi**3,
        L * (xi - 2 * xi**2 + xi**3),
        3 * xi**2 - 2 * xi**3,
        L * (-xi**2 + xi**3),
    ])
    dN = np.array([
        (-6 * xi + 6 * xi**2) / L,
        1 - 4 * xi + 3 * xi**2,
        (6 * xi - 6 * xi**2) / L,
        -2 * xi + 3 * xi**2,
    ])
    return N, dN


# ---------------------------------------------------------------------------
# Carichi in campata nel sistema locale
# ---------------------------------------------------------------------------

@dataclass
class CarichiLocali:
    """Tutti i carichi agenti su un elemento, già in coordinate locali."""
    # distribuiti lineari: (pxi, pxj, pyi, pyj)
    distribuiti: List[Tuple[float, float, float, float]] = field(default_factory=list)
    # forze concentrate: (a, Px, Py)
    forze: List[Tuple[float, float, float]] = field(default_factory=list)
    # coppie concentrate: (a, M)
    coppie: List[Tuple[float, float]] = field(default_factory=list)
    eps0: float = 0.0     # deformazione assiale impressa (termica)
    kappa0: float = 0.0   # curvatura impressa (termica a farfalla)

    def vuoto(self) -> bool:
        return not (self.distribuiti or self.forze or self.coppie
                    or self.eps0 or self.kappa0)

    def ha_trasversali(self) -> bool:
        return (any(abs(d[2]) + abs(d[3]) > 0 for d in self.distribuiti)
                or any(abs(f[2]) > 0 for f in self.forze)
                or bool(self.coppie) or bool(self.kappa0))


def forze_equivalenti(L: float, EI: float, EA: float, cl: CarichiLocali) -> np.ndarray:
    """Forze nodali equivalenti (locali, trave doppiamente incastrata)."""
    f = np.zeros(6)
    for pxi, pxj, pyi, pyj in cl.distribuiti:
        for g, w in zip(_GP, _GW):
            px = pxi + (pxj - pxi) * g
            py = pyi + (pyj - pyi) * g
            N, _ = _hermite(g, L)
            f[0] += w * L * (1 - g) * px
            f[3] += w * L * g * px
            f[[1, 2, 4, 5]] += w * L * N * py
    for a, Px, Py in cl.forze:
        xi = a / L
        N, _ = _hermite(xi, L)
        f[0] += (1 - xi) * Px
        f[3] += xi * Px
        f[[1, 2, 4, 5]] += N * Py
    for a, M0 in cl.coppie:
        _, dN = _hermite(a / L, L)
        f[[1, 2, 4, 5]] += dN * M0
    if cl.eps0:
        f[0] -= EA * cl.eps0
        f[3] += EA * cl.eps0
    if cl.kappa0:
        f[2] -= EI * cl.kappa0
        f[5] += EI * cl.kappa0
    return f


# ---------------------------------------------------------------------------
# Condensazione statica degli svincoli
# ---------------------------------------------------------------------------

def gdl_svincolati(el: Trave, nodi: List[Nodo]) -> List[int]:
    rel = []
    if el.asta or el.cerniera_i or nodi[el.nodo_i].cerniera:
        rel.append(2)
    if el.asta or el.cerniera_j or nodi[el.nodo_j].cerniera:
        rel.append(5)
    return rel


def condensa(K: np.ndarray, f: np.ndarray, rel: List[int]):
    """Condensa i GDL ``rel``: restituisce (Kc, fc) 6×6 / 6 con zeri su rel."""
    if not rel:
        return K.copy(), f.copy()
    ret = [k for k in range(6) if k not in rel]
    Krr = K[np.ix_(rel, rel)]
    Kar = K[np.ix_(ret, rel)]
    X = np.linalg.solve(Krr, Kar.T)          # Krr⁻¹ Kra
    Kc = np.zeros((6, 6)); fc = np.zeros(6)
    Kc[np.ix_(ret, ret)] = K[np.ix_(ret, ret)] - Kar @ X
    fc[ret] = f[ret] - X.T @ f[rel]
    return Kc, fc


def recupera_svincolati(K: np.ndarray, f: np.ndarray, u: np.ndarray, rel: List[int]) -> np.ndarray:
    """Ricava le rotazioni di estremità svincolate (momento nullo all'estremo)."""
    u = u.copy()
    if rel:
        ret = [k for k in range(6) if k not in rel]
        Krr = K[np.ix_(rel, rel)]
        u[rel] = np.linalg.solve(Krr, f[rel] - K[np.ix_(rel, ret)] @ u[ret])
    return u


# ---------------------------------------------------------------------------
# Sollecitazioni e deformata lungo l'asse
# ---------------------------------------------------------------------------

def griglia(L: float, cl: CarichiLocali, n: int = 41):
    """Ascisse di campionamento; le discontinuità compaiono due volte (sx, dx)."""
    salti = sorted({a for a, *_ in cl.forze} | {a for a, _ in cl.coppie})
    eps = 1e-9 * L
    xs = [x for x in np.linspace(0.0, L, n)
          if not any(abs(x - a) < 1e-7 * L for a in salti if eps < a < L - eps)]
    lati = [0] * len(xs)
    for a in salti:
        if eps < a < L - eps:
            xs += [a, a]; lati += [-1, 1]
    ordine = sorted(range(len(xs)), key=lambda k: (xs[k], lati[k]))
    x = np.array([xs[k] for k in ordine])
    x_eval = np.array([xs[k] + lati[k] * eps for k in ordine])
    return x, x_eval


def sollecitazioni(L: float, f_end: np.ndarray, cl: CarichiLocali, x_eval: np.ndarray):
    """N, V, M per equilibrio del tratto [0, x] (f_end = forze dei nodi sull'elemento)."""
    X1, Y1, M1 = f_end[0], f_end[1], f_end[2]
    x = np.asarray(x_eval, dtype=float)
    N = -X1 * np.ones_like(x)
    V = Y1 * np.ones_like(x)
    M = -M1 + Y1 * x
    for pxi, pxj, pyi, pyj in cl.distribuiti:
        dpx, dpy = (pxj - pxi) / L, (pyj - pyi) / L
        N -= pxi * x + dpx * x**2 / 2
        V += pyi * x + dpy * x**2 / 2
        M += pyi * x**2 / 2 + dpy * x**3 / 6
    for a, Px, Py in cl.forze:
        h = (x > a).astype(float)
        N -= Px * h
        V += Py * h
        M += Py * (x - a) * h
    for a, M0 in cl.coppie:
        M -= M0 * (x > a)
    return N, V, M


def _cumtrapz(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


def deformata_locale(L, EI, EA, u_loc, f_end, cl: CarichiLocali, n: int = 61):
    """Spostamenti locali (u, v) lungo l'asse, integrando curvatura e deformazione assiale.

    L'integrazione avviene su una griglia fine (errore relativo < 1e-6) e il
    risultato è poi campionato in ``n`` punti.
    """
    xo, _ = griglia(L, cl, n)
    x, xe = griglia(L, cl, 2001)
    N, _, M = sollecitazioni(L, f_end, cl, xe)
    kappa = M / EI + cl.kappa0
    eps = N / EA + cl.eps0
    phi = u_loc[2] + _cumtrapz(kappa, x)
    v = u_loc[1] + _cumtrapz(phi, x)
    u = u_loc[0] + _cumtrapz(eps, x)
    # correzione dell'errore di quadratura per chiudere esattamente sul nodo j
    v += (u_loc[4] - v[-1]) * x / L
    u += (u_loc[3] - u[-1]) * x / L
    return xo, np.interp(xo, x, u), np.interp(xo, x, v)


# ---------------------------------------------------------------------------
# Compatibilità v1
# ---------------------------------------------------------------------------

def k_globale(el: Trave, nodi: List[Nodo]):
    L, theta = _lunghezza_theta(el, nodi)
    T = matrice_rotazione(theta)
    i, j = el.nodo_i, el.nodo_j
    return T.T @ k_locale(L, el.EI, el.EA) @ T, [3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2]
