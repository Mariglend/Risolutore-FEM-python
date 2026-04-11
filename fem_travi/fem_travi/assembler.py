"""
fem_travi/assembler.py
======================
Assemblaggio e soluzione del sistema FEM.

Flusso:
    1. k_locale(el)          → K_loc 6×6
    2. T(theta)               → matrice rotazione 6×6
    3. k_globale(el, nodi)   → K_glob = Tᵀ K_loc T
    4. assembla_K(travi, nodi) → K_globale n×n
    5. assembla_F(carichi, travi, nodi) → F_globale n
    6. applica_BC(K, F, vincoli) → K_ridotta, F_ridotta
    7. risolvi(K_rid, F_rid)  → U_liberi
    8. calcola_reazioni(K, U, F) → R
"""

from __future__ import annotations

import math
from typing import List, Tuple

import numpy as np

from .core import Carico, Nodo, Trave, Vincolo


# ---------------------------------------------------------------------------
# Geometria elemento
# ---------------------------------------------------------------------------

def _lunghezza_theta(el: Trave, nodi: List[Nodo]) -> Tuple[float, float]:
    """Restituisce (L, theta) per l'elemento, calcolando theta se non fornito."""
    ni, nj = nodi[el.nodo_i], nodi[el.nodo_j]
    dx, dy = nj.x - ni.x, nj.y - ni.y
    L = math.hypot(dx, dy)
    if L < 1e-12:
        raise ValueError(
            f"Trave {el.id}: i nodi {el.nodo_i} e {el.nodo_j} coincidono."
        )
    theta = el.theta if el.theta is not None else math.atan2(dy, dx)
    return L, theta


# ---------------------------------------------------------------------------
# Matrici elemento
# ---------------------------------------------------------------------------

def k_locale(L: float, EI: float, EA: float) -> np.ndarray:
    """
    Matrice di rigidezza locale 6×6 per elemento trave Euler-Bernoulli.

    GDL ordine: [u_i, v_i, φ_i, u_j, v_j, φ_j]  (sistema locale)

    References
    ----------
    Cook et al., "Concepts and Applications of Finite Element Analysis",
    4th ed., Chapter 4.
    """
    a = EA / L
    b = 12 * EI / L**3
    c = 6  * EI / L**2
    d = 4  * EI / L
    e = 2  * EI / L

    K = np.zeros((6, 6))
    # Contributo assiale
    K[0, 0] =  a;  K[0, 3] = -a
    K[3, 0] = -a;  K[3, 3] =  a
    # Contributo flessionale
    K[1, 1] =  b;  K[1, 2] =  c;  K[1, 4] = -b;  K[1, 5] =  c
    K[2, 1] =  c;  K[2, 2] =  d;  K[2, 4] = -c;  K[2, 5] =  e
    K[4, 1] = -b;  K[4, 2] = -c;  K[4, 4] =  b;  K[4, 5] = -c
    K[5, 1] =  c;  K[5, 2] =  e;  K[5, 4] = -c;  K[5, 5] =  d
    return K


def matrice_rotazione(theta: float) -> np.ndarray:
    """
    Matrice di rotazione 6×6 da sistema locale a globale.

    T tale che:  u_loc = T · u_glob
    K_glob = Tᵀ · K_loc · T
    """
    c, s = math.cos(theta), math.sin(theta)
    T = np.zeros((6, 6))
    T[0, 0] =  c;  T[0, 1] =  s
    T[1, 0] = -s;  T[1, 1] =  c
    T[2, 2] =  1.0
    T[3, 3] =  c;  T[3, 4] =  s
    T[4, 3] = -s;  T[4, 4] =  c
    T[5, 5] =  1.0
    return T


def k_globale(el: Trave, nodi: List[Nodo]) -> Tuple[np.ndarray, List[int]]:
    """
    Restituisce (K_glob 6×6, dofs) per l'elemento.

    dofs : lista dei 6 GDL globali dell'elemento
           [3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2]
    """
    L, theta = _lunghezza_theta(el, nodi)
    Kl = k_locale(L, el.EI, el.EA)
    T  = matrice_rotazione(theta)
    Kg = T.T @ Kl @ T
    i, j = el.nodo_i, el.nodo_j
    dofs = [3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2]
    return Kg, dofs


# ---------------------------------------------------------------------------
# Assemblaggio K globale
# ---------------------------------------------------------------------------

def assembla_K(travi: List[Trave], nodi: List[Nodo]) -> np.ndarray:
    """
    Assembla la matrice di rigidezza globale K (ndof × ndof).

    ndof = 3 * len(nodi)
    """
    ndof = 3 * len(nodi)
    K = np.zeros((ndof, ndof))
    for el in travi:
        Kg, dofs = k_globale(el, nodi)
        for a, da in enumerate(dofs):
            for b, db in enumerate(dofs):
                K[da, db] += Kg[a, b]
    return K


# ---------------------------------------------------------------------------
# Forze nodali equivalenti per carichi distribuiti
# ---------------------------------------------------------------------------

def _forze_equivalenti_locali(c: Carico, L: float) -> np.ndarray:
    """
    Vettore 6×1 delle forze nodali equivalenti in coordinate LOCALI
    per un carico distribuito sulla trave.

    Formule da Euler-Bernoulli (forze nella direzione trasversale locale v).
    """
    fne = np.zeros(6)
    q = c.val1

    if c.type == "uniforme":
        # val > 0 = verso il basso; forze nodali equiv. negative (verso il basso)
        fne[1] = -q * L / 2
        fne[2] = -q * L**2 / 12
        fne[4] = -q * L / 2
        fne[5] =  q * L**2 / 12

    elif c.type == "triangolare_sx":
        # massimo q all'inizio (nodo i), zero alla fine (nodo j)
        fne[1] = -7 * q * L / 20
        fne[2] = -   q * L**2 / 20
        fne[4] = -3 * q * L / 20
        fne[5] =     q * L**2 / 30

    elif c.type == "triangolare_dx":
        # zero all'inizio (nodo i), massimo q alla fine (nodo j)
        fne[1] = -3 * q * L / 20
        fne[2] = -   q * L**2 / 30
        fne[4] = -7 * q * L / 20
        fne[5] =     q * L**2 / 20

    return fne


# ---------------------------------------------------------------------------
# Assemblaggio vettore F
# ---------------------------------------------------------------------------

def assembla_F(
    carichi: List[Carico],
    travi: List[Trave],
    nodi: List[Nodo],
) -> np.ndarray:
    """
    Assembla il vettore dei carichi nodali equivalenti F (ndof,).

    Convenzione carichi (coerente con ingegneria strutturale):
    nodo_fy: val > 0 → verso il basso  (F[uy] -= val)
    nodo_fx: val > 0 → verso sinistra  (F[ux] -= val)
    nodo_m:  val > 0 → orario           (F[phi] -= val)
    Distribuiti: val > 0 → verso l'alto (perp. asse locale, forze equiv.)
    """
    ndof = 3 * len(nodi)
    F = np.zeros(ndof)

    for c in carichi:
        if c.type == "nodo_fy":
            F[c.obj * 3 + 1] -= c.val1   # val>0 = verso il basso
        elif c.type == "nodo_fx":
            F[c.obj * 3]     -= c.val1   # val>0 = verso sinistra
        elif c.type == "nodo_m":
            F[c.obj * 3 + 2] -= c.val1   # val>0 = orario
        else:
            # carico distribuito su trave
            el = travi[c.obj]
            L, theta = _lunghezza_theta(el, nodi)
            fne_loc = _forze_equivalenti_locali(c, L)
            T = matrice_rotazione(theta)
            fne_glob = T.T @ fne_loc          # rotazione al sistema globale
            i, j = el.nodo_i, el.nodo_j
            dofs = [3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2]
            for a, da in enumerate(dofs):
                F[da] += fne_glob[a]

    return F


# ---------------------------------------------------------------------------
# Condizioni al contorno
# ---------------------------------------------------------------------------

def applica_BC(
    K: np.ndarray,
    F: np.ndarray,
    vincoli: List[Vincolo],
) -> Tuple[np.ndarray, np.ndarray, List[int], List[int]]:
    """
    Elimina i GDL vincolati e restituisce il sistema ridotto.

    Restituisce
    -----------
    K_rid : ndarray (nf × nf)
    F_rid : ndarray (nf,)
    free  : lista degli indici GDL liberi
    fixed : lista degli indici GDL vincolati
    """
    ndof = K.shape[0]
    is_fixed = np.zeros(ndof, dtype=bool)

    for v in vincoli:
        n = v.nodo
        if v.ux_fisso:  is_fixed[3*n]   = True
        if v.uy_fisso:  is_fixed[3*n+1] = True
        if v.phi_fisso: is_fixed[3*n+2] = True

    free  = [i for i in range(ndof) if not is_fixed[i]]
    fixed = [i for i in range(ndof) if     is_fixed[i]]

    K_rid = K[np.ix_(free, free)]
    F_rid = F[free]
    return K_rid, F_rid, free, fixed


# ---------------------------------------------------------------------------
# Soluzione
# ---------------------------------------------------------------------------

def risolvi(
    K_rid: np.ndarray,
    F_rid: np.ndarray,
) -> np.ndarray:
    """
    Risolve K_rid · U_f = F_rid con numpy.linalg.solve (LU con pivot).

    Lancia numpy.linalg.LinAlgError se la matrice è singolare.
    """
    if K_rid.size == 0:
        return np.array([])
    return np.linalg.solve(K_rid, F_rid)


def calcola_reazioni(
    K: np.ndarray,
    U: np.ndarray,
    F: np.ndarray,
) -> np.ndarray:
    """
    Reazioni vincolari: R = K · U − F

    Solo i componenti corrispondenti ai GDL vincolati sono significativi.
    """
    return K @ U - F
