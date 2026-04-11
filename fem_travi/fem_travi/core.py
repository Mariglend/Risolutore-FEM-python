"""
fem_travi/core.py
=================
Modello FEM per travi piane 2D (Euler-Bernoulli).

GDL per nodo: [ux, uy, φ]  →  ndof = 3 * n_nodi
Elemento trave: 6 GDL locali, matrice di rigidezza 6×6,
rotazione con matrice T per travi inclinate di θ gradi.

Convenzioni di segno
--------------------
- Forze nodali positive → verso l'alto (uy+) o verso destra (ux+)
- Carichi distribuiti positivi → verso l'alto (perpendicolare all'asse locale)
- Momenti positivi → antiorari (φ+)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Dataclass di input
# ---------------------------------------------------------------------------

@dataclass
class Nodo:
    """Nodo della struttura con coordinate globali (m)."""
    x: float
    y: float
    id: int = field(default=-1, repr=False)

    def __post_init__(self):
        self.x = float(self.x)
        self.y = float(self.y)


@dataclass
class Trave:
    """
    Elemento trave Euler-Bernoulli.

    Parametri
    ----------
    nodo_i, nodo_j : int
        Indici dei nodi terminali.
    EI : float
        Rigidezza flessionale [kN·m²].
    EA : float
        Rigidezza assiale [kN].
    theta : float | None
        Angolo di inclinazione [rad]. Se None viene calcolato
        automaticamente dai nodi.
    """
    nodo_i: int
    nodo_j: int
    EI: float
    EA: float
    theta: Optional[float] = None
    id: int = field(default=-1, repr=False)


@dataclass
class Vincolo:
    """
    Condizioni al contorno su un nodo.

    ux_fisso, uy_fisso, phi_fisso : bool
        True  → GDL bloccato (spostamento imposto = 0)
        False → GDL libero
    Combinazioni tipiche:
        incastro  → (True, True, True)
        cerniera  → (True, True, False)
        carrello  → (False, True, False)
    """
    nodo: int
    ux_fisso: bool = True
    uy_fisso: bool = True
    phi_fisso: bool = True


@dataclass
class Carico:
    """
    Carico applicato alla struttura.

    type : str
        'nodo_fx'        – forza orizzontale su nodo [kN]
        'nodo_fy'        – forza verticale su nodo [kN]  (↓ positivo)
        'nodo_m'         – momento su nodo [kN·m]  (antiorario positivo)
        'uniforme'       – carico distribuito uniforme su trave [kN/m]
        'triangolare_sx' – triangolare con massimo all'inizio [kN/m]
        'triangolare_dx' – triangolare con massimo alla fine  [kN/m]
    obj : int
        Indice nodo (per carichi nodali) o indice trave (per distribuiti).
    val1 : float
        Valore principale del carico.
    val2 : float
        Secondo valore (usato solo per 'uniforme' se si vuole trapezioidale
        in futuro; ignorato dagli altri tipi).
    """
    type: str
    obj: int
    val1: float
    val2: float = 0.0
