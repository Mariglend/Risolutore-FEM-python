"""
fem_travi
=========
Solver FEM per travi piane 2D (Euler-Bernoulli).

Esportazioni principali
-----------------------
    Struttura   – contenitore principale
    Nodo        – nodo con coordinate (x, y)
    Trave       – elemento trave (nodo_i, nodo_j, EI, EA, theta)
    Vincolo     – condizione al contorno su un nodo
    Carico      – carico nodale o distribuito
    Risultato   – output della soluzione (U, R, K, F)

    plot_struttura   – disegna struttura + deformata
    plot_diagrammi   – diagrammi M, V, N
"""

from .core import Carico, Nodo, Trave, Vincolo
from .solver import Risultato, Struttura
from .plotter import plot_struttura, plot_diagrammi

__all__ = [
    "Struttura",
    "Nodo",
    "Trave",
    "Vincolo",
    "Carico",
    "Risultato",
    "plot_struttura",
    "plot_diagrammi",
]

__version__ = "1.0.0"
__author__ = "fem_travi"
