"""
fem_travi
=========
Solver FEM per travi e telai piani 2D (Euler-Bernoulli).

    from fem_travi import Struttura
    s = Struttura()
    a, b = s.nodo(0, 0), s.nodo(6, 0)
    t = s.trave(a, b, EI=1e4)
    s.cerniera(a); s.carrello(b)
    s.carico_distribuito(t, -10)
    s.risolvi().stampa()

Interfaccia grafica:  ``fem-travi``  (oppure ``python -m fem_travi``)
"""

from .core import (
    Carico,
    CaricoDistribuito,
    CaricoTermico,
    CoppiaInCampata,
    CoppiaNodale,
    ForzaInCampata,
    ForzaNodale,
    Nodo,
    Trave,
    Vincolo,
)
from .solver import Risultato, Struttura, StrutturaLabile
from .io import carica_json, esporta_csv, risolvi_dict, struttura_da_dict


def plot_struttura(*a, **k):
    from .plotter import plot_struttura as f
    return f(*a, **k)


def plot_diagrammi(*a, **k):
    from .plotter import plot_diagrammi as f
    return f(*a, **k)


__all__ = [
    "Struttura", "Risultato", "StrutturaLabile",
    "Nodo", "Trave", "Vincolo",
    "ForzaNodale", "CoppiaNodale", "CaricoDistribuito", "ForzaInCampata",
    "CoppiaInCampata", "CaricoTermico", "Carico",
    "struttura_da_dict", "risolvi_dict", "carica_json", "esporta_csv",
    "plot_struttura", "plot_diagrammi",
]

__version__ = "2.0.0"
