"""
examples/esempi_v2.py
=====================
Gli stessi casi da esame scritti con l'API v2 (convenzione: y verso l'alto,
carichi verso il basso negativi).

    python examples/esempi_v2.py            # stampa e grafici
    python examples/esempi_v2.py --no-plot  # solo testo
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from fem_travi import Struttura, plot_diagrammi, plot_struttura


def trave_gerber():
    """Incastro – cerniera interna – carrello, q sulla seconda campata."""
    s = Struttura()
    a, b, c = s.nodo(0, 0), s.nodo(3, 0), s.nodo(6, 0)
    s.trave(a, b)
    t2 = s.trave(b, c)
    s.cerniera_interna(b)
    s.incastro(a)
    s.carrello(c)
    s.carico_distribuito(t2, -10)            # 10 kN/m verso il basso
    return s


def portale_con_termico():
    """Portale incastrato con trave caldata a farfalla e forza orizzontale."""
    s = Struttura()
    n = [s.nodo(0, 0), s.nodo(0, 4), s.nodo(6, 4), s.nodo(6, 0)]
    s.trave(n[0], n[1], EI=2e4)
    t = s.trave(n[1], n[2], EI=3e4)
    s.trave(n[2], n[3], EI=2e4)
    s.incastro(n[0])
    s.incastro(n[3])
    s.forza(n[1], Fx=15)
    s.carico_termico(t, alpha=1.2e-5, dT_farfalla=20, h=0.4)
    return s


def arco_tre_cerniere():
    """Arco a tre cerniere con carico riferito alla proiezione orizzontale."""
    s = Struttura()
    a, c, b = s.nodo(0, 0), s.nodo(4, 2), s.nodo(8, 0)
    t1, t2 = s.trave(a, c), s.trave(c, b)
    s.cerniera(a)
    s.cerniera(b)
    s.cerniera_interna(c)
    for t in (t1, t2):
        s.carico_distribuito(t, -10, proiezione=True)
    return s


def trave_su_carrello_inclinato():
    """Trave appoggiata con carrello su piano inclinato e forza in campata."""
    s = Struttura()
    a, b = s.nodo(0, 0), s.nodo(6, 0)
    t = s.trave(a, b)
    s.cerniera(a)
    s.carrello(b, angolo=30)
    s.forza_in_campata(t, a=2.0, Fy=-20)
    return s


if __name__ == "__main__":
    grafici = "--no-plot" not in sys.argv
    for f in (trave_gerber, portale_con_termico, arco_tre_cerniere, trave_su_carrello_inclinato):
        print(f"\n### {f.__doc__.strip()}")
        s = f()
        r = s.risolvi()
        r.stampa()
        if grafici:
            plot_struttura(s, r, titolo=f.__name__)
            plot_diagrammi(s, r)
