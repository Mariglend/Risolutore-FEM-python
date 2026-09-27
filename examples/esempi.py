"""
examples/esempi.py
==================
Esempi pronti all'uso per testare il solver.

Esegui un esempio specifico:
    python examples/esempi.py --esempio 1
    python examples/esempi.py --esempio 2
    python examples/esempi.py --esempio 3
    python examples/esempi.py --esempio 4
    python examples/esempi.py --tutti
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from fem_travi import Struttura, Nodo, Trave, Vincolo, Carico
from fem_travi import plot_struttura, plot_diagrammi


# ---------------------------------------------------------------------------
# Esempio 1: trave semplicemente appoggiata con carico uniforme
# ---------------------------------------------------------------------------

def esempio_1():
    """
    Trave isostatica appoggiata agli estremi con carico uniforme q = 10 kN/m.

    N1 ─────────────── N2
    (cerniera)         (carrello)
         q = 10 kN/m

    Soluzione analitica:
        Ry_N1 = Ry_N2 = q*L/2 = 30 kN
        M_max = q*L²/8 = 45 kNm  (centro)
    """
    print("\n" + "="*60)
    print("  Esempio 1: trave appoggiata + carico uniforme")
    print("="*60)

    s = Struttura()
    s.aggiungi_nodo(Nodo(0, 0))   # N1
    s.aggiungi_nodo(Nodo(6, 0))   # N2

    s.aggiungi_trave(Trave(0, 1, EI=10_000, EA=100_000))

    s.aggiungi_vincolo(Vincolo(0, ux_fisso=True,  uy_fisso=True, phi_fisso=False))  # cerniera
    s.aggiungi_vincolo(Vincolo(1, ux_fisso=False, uy_fisso=True, phi_fisso=False))  # carrello

    s.aggiungi_carico(Carico("uniforme", obj=0, val1=10.0))

    print(s.info())
    r = s.risolvi()
    r.stampa()

    print(f"\n  Soluzione analitica: Ry = {10*6/2:.1f} kN  M_max = {10*6**2/8:.1f} kNm")

    plot_struttura(s, r, titolo="Esempio 1 — trave appoggiata + q uniforme")
    plot_diagrammi(s, r)
    return r


# ---------------------------------------------------------------------------
# Esempio 2: trave a sbalzo con forza puntuale
# ---------------------------------------------------------------------------

def esempio_2():
    """
    Trave a sbalzo (incastro a sinistra) con forza puntuale in mezzeria
    e forza concentrata all'estremo.

    N1 ═══════════════ N2 ──────────────── N3
    (incastro)         F=20kN ↓             F=10kN ↓

    L1 = L2 = 3 m
    """
    print("\n" + "="*60)
    print("  Esempio 2: trave a sbalzo con forze puntuali")
    print("="*60)

    s = Struttura()
    s.aggiungi_nodo(Nodo(0, 0))   # N1 — incastro
    s.aggiungi_nodo(Nodo(3, 0))   # N2 — nodo intermedio
    s.aggiungi_nodo(Nodo(6, 0))   # N3 — estremo libero

    s.aggiungi_trave(Trave(0, 1, EI=15_000, EA=120_000))  # T1
    s.aggiungi_trave(Trave(1, 2, EI=15_000, EA=120_000))  # T2

    s.aggiungi_vincolo(Vincolo(0, ux_fisso=True, uy_fisso=True, phi_fisso=True))  # incastro

    s.aggiungi_carico(Carico("nodo_fy", obj=1, val1=20.0))   # 20 kN su N2
    s.aggiungi_carico(Carico("nodo_fy", obj=2, val1=10.0))   # 10 kN su N3

    print(s.info())
    r = s.risolvi()
    r.stampa()

    plot_struttura(s, r, titolo="Esempio 2 — trave a sbalzo + forze puntuali")
    plot_diagrammi(s, r)
    return r


# ---------------------------------------------------------------------------
# Esempio 3: telaio a L (trave + colonna)
# ---------------------------------------------------------------------------

def esempio_3():
    """
    Telaio piano a L: colonna verticale + trave orizzontale.

            N3 ──────── N2
            |
            |  F=15kN ←
            |
            N1
         (incastro)

    Colonna: N1(0,0) → N3(0,4)
    Trave:   N3(0,4) → N2(5,4)
    """
    print("\n" + "="*60)
    print("  Esempio 3: telaio a L (colonna + trave)")
    print("="*60)

    s = Struttura()
    s.aggiungi_nodo(Nodo(0, 0))   # N1 — base colonna
    s.aggiungi_nodo(Nodo(5, 4))   # N2 — estremo trave
    s.aggiungi_nodo(Nodo(0, 4))   # N3 — nodo di giunzione

    s.aggiungi_trave(Trave(0, 2, EI=20_000, EA=150_000))  # T1 — colonna
    s.aggiungi_trave(Trave(2, 1, EI=20_000, EA=150_000))  # T2 — trave orizzontale

    s.aggiungi_vincolo(Vincolo(0, ux_fisso=True, uy_fisso=True, phi_fisso=True))  # incastro N1
    s.aggiungi_vincolo(Vincolo(1, ux_fisso=False, uy_fisso=True, phi_fisso=False)) # carrello N2

    s.aggiungi_carico(Carico("nodo_fx", obj=2, val1=15.0))  # forza orizzontale su N3
    s.aggiungi_carico(Carico("uniforme", obj=1, val1=8.0))  # q su trave orizzontale

    print(s.info())
    r = s.risolvi()
    r.stampa()

    plot_struttura(s, r, titolo="Esempio 3 — telaio a L")
    plot_diagrammi(s, r)
    return r


# ---------------------------------------------------------------------------
# Esempio 4: trave su due campate con carico triangolare
# ---------------------------------------------------------------------------

def esempio_4():
    """
    Trave continua su tre appoggi con carico triangolare sulla seconda campata.

    N1 ──────────── N2 ──────────── N3
    (cerniera)    (cerniera)      (carrello)
                  ↑ q triangolare ↑
                  max a sx = 12 kN/m

    L1 = L2 = 4 m
    """
    print("\n" + "="*60)
    print("  Esempio 4: trave continua + carico triangolare")
    print("="*60)

    s = Struttura()
    s.aggiungi_nodo(Nodo(0, 0))   # N1
    s.aggiungi_nodo(Nodo(4, 0))   # N2
    s.aggiungi_nodo(Nodo(8, 0))   # N3

    s.aggiungi_trave(Trave(0, 1, EI=8_000, EA=80_000))   # T1 — prima campata
    s.aggiungi_trave(Trave(1, 2, EI=8_000, EA=80_000))   # T2 — seconda campata

    s.aggiungi_vincolo(Vincolo(0, ux_fisso=True,  uy_fisso=True, phi_fisso=False))  # cerniera
    s.aggiungi_vincolo(Vincolo(1, ux_fisso=False, uy_fisso=True, phi_fisso=False))  # cerniera interna
    s.aggiungi_vincolo(Vincolo(2, ux_fisso=False, uy_fisso=True, phi_fisso=False))  # carrello

    # carico triangolare sulla seconda campata (T2), max a sinistra
    s.aggiungi_carico(Carico("triangolare_sx", obj=1, val1=12.0))

    print(s.info())
    r = s.risolvi()
    r.stampa()

    plot_struttura(s, r, titolo="Esempio 4 — trave continua + carico triangolare")
    plot_diagrammi(s, r)
    return r


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Esempi solver FEM travi")
    parser.add_argument("--esempio", type=int, choices=[1, 2, 3, 4],
                        help="Numero dell'esempio da eseguire")
    parser.add_argument("--tutti", action="store_true",
                        help="Esegui tutti gli esempi")
    parser.add_argument("--no-plot", action="store_true",
                        help="Disabilita i grafici (solo output testo)")
    args = parser.parse_args()

    if args.no_plot:
        import fem_travi.plotter as _p
        _p.plot_struttura = lambda *a, **k: None
        _p.plot_diagrammi = lambda *a, **k: None
        plot_struttura = lambda *a, **k: None
        plot_diagrammi = lambda *a, **k: None

    esegui = {
        1: esempio_1,
        2: esempio_2,
        3: esempio_3,
        4: esempio_4,
    }

    if args.tutti:
        for fn in esegui.values():
            fn()
    elif args.esempio:
        esegui[args.esempio]()
    else:
        parser.print_help()
        print("\nEsempio rapido (--no-plot per saltare i grafici):")
        print("  python examples/esempi.py --esempio 1")
