"""
fem_travi/solver.py
===================
API di alto livello per il solver FEM.

Esempio d'uso
-------------
    from fem_travi import Struttura, Nodo, Trave, Vincolo, Carico

    s = Struttura()
    s.aggiungi_nodo(Nodo(0, 0))
    s.aggiungi_nodo(Nodo(6, 0))
    s.aggiungi_trave(Trave(0, 1, EI=10000, EA=100000))
    s.aggiungi_vincolo(Vincolo(0, ux_fisso=True, uy_fisso=True, phi_fisso=True))
    s.aggiungi_vincolo(Vincolo(1, ux_fisso=False, uy_fisso=True, phi_fisso=False))
    s.aggiungi_carico(Carico('uniforme', obj=0, val1=10.0))

    risultato = s.risolvi()
    risultato.stampa()
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from .assembler import (
    applica_BC,
    assembla_F,
    assembla_K,
    calcola_reazioni,
    risolvi as _risolvi,
)
from .core import Carico, Nodo, Trave, Vincolo


# ---------------------------------------------------------------------------
# Risultati
# ---------------------------------------------------------------------------

@dataclass
class Risultato:
    """Contenitore per i risultati della soluzione FEM."""

    U: np.ndarray
    """Vettore degli spostamenti nodali [m, m, rad] per ogni nodo."""

    R: np.ndarray
    """Vettore delle reazioni vincolari [kN, kN, kN·m] per ogni GDL."""

    K: np.ndarray
    """Matrice di rigidezza globale."""

    F: np.ndarray
    """Vettore dei carichi nodali equivalenti."""

    nodi: List[Nodo]
    vincoli: List[Vincolo]
    travi: List[Trave]

    # indici GDL
    free: List[int]
    fixed: List[int]

    def reazioni_vincolari(self) -> dict:
        """
        Dizionario delle reazioni vincolari significative.

        Reazione fisica = R[gi] = (K·U - F)[gi].
        Per i GDL vincolati F[gi]=0, quindi R[gi] = (K·U)[gi]
        = forza che la struttura trasmette al vincolo = reazione del vincolo.

        Convezione segno (positivo nel senso del GDL):
          Ry > 0 → verso l'alto
          Rx > 0 → verso destra
          Mz > 0 → antiorario (momento reazione che bilancia il carico)
        """
        comp = ["Rx", "Ry", "Mz"]
        unita = ["kN", "kN", "kN·m"]
        out = {}
        for v in self.vincoli:
            flags = [v.ux_fisso, v.uy_fisso, v.phi_fisso]
            for k, fisso in enumerate(flags):
                if fisso:
                    gi = v.nodo * 3 + k
                    label = f"{comp[k]}_N{v.nodo + 1}"
                    out[label] = (float(self.R[gi]), unita[k])
        return out

    def spostamenti_nodali(self) -> dict:
        """Dizionario degli spostamenti per ogni nodo."""
        comp = ["ux", "uy", "φ"]
        unita = ["m", "m", "rad"]
        out = {}
        for i in range(len(self.nodi)):
            for k in range(3):
                gi = i * 3 + k
                label = f"{comp[k]}_N{i + 1}"
                out[label] = (self.U[gi], unita[k])
        return out

    def verifica_equilibrio(self, tol: float = 1e-6) -> dict:
        """
        Verifica ΣFx = 0 e ΣFy = 0.

        Identità: R + F = K·U  →  sum(R + F) = sum(K·U)
        Per struttura in equilibrio: sum_x(R+F) = 0 e sum_y(R+F) = 0.
        Equivalente a: sum(reazioni_fisiche) = sum(carichi_nodali).
        """
        # R = K·U - F  →  R + F = K·U
        # Equilibrio globale: sum_x(K·U) = 0, sum_y(K·U) = 0
        KU = self.K @ self.U
        sfx = float(np.sum(KU[0::3]))
        sfy = float(np.sum(KU[1::3]))
        return {
            "ΣFx [kN]": sfx,
            "ΣFy [kN]": sfy,
            "ok": abs(sfx) < tol and abs(sfy) < tol,
        }

    def stampa(self, decimali: int = 4) -> None:
        """Stampa a video un report completo."""
        n_nodi = len(self.nodi)
        ndof = 3 * n_nodi
        print("=" * 60)
        print("  RISULTATI FEM — TRAVI PIANE 2D")
        print("=" * 60)
        print(f"  Nodi: {n_nodi}   GDL totali: {ndof}")
        print(f"  GDL liberi: {len(self.free)}   GDL vincolati: {len(self.fixed)}")
        hi = len(self.fixed) - 3
        print(f"  Grado di iperstaticità: {max(hi, 0)}")
        print()

        print("── Spostamenti nodali ──────────────────────────")
        for label, (val, unit) in self.spostamenti_nodali().items():
            print(f"  {label:12s} = {val:+.{decimali}e}  {unit}")
        print()

        print("── Reazioni vincolari ──────────────────────────")
        for label, (val, unit) in self.reazioni_vincolari().items():
            print(f"  {label:12s} = {val:+.{decimali}f}  {unit}")
        print()

        eq = self.verifica_equilibrio()
        ok_str = "OK ✓" if eq["ok"] else "ERRORE ✗"
        print("── Verifica equilibrio ─────────────────────────")
        print(f"  ΣFx = {eq['ΣFx [kN]']:+.6f} kN")
        print(f"  ΣFy = {eq['ΣFy [kN]']:+.6f} kN")
        print(f"  {ok_str}")
        print("=" * 60)

    def matrice_K_str(self, decimali: int = 2) -> str:
        """Restituisce la matrice K come stringa formattata."""
        n = self.K.shape[0]
        comp = ["ux", "uy", "φ"]
        labels = [f"{comp[i%3]}{i//3+1}" for i in range(n)]
        header = f"{'':8s}" + "".join(f"{l:>12s}" for l in labels)
        lines = [header]
        for i, row in enumerate(self.K):
            vals = "".join(f"{v:+12.{decimali}e}" for v in row)
            lines.append(f"{labels[i]:8s}{vals}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Struttura
# ---------------------------------------------------------------------------

class Struttura:
    """
    Contenitore principale della struttura.

    Gestisce nodi, travi, vincoli e carichi, espone il metodo `risolvi()`.
    """

    def __init__(self):
        self.nodi:    List[Nodo]    = []
        self.travi:   List[Trave]   = []
        self.vincoli: List[Vincolo] = []
        self.carichi: List[Carico]  = []

    # -- builder methods -------------------------------------------------------

    def aggiungi_nodo(self, nodo: Nodo) -> int:
        """Aggiunge un nodo e restituisce il suo indice."""
        nodo.id = len(self.nodi)
        self.nodi.append(nodo)
        return nodo.id

    def aggiungi_trave(self, trave: Trave) -> int:
        """Aggiunge una trave e restituisce il suo indice."""
        self._valida_trave(trave)
        trave.id = len(self.travi)
        self.travi.append(trave)
        return trave.id

    def aggiungi_vincolo(self, vincolo: Vincolo) -> None:
        """Aggiunge o sostituisce il vincolo su un nodo."""
        self._valida_nodo_idx(vincolo.nodo, "vincolo")
        self.vincoli = [v for v in self.vincoli if v.nodo != vincolo.nodo]
        self.vincoli.append(vincolo)

    def aggiungi_carico(self, carico: Carico) -> None:
        """Aggiunge un carico alla struttura."""
        self._valida_carico(carico)
        self.carichi.append(carico)

    # -- validazione -----------------------------------------------------------

    def _valida_nodo_idx(self, idx: int, contesto: str = "") -> None:
        if not (0 <= idx < len(self.nodi)):
            raise IndexError(
                f"Indice nodo {idx} non valido ({contesto}). "
                f"Nodi presenti: {len(self.nodi)}"
            )

    def _valida_trave(self, t: Trave) -> None:
        self._valida_nodo_idx(t.nodo_i, "trave nodo_i")
        self._valida_nodo_idx(t.nodo_j, "trave nodo_j")
        if t.nodo_i == t.nodo_j:
            raise ValueError("Trave con nodo_i == nodo_j.")
        if t.EI <= 0:
            raise ValueError(f"EI deve essere > 0, trovato {t.EI}.")
        if t.EA <= 0:
            raise ValueError(f"EA deve essere > 0, trovato {t.EA}.")

    def _valida_carico(self, c: Carico) -> None:
        tipi_nodo = {"nodo_fx", "nodo_fy", "nodo_m"}
        tipi_trave = {"uniforme", "triangolare_sx", "triangolare_dx"}
        if c.type in tipi_nodo:
            self._valida_nodo_idx(c.obj, f"carico {c.type}")
        elif c.type in tipi_trave:
            if not (0 <= c.obj < len(self.travi)):
                raise IndexError(
                    f"Indice trave {c.obj} non valido per carico {c.type}."
                )
        else:
            raise ValueError(f"Tipo carico sconosciuto: '{c.type}'")

    # -- solve -----------------------------------------------------------------

    def risolvi(self) -> Risultato:
        """
        Assembla e risolve il sistema FEM.

        Restituisce un oggetto `Risultato` con spostamenti, reazioni
        e la matrice K completa.

        Lancia
        ------
        ValueError  se la struttura è incompleta
        numpy.linalg.LinAlgError  se il sistema è singolare (struttura labile)
        """
        if len(self.nodi) < 2:
            raise ValueError("Servono almeno 2 nodi.")
        if len(self.travi) < 1:
            raise ValueError("Servono almeno 1 trave.")
        if len(self.vincoli) < 1:
            raise ValueError("Servono almeno 1 vincolo.")

        K = assembla_K(self.travi, self.nodi)
        F = assembla_F(self.carichi, self.travi, self.nodi)
        K_rid, F_rid, free, fixed = applica_BC(K, F, self.vincoli)

        if len(free) == 0:
            raise ValueError(
                "Struttura completamente bloccata: nessun GDL libero."
            )

        U_free = _risolvi(K_rid, F_rid)

        ndof = K.shape[0]
        U = np.zeros(ndof)
        for li, gi in enumerate(free):
            U[gi] = U_free[li]

        R = calcola_reazioni(K, U, F)

        return Risultato(
            U=U, R=R, K=K, F=F,
            nodi=self.nodi,
            vincoli=self.vincoli,
            travi=self.travi,
            free=free,
            fixed=fixed,
        )

    # -- informazioni ----------------------------------------------------------

    def info(self) -> str:
        lines = [
            f"Struttura FEM 2D",
            f"  Nodi:    {len(self.nodi)}",
            f"  Travi:   {len(self.travi)}",
            f"  Vincoli: {len(self.vincoli)} nodi vincolati",
            f"  Carichi: {len(self.carichi)}",
            f"  GDL totali: {3 * len(self.nodi)}",
        ]
        return "\n".join(lines)
