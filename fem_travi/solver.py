"""
fem_travi/solver.py
===================
API di alto livello.

Esempio
-------
    from fem_travi import Struttura

    s = Struttura()
    a = s.nodo(0, 0)
    b = s.nodo(6, 0)
    t = s.trave(a, b, EI=1e4, EA=1e6)
    s.cerniera(a)
    s.carrello(b)
    s.carico_distribuito(t, -10)          # 10 kN/m verso il basso

    r = s.risolvi()
    r.stampa()
    r.sollecitazioni(t)["M"].max()        # 45 kN·m
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from .assembler import (
    CarichiLocali,
    _lunghezza_theta,
    condensa,
    deformata_locale,
    forze_equivalenti,
    gdl_svincolati,
    griglia,
    k_locale,
    matrice_rotazione,
    recupera_svincolati,
    sollecitazioni,
)
from .core import (
    TIPI_CARICO_LEGACY_NODO,
    TIPI_CARICO_LEGACY_TRAVE,
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


class StrutturaLabile(np.linalg.LinAlgError):
    """La struttura ammette un cinematismo. ``meccanismo`` è il vettore (ndof,)
    degli spostamenti del cinematismo (normalizzato), utile per disegnarlo."""

    def __init__(self, msg: str, meccanismo: Optional[np.ndarray] = None,
                 grado_labilita: int = 0):
        super().__init__(msg)
        self.meccanismo = meccanismo
        self.grado_labilita = grado_labilita


# ---------------------------------------------------------------------------
# Risultato
# ---------------------------------------------------------------------------

@dataclass
class Risultato:
    U: np.ndarray            # spostamenti nodali globali [m, m, rad] (3 per nodo)
    R: np.ndarray            # reazioni vincolari globali (3 per nodo)
    K: np.ndarray            # rigidezza globale assemblata (senza vincoli)
    F: np.ndarray            # carichi nodali equivalenti
    nodi: List[Nodo]
    vincoli: List[Vincolo]
    travi: List[Trave]
    free: List[int]
    fixed: List[int]
    grado_iperstaticita: int = 0
    grado_labilita: int = 0
    _struttura: "Struttura" = field(default=None, repr=False)
    _cache: Dict[int, dict] = field(default_factory=dict, repr=False)

    # -- nodi -----------------------------------------------------------------

    def spostamento(self, nodo: int):
        """(ux, uy, φ) del nodo."""
        return tuple(float(v) for v in self.U[3*nodo:3*nodo+3])

    def reazioni(self) -> Dict[int, dict]:
        """Reazioni per nodo vincolato: componenti globali e negli assi del vincolo."""
        out = {}
        for v in self.vincoli:
            Rx, Ry, M = (float(c) for c in self.R[3*v.nodo:3*v.nodo+3])
            (cx, sx), (cy, sy) = v.assi
            out[v.nodo] = {
                "Rx": Rx, "Ry": Ry, "M": M,
                "R_parallela": Rx * cx + Ry * sx,
                "R_normale": Rx * cy + Ry * sy,
                "tipo": v.tipo,
            }
        return out

    # -- elementi -------------------------------------------------------------

    def _elem(self, t: int) -> dict:
        if t not in self._cache:
            self._cache[t] = self._struttura._risultati_elemento(t, self.U)
        return self._cache[t]

    def forze_estremita(self, t: int) -> np.ndarray:
        """Forze dei nodi sull'elemento in coordinate locali [X_i, Y_i, M_i, X_j, Y_j, M_j]."""
        return self._elem(t)["f_end"].copy()

    def sollecitazioni(self, t: int, n: int = 41) -> dict:
        """Diagrammi N, V, M lungo la trave ``t`` (x dal nodo i)."""
        e = self._elem(t)
        x, xe = griglia(e["L"], e["cl"], n)
        N, V, M = sollecitazioni(e["L"], e["f_end"], e["cl"], xe)
        return {"x": x, "N": N, "V": V, "M": M}

    def deformata(self, t: int, scala: float = 1.0, n: int = 61) -> dict:
        """Deformata della trave in coordinate globali (spostamenti amplificati di ``scala``)."""
        e = self._elem(t)
        x, u, v = deformata_locale(e["L"], e["EI"], e["EA"], e["u_loc"], e["f_end"], e["cl"], n)
        c, s = math.cos(e["theta"]), math.sin(e["theta"])
        ni = self.nodi[self.travi[t].nodo_i]
        X = ni.x + c * x + scala * (c * u - s * v)
        Y = ni.y + s * x + scala * (s * u + c * v)
        return {"x": x, "X": X, "Y": Y, "u": u, "v": v}

    def estremi(self, t: int) -> dict:
        """Valori massimi e minimi di N, V, M sulla trave con la posizione."""
        d = self.sollecitazioni(t, n=201)
        out = {}
        for k in ("N", "V", "M"):
            a = d[k]
            out[k] = {"max": float(a.max()), "x_max": float(d["x"][a.argmax()]),
                      "min": float(a.min()), "x_min": float(d["x"][a.argmin()])}
        return out

    # -- verifica -------------------------------------------------------------

    def verifica_equilibrio(self, tol: float = 1e-6) -> dict:
        """Equilibrio globale: carichi + reazioni = 0 (forze e momento rispetto all'origine)."""
        tot = self.F + self.R
        sfx = float(tot[0::3].sum())
        sfy = float(tot[1::3].sum())
        sm = float(sum(n.x * tot[3*k+1] - n.y * tot[3*k] + tot[3*k+2]
                       for k, n in enumerate(self.nodi)))
        scala = max(1.0, float(np.abs(self.F).max(initial=0)), float(np.abs(self.R).max(initial=0)))
        ok = max(abs(sfx), abs(sfy), abs(sm)) < tol * scala
        return {"ΣFx [kN]": sfx, "ΣFy [kN]": sfy, "ΣM [kN·m]": sm, "ok": ok}

    # -- compatibilità v1 -----------------------------------------------------

    def reazioni_vincolari(self) -> dict:
        comp = ["Rx", "Ry", "Mz"]; unita = ["kN", "kN", "kN·m"]
        out = {}
        for v in self.vincoli:
            flags = [v.ux_fisso or v.kx > 0, v.uy_fisso or v.ky > 0, v.phi_fisso or v.kphi > 0]
            r = self.reazioni()[v.nodo]
            if v.angolo:
                vals = [r["R_parallela"], r["R_normale"], r["M"]]
                nomi = ["Rx'", "Ry'", "Mz"]
            else:
                vals = [r["Rx"], r["Ry"], r["M"]]
                nomi = comp
            for k, f in enumerate(flags):
                if f:
                    out[f"{nomi[k]}_N{v.nodo + 1}"] = (vals[k], unita[k])
        return out

    def spostamenti_nodali(self) -> dict:
        comp = ["ux", "uy", "φ"]; unita = ["m", "m", "rad"]
        return {f"{comp[k]}_N{i + 1}": (float(self.U[3*i+k]), unita[k])
                for i in range(len(self.nodi)) for k in range(3)}

    def stampa(self, decimali: int = 4) -> None:
        print("=" * 60)
        print("  RISULTATI FEM — TRAVI PIANE 2D")
        print("=" * 60)
        print(f"  Nodi: {len(self.nodi)}   Elementi: {len(self.travi)}")
        print(f"  Grado di iperstaticità: {self.grado_iperstaticita}")
        print()
        print("── Spostamenti nodali ──────────────────────────")
        for label, (val, unit) in self.spostamenti_nodali().items():
            print(f"  {label:12s} = {val:+.{decimali}e}  {unit}")
        print()
        print("── Reazioni vincolari ──────────────────────────")
        for label, (val, unit) in self.reazioni_vincolari().items():
            print(f"  {label:12s} = {val:+.{decimali}f}  {unit}")
        print()
        print("── Sollecitazioni (estremi per trave) ──────────")
        for t in range(len(self.travi)):
            e = self.estremi(t)
            print(f"  T{t+1}: M ∈ [{e['M']['min']:+.{decimali}f}, {e['M']['max']:+.{decimali}f}]"
                  f"  V ∈ [{e['V']['min']:+.{decimali}f}, {e['V']['max']:+.{decimali}f}]"
                  f"  N ∈ [{e['N']['min']:+.{decimali}f}, {e['N']['max']:+.{decimali}f}]")
        print()
        eq = self.verifica_equilibrio()
        print("── Verifica equilibrio ─────────────────────────")
        print(f"  ΣFx = {eq['ΣFx [kN]']:+.2e}  ΣFy = {eq['ΣFy [kN]']:+.2e}  ΣM = {eq['ΣM [kN·m]']:+.2e}")
        print(f"  {'OK ✓' if eq['ok'] else 'ERRORE ✗'}")
        print("=" * 60)

    def matrice_K_str(self, decimali: int = 2) -> str:
        n = self.K.shape[0]
        comp = ["ux", "uy", "φ"]
        labels = [f"{comp[i % 3]}{i // 3 + 1}" for i in range(n)]
        lines = [f"{'':8s}" + "".join(f"{l:>12s}" for l in labels)]
        for i, row in enumerate(self.K):
            lines.append(f"{labels[i]:8s}" + "".join(f"{v:+12.{decimali}e}" for v in row))
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Struttura
# ---------------------------------------------------------------------------

class Struttura:
    """Contenitore del modello: nodi, travi, vincoli, carichi."""

    def __init__(self):
        self.nodi: List[Nodo] = []
        self.travi: List[Trave] = []
        self.vincoli: List[Vincolo] = []
        self.carichi: list = []

    # -- costruzione rapida ---------------------------------------------------

    def nodo(self, x: float, y: float, cerniera: bool = False) -> int:
        """Aggiunge un nodo (o restituisce quello già presente nelle stesse coordinate)."""
        for n in self.nodi:
            if abs(n.x - x) < 1e-9 and abs(n.y - y) < 1e-9:
                n.cerniera = n.cerniera or cerniera
                return n.id
        return self.aggiungi_nodo(Nodo(x, y, cerniera))

    def trave(self, i: int, j: int, EI: float = 1e4, EA: float = 1e6,
              cerniera_i: bool = False, cerniera_j: bool = False) -> int:
        return self.aggiungi_trave(Trave(i, j, EI=EI, EA=EA,
                                         cerniera_i=cerniera_i, cerniera_j=cerniera_j))

    def asta(self, i: int, j: int, EA: float = 1e6) -> int:
        """Elemento reticolare (solo sforzo normale)."""
        return self.aggiungi_trave(Trave(i, j, EI=1.0, EA=EA, asta=True))

    def vincolo(self, nodo: int, tipo: str, angolo: float = 0.0, **kw) -> None:
        self.aggiungi_vincolo(Vincolo.da_tipo(nodo, tipo, angolo, **kw))

    def incastro(self, nodo: int, **kw):       self.vincolo(nodo, "incastro", **kw)
    def cerniera(self, nodo: int, **kw):       self.vincolo(nodo, "cerniera", **kw)

    def carrello(self, nodo: int, angolo: float = 0.0, **kw):
        """Carrello (= pendolo): impedisce lo spostamento normale al piano di
        scorrimento, inclinato di ``angolo`` gradi."""
        self.vincolo(nodo, "carrello", angolo, **kw)

    def doppio_pendolo(self, nodo: int, angolo: float = 0.0, **kw):
        """Doppio pendolo (bipendolo): scorre lungo la direzione ``angolo``,
        impedisce la traslazione normale e la rotazione."""
        self.vincolo(nodo, "doppio_pendolo", angolo, **kw)

    def molla(self, nodo: int, kx: float = 0.0, ky: float = 0.0, kphi: float = 0.0):
        self.vincolo(nodo, "molla", kx=kx, ky=ky, kphi=kphi)

    def cerniera_interna(self, nodo: int):
        self.nodi[nodo].cerniera = True

    def forza(self, nodo: int, Fx: float = 0.0, Fy: float = 0.0):
        self.aggiungi_carico(ForzaNodale(nodo, Fx, Fy))

    def coppia(self, nodo: int, M: float):
        self.aggiungi_carico(CoppiaNodale(nodo, M))

    def carico_distribuito(self, trave: int, qi: float, qj: Optional[float] = None,
                           direzione: str = "y", proiezione: bool = False):
        self.aggiungi_carico(CaricoDistribuito(trave, qi, qj, direzione, proiezione))

    def forza_in_campata(self, trave: int, a: float, Fx: float = 0.0, Fy: float = 0.0,
                         locale: bool = False):
        self.aggiungi_carico(ForzaInCampata(trave, a, Fx, Fy, locale))

    def coppia_in_campata(self, trave: int, a: float, M: float):
        self.aggiungi_carico(CoppiaInCampata(trave, a, M))

    def carico_termico(self, trave: int, alpha: float = 1.2e-5, dT: float = 0.0,
                       dT_farfalla: float = 0.0, h: float = 1.0):
        self.aggiungi_carico(CaricoTermico(trave, alpha, dT, dT_farfalla, h))

    def cedimento(self, nodo: int, x: float = 0.0, y: float = 0.0, phi: float = 0.0):
        """Spostamenti imposti sul vincolo del nodo (negli assi del vincolo)."""
        for v in self.vincoli:
            if v.nodo == nodo:
                v.cedimento_x, v.cedimento_y, v.cedimento_phi = x, y, phi
                return
        raise ValueError(f"Il nodo {nodo + 1} non ha vincoli: impossibile imporre un cedimento.")

    # -- API v1 ---------------------------------------------------------------

    def aggiungi_nodo(self, nodo: Nodo) -> int:
        nodo.id = len(self.nodi)
        self.nodi.append(nodo)
        return nodo.id

    def aggiungi_trave(self, trave: Trave) -> int:
        self._valida_trave(trave)
        trave.id = len(self.travi)
        self.travi.append(trave)
        return trave.id

    def aggiungi_vincolo(self, vincolo: Vincolo) -> None:
        self._valida_nodo_idx(vincolo.nodo, "vincolo")
        self.vincoli = [v for v in self.vincoli if v.nodo != vincolo.nodo]
        self.vincoli.append(vincolo)

    def aggiungi_carico(self, carico) -> None:
        if isinstance(carico, Carico):
            self._valida_carico_legacy(carico)
            carico = carico.converti()
        self._valida_carico(carico)
        self.carichi.append(carico)

    # -- validazione ----------------------------------------------------------

    def _valida_nodo_idx(self, idx: int, contesto: str = "") -> None:
        if not (0 <= idx < len(self.nodi)):
            raise IndexError(f"Indice nodo {idx} non valido ({contesto}). Nodi presenti: {len(self.nodi)}")

    def _valida_trave_idx(self, idx: int, contesto: str = "") -> None:
        if not (0 <= idx < len(self.travi)):
            raise IndexError(f"Indice trave {idx} non valido ({contesto}). Travi presenti: {len(self.travi)}")

    def _valida_trave(self, t: Trave) -> None:
        self._valida_nodo_idx(t.nodo_i, "trave nodo_i")
        self._valida_nodo_idx(t.nodo_j, "trave nodo_j")
        if t.nodo_i == t.nodo_j:
            raise ValueError("Trave con nodo_i == nodo_j.")
        if t.EI <= 0:
            raise ValueError(f"EI deve essere > 0, trovato {t.EI}.")
        if t.EA <= 0:
            raise ValueError(f"EA deve essere > 0, trovato {t.EA}.")
        ni, nj = self.nodi[t.nodo_i], self.nodi[t.nodo_j]
        if math.hypot(nj.x - ni.x, nj.y - ni.y) < 1e-12:
            raise ValueError(f"I nodi {t.nodo_i + 1} e {t.nodo_j + 1} coincidono.")

    def _valida_carico_legacy(self, c: Carico) -> None:
        if c.type in TIPI_CARICO_LEGACY_NODO:
            self._valida_nodo_idx(c.obj, f"carico {c.type}")
        elif c.type in TIPI_CARICO_LEGACY_TRAVE:
            self._valida_trave_idx(c.obj, f"carico {c.type}")
        else:
            raise ValueError(f"Tipo carico sconosciuto: '{c.type}'")

    def _valida_carico(self, c) -> None:
        if isinstance(c, (ForzaNodale, CoppiaNodale)):
            self._valida_nodo_idx(c.nodo, "carico nodale")
            return
        if not isinstance(c, (CaricoDistribuito, ForzaInCampata, CoppiaInCampata, CaricoTermico)):
            raise ValueError(f"Tipo carico sconosciuto: {type(c).__name__}")
        self._valida_trave_idx(c.trave, type(c).__name__)
        el = self.travi[c.trave]
        L, _ = _lunghezza_theta(el, self.nodi)
        if isinstance(c, (ForzaInCampata, CoppiaInCampata)) and not (0 <= c.a <= L + 1e-9):
            raise ValueError(f"Posizione a={c.a} fuori dalla trave {c.trave + 1} (L={L:.3f} m).")
        if el.asta:
            trasv = (
                isinstance(c, CoppiaInCampata)
                or (isinstance(c, CaricoTermico) and c.dT_farfalla)
                or (isinstance(c, CaricoDistribuito) and not
                    (c.direzione == "assiale"))
                or isinstance(c, ForzaInCampata)
            )
            if trasv:
                raise ValueError(
                    f"L'elemento {c.trave + 1} è un'asta: accetta solo carichi nodali, "
                    "assiali o termici uniformi. Usa una trave, oppure spezza l'asta con un nodo."
                )

    # -- carichi locali per elemento ------------------------------------------

    def _carichi_locali(self, t: int) -> CarichiLocali:
        el = self.travi[t]
        L, theta = _lunghezza_theta(el, self.nodi)
        c, s = math.cos(theta), math.sin(theta)
        cl = CarichiLocali()
        for q in self.carichi:
            if getattr(q, "trave", None) != t:
                continue
            if isinstance(q, CaricoDistribuito):
                if q.direzione == "perp":
                    cl.distribuiti.append((0.0, 0.0, q.qi, q.qj))
                elif q.direzione == "assiale":
                    cl.distribuiti.append((q.qi, q.qj, 0.0, 0.0))
                else:
                    gx, gy = (1.0, 0.0) if q.direzione == "x" else (0.0, 1.0)
                    f = 1.0
                    if q.proiezione:
                        f = abs(s) if q.direzione == "x" else abs(c)
                    ax, ay = f * (c * gx + s * gy), f * (-s * gx + c * gy)
                    cl.distribuiti.append((ax * q.qi, ax * q.qj, ay * q.qi, ay * q.qj))
            elif isinstance(q, ForzaInCampata):
                if q.locale:
                    cl.forze.append((q.a, q.Fx, q.Fy))
                else:
                    cl.forze.append((q.a, c * q.Fx + s * q.Fy, -s * q.Fx + c * q.Fy))
            elif isinstance(q, CoppiaInCampata):
                cl.coppie.append((q.a, q.M))
            elif isinstance(q, CaricoTermico):
                cl.eps0 += q.alpha * q.dT
                if q.dT_farfalla:
                    cl.kappa0 += q.alpha * q.dT_farfalla / q.h
        return cl

    def _dati_elemento(self, t: int) -> dict:
        el = self.travi[t]
        L, theta = _lunghezza_theta(el, self.nodi)
        cl = self._carichi_locali(t)
        K = k_locale(L, el.EI, el.EA)
        f = forze_equivalenti(L, el.EI, el.EA, cl)
        rel = gdl_svincolati(el, self.nodi)
        Kc, fc = condensa(K, f, rel)
        T = matrice_rotazione(theta)
        i, j = el.nodo_i, el.nodo_j
        return dict(L=L, theta=theta, cl=cl, K=K, f=f, rel=rel, Kc=Kc, fc=fc, T=T,
                    EI=el.EI, EA=el.EA, dofs=[3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2])

    def _risultati_elemento(self, t: int, U: np.ndarray) -> dict:
        d = self._dati_elemento(t)
        u_loc = d["T"] @ U[d["dofs"]]
        u_loc = recupera_svincolati(d["K"], d["f"], u_loc, d["rel"])
        f_end = d["K"] @ u_loc - d["f"]
        f_end[d["rel"]] = 0.0
        d.update(u_loc=u_loc, f_end=f_end)
        return d

    # -- assemblaggio ---------------------------------------------------------

    def _assembla(self):
        ndof = 3 * len(self.nodi)
        K = np.zeros((ndof, ndof))
        F = np.zeros(ndof)
        for t in range(len(self.travi)):
            d = self._dati_elemento(t)
            Kg = d["T"].T @ d["Kc"] @ d["T"]
            fg = d["T"].T @ d["fc"]
            K[np.ix_(d["dofs"], d["dofs"])] += Kg
            F[d["dofs"]] += fg
        for q in self.carichi:
            if isinstance(q, ForzaNodale):
                F[3*q.nodo] += q.Fx
                F[3*q.nodo+1] += q.Fy
            elif isinstance(q, CoppiaNodale):
                F[3*q.nodo+2] += q.M
        return K, F

    def _vincoli_matrici(self, ndof: int):
        """Righe di vincolo C·u = d e matrice delle molle."""
        C, d = [], []
        Ks = np.zeros((ndof, ndof))
        for v in self.vincoli:
            n = v.nodo
            ex, ey = v.assi
            dirs = [(ex, v.ux_fisso, v.kx, v.cedimento_x),
                    (ey, v.uy_fisso, v.ky, v.cedimento_y)]
            for (a, b), fisso, k, ced in dirs:
                vec = np.zeros(ndof); vec[3*n] = a; vec[3*n+1] = b
                if fisso:
                    C.append(vec); d.append(ced)
                elif k > 0:
                    Ks += k * np.outer(vec, vec)
            if v.phi_fisso:
                vec = np.zeros(ndof); vec[3*n+2] = 1.0
                C.append(vec); d.append(v.cedimento_phi)
            elif v.kphi > 0:
                Ks[3*n+2, 3*n+2] += v.kphi
        C = np.array(C) if C else np.zeros((0, ndof))
        return C, np.array(d, dtype=float), Ks

    # -- analisi cinematica / statica -----------------------------------------

    def gradi(self, C=None, Ks=None, attivi=None):
        """(grado di iperstaticità, grado di labilità) dal rango della matrice di equilibrio."""
        ndof = 3 * len(self.nodi)
        if C is None:
            K, _ = self._assembla()
            C, _, Ks = self._vincoli_matrici(ndof)
            attivi = self._gdl_attivi(K + Ks, C)
        cols = []
        for t, el in enumerate(self.travi):
            L, theta = _lunghezza_theta(el, self.nodi)
            T = matrice_rotazione(theta)
            rel = gdl_svincolati(el, self.nodi)
            dofs = [3*el.nodo_i, 3*el.nodo_i+1, 3*el.nodo_i+2,
                    3*el.nodo_j, 3*el.nodo_j+1, 3*el.nodo_j+2]
            base = [np.array([-1, 0, 0, 1, 0, 0.0])]
            if 2 not in rel:
                base.append(np.array([0, 1/L, 1, 0, -1/L, 0.0]))
            if 5 not in rel:
                base.append(np.array([0, 1/L, 0, 0, -1/L, 1.0]))
            for b in base:
                col = np.zeros(ndof)
                col[dofs] = T.T @ b
                cols.append(col)
        for row in C:
            cols.append(row)
        for k in range(ndof):
            if Ks is not None and np.any(Ks[k] != 0):
                col = np.zeros(ndof); col[k] = 1.0
                cols.append(col)
        if not cols:
            return 0, int(len(attivi))
        A = np.array(cols).T[attivi]
        rango = np.linalg.matrix_rank(A, tol=1e-9 * max(1.0, np.abs(A).max()))
        return int(A.shape[1] - rango), int(A.shape[0] - rango)

    @staticmethod
    def _gdl_attivi(Ktot, C):
        """Esclude i GDL senza alcuna rigidezza né vincolo (rotazioni dei nodi
        di sole aste o di cerniere interne)."""
        ndof = Ktot.shape[0]
        scala = max(1.0, np.abs(Ktot).max(initial=0))
        att = []
        for k in range(ndof):
            if np.abs(Ktot[k]).max() > 1e-14 * scala or (C.shape[0] and np.abs(C[:, k]).max() > 0):
                att.append(k)
        return att

    # -- soluzione ------------------------------------------------------------

    def risolvi(self) -> Risultato:
        """Assembla e risolve. Solleva ``StrutturaLabile`` se la struttura è labile."""
        if len(self.nodi) < 2:
            raise ValueError("Servono almeno 2 nodi.")
        if len(self.travi) < 1:
            raise ValueError("Serve almeno 1 trave.")
        if len(self.vincoli) < 1:
            raise ValueError("Serve almeno 1 vincolo.")

        ndof = 3 * len(self.nodi)
        K, F = self._assembla()
        C, d, Ks = self._vincoli_matrici(ndof)
        Kt = K + Ks
        att = self._gdl_attivi(Kt, C)
        inattivi = sorted(set(range(ndof)) - set(att))

        for k in inattivi:
            if abs(F[k]) > 1e-12:
                n = k // 3 + 1
                raise StrutturaLabile(
                    f"Carico applicato a un GDL senza rigidezza (nodo {n}): "
                    "per esempio una coppia su una cerniera o su un nodo di sole aste.",
                    grado_labilita=1)

        i_grado, l_grado = self.gradi(C, Ks, att)

        Ka = Kt[np.ix_(att, att)]
        Fa = F[att]
        Ca = C[:, att]

        # soluzione particolare dei vincoli e base del nucleo
        if Ca.shape[0]:
            Uc, sc, Vt = np.linalg.svd(Ca)
            rc = int((sc > 1e-10 * max(1.0, sc.max(initial=0))).sum())
            up = np.linalg.pinv(Ca) @ d
            if np.abs(Ca @ up - d).max(initial=0) > 1e-9 * max(1.0, np.abs(d).max(initial=0)):
                raise ValueError("Cedimenti vincolari incompatibili fra loro.")
            Z = Vt[rc:].T
        else:
            up = np.zeros(len(att))
            Z = np.eye(len(att))

        Kz = Z.T @ Ka @ Z
        Fz = Z.T @ (Fa - Ka @ up)
        if Kz.size:
            w, V = np.linalg.eigh(Kz)
            wmax = max(abs(w).max(), 1e-300)
            piccoli = w < 1e-10 * wmax
            if piccoli.any():
                mecc = np.zeros(ndof)
                mecc[att] = Z @ V[:, 0]
                mecc /= max(np.abs(mecc).max(), 1e-300)
                nm = int(piccoli.sum())
                raise StrutturaLabile(
                    f"Struttura labile ({nm} {'cinematismo' if nm == 1 else 'cinematismi indipendenti'}): "
                    "controlla vincoli e cerniere.",
                    meccanismo=mecc, grado_labilita=int(piccoli.sum()))
            y = V @ ((V.T @ Fz) / w)
            ua = up + Z @ y
        else:
            ua = up

        U = np.zeros(ndof)
        U[att] = ua
        R = K @ U - F                      # forze dei vincoli (e delle molle) sulla struttura
        R[inattivi] = 0.0
        for k in range(ndof):
            if not (C.shape[0] and np.abs(C[:, k]).max() > 0) and not np.any(Ks[k]):
                R[k] = 0.0                 # GDL liberi: rumore numerico

        fixed = sorted({int(np.argmax(np.abs(r))) for r in C if np.count_nonzero(r) == 1})
        free = [k for k in att if k not in fixed]
        return Risultato(U=U, R=R, K=K, F=F, nodi=self.nodi, vincoli=self.vincoli,
                         travi=self.travi, free=free, fixed=fixed,
                         grado_iperstaticita=i_grado, grado_labilita=l_grado,
                         _struttura=self)

    # -- info -----------------------------------------------------------------

    def info(self) -> str:
        return "\n".join([
            "Struttura FEM 2D",
            f"  Nodi:    {len(self.nodi)}",
            f"  Travi:   {len(self.travi)}",
            f"  Vincoli: {len(self.vincoli)} nodi vincolati",
            f"  Carichi: {len(self.carichi)}",
            f"  GDL totali: {3 * len(self.nodi)}",
        ])
