"""
fem_travi/core.py
=================
Modello dati del solver: nodi, elementi, vincoli e carichi.

Convenzione di segno (unica in tutta la libreria)
-------------------------------------------------
Assi globali: x verso destra, y verso l'alto, rotazioni antiorarie positive.

- Forze e carichi: componenti positive nel verso degli assi globali.
  Un carico di 10 kN verso il basso si scrive  Fy = -10.
- Coppie: positive se antiorarie.
- Carichi distribuiti:  direzione "y"/"x" = componenti globali (per metro
  di trave); direzione "perp" = normale all'asse della trave, positiva
  verso l'asse y locale (per una trave da sinistra a destra: verso l'alto).
- Sollecitazioni: N > 0 trazione; M > 0 se tende le fibre inferiori
  (lato -y locale); V = dM/dx.

Unità consigliate: kN, m (EI in kN·m², EA in kN).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Geometria
# ---------------------------------------------------------------------------

@dataclass
class Nodo:
    """Nodo con coordinate globali [m].

    ``cerniera=True`` rende il nodo una cerniera interna: tutte le travi
    che vi arrivano sono svincolate a rotazione in quel punto.
    """
    x: float
    y: float
    cerniera: bool = False
    id: int = field(default=-1, repr=False)

    def __post_init__(self):
        self.x = float(self.x)
        self.y = float(self.y)


@dataclass
class Trave:
    """Elemento trave Euler-Bernoulli (o asta, se ``asta=True``).

    nodo_i, nodo_j : indici dei nodi estremi
    EI : rigidezza flessionale [kN·m²]
    EA : rigidezza assiale [kN]
    cerniera_i, cerniera_j : svincolo a rotazione all'estremo (cerniera interna)
    asta : elemento reticolare (solo sforzo normale, entrambi gli estremi
           incernierati, nessun carico trasversale in campata)
    theta : parametro della v1, ignorato (l'angolo si ricava dai nodi)
    """
    nodo_i: int
    nodo_j: int
    EI: float = 1.0e4
    EA: float = 1.0e6
    theta: Optional[float] = None
    cerniera_i: bool = False
    cerniera_j: bool = False
    asta: bool = False
    id: int = field(default=-1, repr=False)


# ---------------------------------------------------------------------------
# Vincoli
# ---------------------------------------------------------------------------

#: GDL bloccati (u, v, φ) negli assi del vincolo per ciascun tipo.
#: Gli assi del vincolo sono quelli globali ruotati di ``angolo``.
TIPI_VINCOLO = {
    "incastro":        (True,  True,  True),
    "cerniera":        (True,  True,  False),
    "carrello":        (False, True,  False),   # scorre lungo x' (= pendolo)
    "doppio_pendolo":  (False, True,  True),    # scorre lungo x', non ruota
    "bloccarotazione": (False, False, True),
    "molla":           (False, False, False),
}


@dataclass
class Vincolo:
    """Vincolo esterno su un nodo.

    ux_fisso, uy_fisso, phi_fisso : GDL bloccati negli assi del vincolo
    angolo : rotazione antioraria [gradi] degli assi del vincolo rispetto
             agli assi globali (es. carrello su piano inclinato di 30° →
             ``angolo=30``)
    kx, ky, kphi : molle elastiche [kN/m, kN/m, kN·m/rad] negli assi del vincolo
    cedimento_x, cedimento_y, cedimento_phi : spostamenti imposti sui GDL
             bloccati [m, m, rad] (cedimenti vincolari), negli assi del vincolo
    tipo : etichetta descrittiva (vedi ``TIPI_VINCOLO``)
    """
    nodo: int
    ux_fisso: bool = True
    uy_fisso: bool = True
    phi_fisso: bool = True
    angolo: float = 0.0
    kx: float = 0.0
    ky: float = 0.0
    kphi: float = 0.0
    cedimento_x: float = 0.0
    cedimento_y: float = 0.0
    cedimento_phi: float = 0.0
    tipo: str = ""

    @classmethod
    def da_tipo(cls, nodo: int, tipo: str, angolo: float = 0.0, **kw) -> "Vincolo":
        if tipo not in TIPI_VINCOLO:
            raise ValueError(
                f"Tipo vincolo sconosciuto: '{tipo}'. "
                f"Tipi validi: {', '.join(TIPI_VINCOLO)}"
            )
        ux, uy, phi = TIPI_VINCOLO[tipo]
        return cls(nodo, ux, uy, phi, angolo=angolo, tipo=tipo, **kw)

    @property
    def assi(self):
        """Versori (x', y') degli assi del vincolo in coordinate globali."""
        a = math.radians(self.angolo)
        c, s = math.cos(a), math.sin(a)
        return (c, s), (-s, c)


# ---------------------------------------------------------------------------
# Carichi (API nuova)
# ---------------------------------------------------------------------------

@dataclass
class ForzaNodale:
    """Forza concentrata su un nodo, componenti globali [kN]."""
    nodo: int
    Fx: float = 0.0
    Fy: float = 0.0


@dataclass
class CoppiaNodale:
    """Coppia concentrata su un nodo [kN·m], antioraria positiva."""
    nodo: int
    M: float


@dataclass
class CaricoDistribuito:
    """Carico distribuito lineare (uniforme, triangolare, trapezio) su tutta la trave.

    qi, qj : intensità all'estremo i e j [kN/m] (qj=None → uniforme)
    direzione : "y" | "x" (globali) oppure "perp" | "assiale" (locali)
    proiezione : solo per "y"/"x": intensità riferita alla proiezione
                 orizzontale/verticale della trave (es. neve su falda)
    """
    trave: int
    qi: float
    qj: Optional[float] = None
    direzione: str = "y"
    proiezione: bool = False

    def __post_init__(self):
        if self.qj is None:
            self.qj = self.qi
        if self.direzione not in ("x", "y", "perp", "assiale"):
            raise ValueError(f"Direzione carico non valida: '{self.direzione}'")


@dataclass
class ForzaInCampata:
    """Forza concentrata dentro la trave, a distanza ``a`` [m] dal nodo i.

    Fx, Fy globali (o locali assiale/trasversale se ``locale=True``).
    """
    trave: int
    a: float
    Fx: float = 0.0
    Fy: float = 0.0
    locale: bool = False


@dataclass
class CoppiaInCampata:
    """Coppia concentrata dentro la trave a distanza ``a`` [m] dal nodo i."""
    trave: int
    a: float
    M: float


@dataclass
class CaricoTermico:
    """Variazione termica sulla trave.

    alpha : coefficiente di dilatazione [1/°C]
    dT : variazione uniforme [°C] (allungamento se > 0)
    dT_farfalla : T_inferiore − T_superiore [°C] (lato inferiore = −y locale)
    h : altezza della sezione [m] (serve solo se dT_farfalla ≠ 0)
    """
    trave: int
    alpha: float = 1.2e-5
    dT: float = 0.0
    dT_farfalla: float = 0.0
    h: float = 1.0


# ---------------------------------------------------------------------------
# Carico della v1 (mantenuto per compatibilità)
# ---------------------------------------------------------------------------

TIPI_CARICO_LEGACY_NODO = {"nodo_fx", "nodo_fy", "nodo_m"}
TIPI_CARICO_LEGACY_TRAVE = {"uniforme", "triangolare_sx", "triangolare_dx"}


@dataclass
class Carico:
    """Carico con la convenzione della v1 (deprecato, usare le classi sopra).

    'nodo_fy' val>0 verso il basso · 'nodo_fx' val>0 verso sinistra ·
    'nodo_m' val>0 orario · distribuiti val>0 verso il basso (normale all'asse).
    """
    type: str
    obj: int
    val1: float
    val2: float = 0.0

    def converti(self):
        """Traduce nella nuova convenzione."""
        t, o, v = self.type, self.obj, self.val1
        if t == "nodo_fy":
            return ForzaNodale(o, Fy=-v)
        if t == "nodo_fx":
            return ForzaNodale(o, Fx=-v)
        if t == "nodo_m":
            return CoppiaNodale(o, M=-v)
        if t == "uniforme":
            return CaricoDistribuito(o, -v, -v, direzione="perp")
        if t == "triangolare_sx":
            return CaricoDistribuito(o, -v, 0.0, direzione="perp")
        if t == "triangolare_dx":
            return CaricoDistribuito(o, 0.0, -v, direzione="perp")
        raise ValueError(f"Tipo carico sconosciuto: '{t}'")
