"""
fem_travi/io.py
===============
Modello in formato JSON (lo stesso usato dall'interfaccia) ed esportazione CSV.

Formato del modello
-------------------
{
  "nodi":     [{"x": 0, "y": 0, "cerniera": false}, ...],
  "elementi": [{"i": 0, "j": 1, "tipo": "trave"|"asta", "EI": 1e4, "EA": 1e6,
                "inestensibile": false, "cerniera_i": false, "cerniera_j": false}, ...],
  "vincoli":  [{"nodo": 0, "tipo": "incastro"|"cerniera"|"carrello"|
                "doppio_pendolo"|"bloccarotazione"|"molla"|"personalizzato",
                "angolo": 0, "ux": true, "uy": true, "phi": true,
                "kx": 0, "ky": 0, "kphi": 0,
                "cedimento": {"x": 0, "y": 0, "phi": 0}}, ...],
  "carichi":  [{"tipo": "forza", "nodo": 1, "Fx": 0, "Fy": -10},
               {"tipo": "coppia", "nodo": 1, "M": 5},
               {"tipo": "distribuito", "elem": 0, "qi": -10, "qj": -10,
                "direzione": "y", "proiezione": false},
               {"tipo": "forza_campata", "elem": 0, "a": 2, "Fx": 0, "Fy": -20},
               {"tipo": "coppia_campata", "elem": 0, "a": 2, "M": 10},
               {"tipo": "termico", "elem": 0, "alpha": 1.2e-5, "dT": 20,
                "dT_farfalla": 0, "h": 0.5}]
}
"""

from __future__ import annotations

import csv
import io as _io
import json
import math
from typing import Optional

import numpy as np

from .core import TIPI_VINCOLO, Nodo, Trave, Vincolo
from .solver import Risultato, Struttura, StrutturaLabile


def _f(v, default=0.0) -> float:
    try:
        out = float(v)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def struttura_da_dict(m: dict) -> Struttura:
    s = Struttura()
    for n in m.get("nodi", []):
        s.aggiungi_nodo(Nodo(_f(n["x"]), _f(n["y"]), bool(n.get("cerniera", False))))
    for e in m.get("elementi", []):
        asta = e.get("tipo", "trave") == "asta"
        EI = _f(e.get("EI"), 1e4)
        EA = _f(e.get("EA"), 1e6)
        if e.get("inestensibile") and not asta:
            EA = EI * 1e6          # assialmente rigida (rapporto sufficiente, ben condizionato)
        s.aggiungi_trave(Trave(int(e["i"]), int(e["j"]),
                               EI=EI if not asta else 1.0,
                               EA=EA,
                               cerniera_i=bool(e.get("cerniera_i")),
                               cerniera_j=bool(e.get("cerniera_j")),
                               asta=asta))
    for v in m.get("vincoli", []):
        tipo = v.get("tipo", "incastro")
        ced = v.get("cedimento") or {}
        kw = dict(angolo=_f(v.get("angolo")),
                  kx=_f(v.get("kx")), ky=_f(v.get("ky")), kphi=_f(v.get("kphi")),
                  cedimento_x=_f(ced.get("x")), cedimento_y=_f(ced.get("y")),
                  cedimento_phi=_f(ced.get("phi")), tipo=tipo)
        if tipo in TIPI_VINCOLO:
            ux, uy, phi = TIPI_VINCOLO[tipo]
        else:
            ux, uy, phi = bool(v.get("ux")), bool(v.get("uy")), bool(v.get("phi"))
        if tipo == "molla":
            ux = uy = phi = False
        s.aggiungi_vincolo(Vincolo(int(v["nodo"]), ux, uy, phi, **kw))
    for c in m.get("carichi", []):
        t = c.get("tipo")
        if t == "forza":
            s.forza(int(c["nodo"]), _f(c.get("Fx")), _f(c.get("Fy")))
        elif t == "coppia":
            s.coppia(int(c["nodo"]), _f(c.get("M")))
        elif t == "distribuito":
            qi = _f(c.get("qi"))
            qj = _f(c.get("qj"), qi) if c.get("qj") is not None else qi
            s.carico_distribuito(int(c["elem"]), qi, qj, c.get("direzione", "y"),
                                 bool(c.get("proiezione")))
        elif t == "forza_campata":
            s.forza_in_campata(int(c["elem"]), _f(c.get("a")), _f(c.get("Fx")), _f(c.get("Fy")),
                               bool(c.get("locale")))
        elif t == "coppia_campata":
            s.coppia_in_campata(int(c["elem"]), _f(c.get("a")), _f(c.get("M")))
        elif t == "termico":
            s.carico_termico(int(c["elem"]), _f(c.get("alpha"), 1.2e-5), _f(c.get("dT")),
                             _f(c.get("dT_farfalla")), _f(c.get("h"), 1.0))
        else:
            raise ValueError(f"Tipo carico sconosciuto: '{t}'")
    return s


def carica_json(percorso: str) -> Struttura:
    with open(percorso, encoding="utf-8") as fh:
        return struttura_da_dict(json.load(fh))


def _l(a) -> list:
    return [round(float(v), 10) for v in a]


def risultato_a_dict(r: Risultato, n_punti: int = 41, scala: float = 1.0) -> dict:
    """Risultati in forma serializzabile (usata dall'interfaccia)."""
    elementi = []
    for t, el in enumerate(r.travi):
        d = r.sollecitazioni(t, n_punti)
        df = r.deformata(t, scala=scala, n=n_punti)
        e = r._elem(t)
        elementi.append({
            "x": _l(d["x"]), "N": _l(d["N"]), "V": _l(d["V"]), "M": _l(d["M"]),
            "def_x": _l(df["x"]), "def_u": _l(df["u"]), "def_v": _l(df["v"]),
            "L": e["L"], "theta": e["theta"],
            "f_end": _l(e["f_end"]),
            "estremi": r.estremi(t),
        })
    reaz = [{"nodo": n, **{k: v for k, v in rv.items()}} for n, rv in r.reazioni().items()]
    eq = r.verifica_equilibrio()
    return {
        "ok": True,
        "U": _l(r.U),
        "reazioni": reaz,
        "elementi": elementi,
        "grado_iperstaticita": r.grado_iperstaticita,
        "grado_labilita": r.grado_labilita,
        "equilibrio": {"Fx": eq["ΣFx [kN]"], "Fy": eq["ΣFy [kN]"], "M": eq["ΣM [kN·m]"], "ok": bool(eq["ok"])},
    }


def risolvi_dict(m: dict, n_punti: int = 41) -> dict:
    """Risolve un modello JSON e restituisce un dict; gli errori diventano messaggi."""
    try:
        s = struttura_da_dict(m)
        r = s.risolvi()
        return risultato_a_dict(r, n_punti)
    except StrutturaLabile as e:
        out = {"ok": False, "errore": str(e), "labile": True,
               "grado_labilita": e.grado_labilita}
        if e.meccanismo is not None:
            out["meccanismo"] = _l(e.meccanismo)
        try:
            out["grado_iperstaticita"], out["grado_labilita"] = s.gradi()
        except Exception:
            pass
        return out
    except (ValueError, IndexError, KeyError, np.linalg.LinAlgError) as e:
        return {"ok": False, "errore": str(e)}


def esporta_csv(r: Risultato, percorso: Optional[str] = None, sep: str = ";",
                decimale: str = ",", n_punti: int = 21) -> str:
    """CSV con spostamenti, reazioni e sollecitazioni. Di default usa ';' e la
    virgola decimale (formato che Excel in italiano apre direttamente)."""
    def num(v):
        s = f"{v:.6g}"
        return s.replace(".", decimale) if decimale != "." else s

    buf = _io.StringIO()
    w = csv.writer(buf, delimiter=sep, lineterminator="\n")
    w.writerow(["SPOSTAMENTI NODALI"])
    w.writerow(["nodo", "x [m]", "y [m]", "ux [m]", "uy [m]", "phi [rad]"])
    for k, n in enumerate(r.nodi):
        ux, uy, ph = r.spostamento(k)
        w.writerow([k + 1, num(n.x), num(n.y), num(ux), num(uy), num(ph)])
    w.writerow([])
    w.writerow(["REAZIONI VINCOLARI"])
    w.writerow(["nodo", "vincolo", "Rx [kN]", "Ry [kN]", "M [kN m]"])
    for n, rv in r.reazioni().items():
        w.writerow([n + 1, rv["tipo"], num(rv["Rx"]), num(rv["Ry"]), num(rv["M"])])
    w.writerow([])
    w.writerow(["SOLLECITAZIONI"])
    w.writerow(["trave", "x [m]", "N [kN]", "V [kN]", "M [kN m]"])
    for t in range(len(r.travi)):
        d = r.sollecitazioni(t, n_punti)
        for x, N, V, M in zip(d["x"], d["N"], d["V"], d["M"]):
            w.writerow([t + 1, num(x), num(N), num(V), num(M)])
    testo = buf.getvalue()
    if percorso:
        with open(percorso, "w", encoding="utf-8-sig", newline="") as fh:
            fh.write(testo)
    return testo
