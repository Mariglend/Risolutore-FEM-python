"""
Test del solver v2 contro soluzioni analitiche note (tabelle di Scienza delle
Costruzioni). Convenzione: y verso l'alto, M > 0 tende le fibre inferiori.

    pytest -v
"""

import json
import math

import numpy as np
import pytest

from fem_travi import Struttura, StrutturaLabile, esporta_csv, risolvi_dict

EI = 1.0e4


def appoggiata(L=6.0, **kw):
    s = Struttura()
    a, b = s.nodo(0, 0), s.nodo(L, 0)
    t = s.trave(a, b, EI=kw.get("EI", EI), EA=kw.get("EA", 1e6))
    s.cerniera(a)
    s.carrello(b)
    return s, t


# ---------------------------------------------------------------------------
# Diagrammi esatti sotto carichi in campata (bug della v1)
# ---------------------------------------------------------------------------

class TestDiagrammi:

    def test_parabola_momento_appoggiata(self):
        L, q = 6.0, 10.0
        s, t = appoggiata(L)
        s.carico_distribuito(t, -q)
        r = s.risolvi()
        d = r.sollecitazioni(t, n=61)
        x = d["x"]
        assert np.allclose(d["M"], q * L * x / 2 - q * x**2 / 2, atol=1e-8)
        assert np.allclose(d["V"], q * L / 2 - q * x, atol=1e-8)
        assert r.estremi(t)["M"]["max"] == pytest.approx(q * L**2 / 8)

    def test_incastro_incastro_uniforme(self):
        L, q = 6.0, 10.0
        s = Struttura()
        a, b = s.nodo(0, 0), s.nodo(L, 0)
        t = s.trave(a, b, EI=EI)
        s.incastro(a); s.incastro(b)
        s.carico_distribuito(t, -q)
        r = s.risolvi()
        d = r.sollecitazioni(t, n=41)
        assert d["M"][0] == pytest.approx(-q * L**2 / 12)
        assert d["M"][20] == pytest.approx(q * L**2 / 24)
        assert d["M"][-1] == pytest.approx(-q * L**2 / 12)
        assert r.grado_iperstaticita == 3

    def test_freccia_deformata_in_campata(self):
        L, q = 6.0, 10.0
        s, t = appoggiata(L)
        s.carico_distribuito(t, -q)
        df = s.risolvi().deformata(t, n=61)
        assert df["v"][30] == pytest.approx(-5 * q * L**4 / (384 * EI), rel=1e-5)

    def test_forza_in_campata(self):
        L, P, a = 6.0, 12.0, 2.0
        s, t = appoggiata(L)
        s.forza_in_campata(t, a, Fy=-P)
        r = s.risolvi()
        assert r.reazioni()[0]["Ry"] == pytest.approx(P * (L - a) / L)
        assert r.estremi(t)["M"]["max"] == pytest.approx(P * a * (L - a) / L)
        d = r.sollecitazioni(t)
        k = list(d["x"]).index(a)
        assert d["V"][k] == pytest.approx(P * (L - a) / L)
        assert d["V"][k + 1] == pytest.approx(-P * a / L)

    def test_coppia_in_campata_salto(self):
        L, M0 = 6.0, 12.0
        s, t = appoggiata(L)
        s.coppia_in_campata(t, 3.0, M0)
        d = s.risolvi().sollecitazioni(t)
        k = list(d["x"]).index(3.0)
        assert d["M"][k] - d["M"][k + 1] == pytest.approx(M0, rel=1e-6)

    def test_trapezio(self):
        L = 6.0
        s, t = appoggiata(L)
        s.carico_distribuito(t, -10, -20)
        r = s.risolvi()
        assert r.reazioni()[0]["Ry"] == pytest.approx(40.0)
        assert r.reazioni()[1]["Ry"] == pytest.approx(50.0)

    def test_carico_per_proiezione_su_falda(self):
        s = Struttura()
        a, b = s.nodo(0, 0), s.nodo(4, 3)
        t = s.trave(a, b)
        s.cerniera(a); s.carrello(b)
        s.carico_distribuito(t, -10, direzione="y", proiezione=True)
        r = s.risolvi()
        assert r.reazioni()[0]["Ry"] + r.reazioni()[1]["Ry"] == pytest.approx(40.0)
        assert r.estremi(t)["M"]["max"] == pytest.approx(10 * 4**2 / 8, rel=1e-6)


# ---------------------------------------------------------------------------
# Vincoli
# ---------------------------------------------------------------------------

class TestVincoli:

    def test_capriata_pratt(self):
        s = Struttura()
        for p in [(0, 0), (2, 0), (4, 0), (6, 0), (1, 2), (3, 2), (5, 2)]:
            s.nodo(*p)
        for i, j in [(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (0, 4), (1, 4),
                     (1, 5), (2, 5), (2, 6), (3, 6)]:
            s.asta(i, j, EA=1e5)
        s.cerniera(0); s.carrello(3)
        for k in (4, 5, 6):
            s.forza(k, Fy=-20)
        r = s.risolvi()
        assert r.reazioni()[0]["Ry"] == pytest.approx(30.0)
        assert r.reazioni()[3]["Ry"] == pytest.approx(30.0)
        assert r.grado_iperstaticita == 0
        # Ritter: corrente superiore 4-5 (momento attorno al nodo (2,0))
        assert r.sollecitazioni(3)["N"][0] == pytest.approx(-(30 * 2 - 20 * 1) / 2)
        # corrente inferiore 1-2 (momento attorno al nodo (3,2))
        assert r.sollecitazioni(1)["N"][0] == pytest.approx((30 * 3 - 20 * 2) / 2)

    def test_trave_gerber(self):
        s = Struttura()
        a, b, c = s.nodo(0, 0), s.nodo(3, 0), s.nodo(6, 0)
        s.trave(a, b); t2 = s.trave(b, c)
        s.cerniera_interna(b)
        s.incastro(a); s.carrello(c)
        s.carico_distribuito(t2, -10)
        r = s.risolvi()
        assert r.grado_iperstaticita == 0
        assert r.reazioni()[2]["Ry"] == pytest.approx(15.0)
        assert r.reazioni()[0]["M"] == pytest.approx(45.0)
        assert r.sollecitazioni(1)["M"][0] == pytest.approx(0.0, abs=1e-9)

    def test_svincolo_su_estremo_trave(self):
        s = Struttura()
        a, b = s.nodo(0, 0), s.nodo(5, 0)
        t = s.trave(a, b, cerniera_j=True)
        s.incastro(a); s.incastro(b)
        s.carico_distribuito(t, -8)
        r = s.risolvi()
        assert r.reazioni()[1]["M"] == pytest.approx(0.0, abs=1e-9)
        assert r.reazioni()[0]["M"] == pytest.approx(8 * 25 / 8)   # incastro-appoggio: qL²/8

    def test_carrello_inclinato_reazione_normale(self):
        s, t = appoggiata(6.0)
        s.carrello(1, angolo=30)
        s.carico_distribuito(t, -10)
        r = s.risolvi()
        rb = r.reazioni()[1]
        assert rb["R_parallela"] == pytest.approx(0.0, abs=1e-8)
        assert rb["Ry"] == pytest.approx(30.0)
        assert rb["Rx"] == pytest.approx(-30.0 * math.tan(math.radians(30)))
        assert r.verifica_equilibrio()["ok"]

    def test_doppio_pendolo(self):
        # incastro + doppio pendolo verticale in punta: M agli estremi = PL/2
        L, P = 4.0, 10.0
        s = Struttura()
        a, b = s.nodo(0, 0), s.nodo(L, 0)
        t = s.trave(a, b)
        s.incastro(a)
        s.doppio_pendolo(b, angolo=90)       # scorre in verticale, non ruota
        s.forza(b, Fy=-P)
        r = s.risolvi()
        assert r.reazioni()[1]["Ry"] == pytest.approx(0.0, abs=1e-9)
        assert abs(r.reazioni()[1]["M"]) == pytest.approx(P * L / 2)
        assert r.U[3 * b + 2] == pytest.approx(0.0, abs=1e-12)

    def test_molla(self):
        L, P, EIs, k = 4.0, 10.0, 5e3, 1e3
        s = Struttura()
        a, b = s.nodo(0, 0), s.nodo(L, 0)
        s.trave(a, b, EI=EIs)
        s.incastro(a); s.molla(b, ky=k); s.forza(b, Fy=-P)
        r = s.risolvi()
        kb = 3 * EIs / L**3
        assert r.U[4] == pytest.approx(-P / (kb + k))
        assert r.reazioni()[1]["Ry"] == pytest.approx(P * k / (kb + k))

    def test_cedimento(self):
        L, d, EIs = 5.0, 0.01, 2e4
        s = Struttura()
        a, b = s.nodo(0, 0), s.nodo(L, 0)
        t = s.trave(a, b, EI=EIs)
        s.incastro(a); s.incastro(b)
        s.cedimento(b, y=-d)
        r = s.risolvi()
        assert r.U[4] == pytest.approx(-d)
        assert r.sollecitazioni(t)["M"][0] == pytest.approx(-6 * EIs * d / L**2)

    def test_labile_con_cinematismo(self):
        s = Struttura()
        a, b = s.nodo(0, 0), s.nodo(6, 0)
        s.trave(a, b)
        s.carrello(a); s.carrello(b)
        with pytest.raises(StrutturaLabile) as e:
            s.risolvi()
        m = e.value.meccanismo
        assert abs(m[0]) == pytest.approx(1.0) and abs(m[3]) == pytest.approx(1.0)
        assert s.gradi() == (0, 1)

    def test_coppia_su_cerniera_labile(self):
        s = Struttura()
        a, b, c = s.nodo(0, 0), s.nodo(2, 0), s.nodo(4, 0)
        s.trave(a, b); s.trave(b, c)
        s.cerniera_interna(b)
        s.incastro(a); s.incastro(c)
        s.coppia(b, 5.0)
        with pytest.raises(StrutturaLabile):
            s.risolvi()


# ---------------------------------------------------------------------------
# Termici
# ---------------------------------------------------------------------------

class TestTermici:

    def test_uniforme_impedito(self):
        s = Struttura()
        a, b = s.nodo(0, 0), s.nodo(5, 0)
        t = s.trave(a, b, EA=1e6)
        s.cerniera(a); s.cerniera(b)
        s.carico_termico(t, alpha=1e-5, dT=30)
        r = s.risolvi()
        assert r.sollecitazioni(t)["N"][3] == pytest.approx(-300.0)

    def test_farfalla_incastrata(self):
        s = Struttura()
        a, b = s.nodo(0, 0), s.nodo(5, 0)
        t = s.trave(a, b, EI=2e4)
        s.incastro(a); s.incastro(b)
        s.carico_termico(t, alpha=1e-5, dT_farfalla=20, h=0.5)
        M = s.risolvi().sollecitazioni(t)["M"]
        assert np.allclose(M, -2e4 * 1e-5 * 20 / 0.5)

    def test_farfalla_isostatica_nessuna_sollecitazione(self):
        s, t = appoggiata(6.0)
        s.carico_termico(t, alpha=1e-5, dT_farfalla=20, h=0.5)
        r = s.risolvi()
        assert np.allclose(r.sollecitazioni(t)["M"], 0.0, atol=1e-9)
        # rotazione agli appoggi: κ L / 2 (concavità verso l'alto)
        assert r.U[2] == pytest.approx(-1e-5 * 20 / 0.5 * 6 / 2)


# ---------------------------------------------------------------------------
# Telai e grado di iperstaticità
# ---------------------------------------------------------------------------

class TestTelai:

    def test_portale_incastrato_equilibrio(self):
        s = Struttura()
        n = [s.nodo(0, 0), s.nodo(0, 4), s.nodo(6, 4), s.nodo(6, 0)]
        s.trave(n[0], n[1]); t = s.trave(n[1], n[2]); s.trave(n[2], n[3])
        s.incastro(n[0]); s.incastro(n[3])
        s.carico_distribuito(t, -20); s.forza(n[1], Fx=15)
        r = s.risolvi()
        assert r.grado_iperstaticita == 3
        assert r.verifica_equilibrio()["ok"]

    def test_arco_tre_cerniere(self):
        # simmetrico, carico uniforme sulla proiezione: spinta H = qL²/(8f)
        L, f, q = 8.0, 2.0, 10.0
        s = Struttura()
        a, c, b = s.nodo(0, 0), s.nodo(L / 2, f), s.nodo(L, 0)
        t1 = s.trave(a, c); t2 = s.trave(c, b)
        s.cerniera(a); s.cerniera(b); s.cerniera_interna(c)
        s.carico_distribuito(t1, -q, proiezione=True)
        s.carico_distribuito(t2, -q, proiezione=True)
        r = s.risolvi()
        assert r.grado_iperstaticita == 0
        assert r.reazioni()[a]["Rx"] == pytest.approx(q * L**2 / (8 * f))
        assert r.reazioni()[a]["Ry"] == pytest.approx(q * L / 2)

    def test_trave_continua_tre_appoggi(self):
        # due campate uguali, q uniforme: M appoggio centrale = -qL²/8
        L, q = 4.0, 10.0
        s = Struttura()
        a, b, c = s.nodo(0, 0), s.nodo(L, 0), s.nodo(2 * L, 0)
        t1 = s.trave(a, b); t2 = s.trave(b, c)
        s.cerniera(a); s.carrello(b); s.carrello(c)
        s.carico_distribuito(t1, -q); s.carico_distribuito(t2, -q)
        r = s.risolvi()
        assert r.grado_iperstaticita == 1
        assert r.sollecitazioni(t1)["M"][-1] == pytest.approx(-q * L**2 / 8)
        assert r.reazioni()[b]["Ry"] == pytest.approx(10 * q * L / 8)


# ---------------------------------------------------------------------------
# JSON e CSV
# ---------------------------------------------------------------------------

class TestIO:

    MODELLO = {
        "nodi": [{"x": 0, "y": 0}, {"x": 6, "y": 0}],
        "elementi": [{"i": 0, "j": 1, "tipo": "trave", "EI": 1e4, "EA": 1e6}],
        "vincoli": [{"nodo": 0, "tipo": "cerniera"}, {"nodo": 1, "tipo": "carrello"}],
        "carichi": [{"tipo": "distribuito", "elem": 0, "qi": -10, "direzione": "y"}],
    }

    def test_risolvi_dict(self):
        out = risolvi_dict(self.MODELLO)
        assert out["ok"]
        assert out["elementi"][0]["estremi"]["M"]["max"] == pytest.approx(45.0)
        json.dumps(out)  # serializzabile

    def test_risolvi_dict_labile(self):
        m = json.loads(json.dumps(self.MODELLO))
        m["vincoli"] = [{"nodo": 0, "tipo": "carrello"}, {"nodo": 1, "tipo": "carrello"}]
        out = risolvi_dict(m)
        assert not out["ok"] and out["labile"] and "meccanismo" in out

    def test_csv(self, tmp_path):
        from fem_travi import struttura_da_dict
        r = struttura_da_dict(self.MODELLO).risolvi()
        p = tmp_path / "r.csv"
        testo = esporta_csv(r, str(p))
        assert "REAZIONI VINCOLARI" in testo and "SOLLECITAZIONI" in testo
        assert "45" in testo
        assert p.exists()
