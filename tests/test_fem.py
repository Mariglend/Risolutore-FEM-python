"""
tests/test_fem.py
=================
Test unitari e di integrazione per il solver FEM.

Convenzioni di segno usate nel solver:
  nodo_fy:  val > 0 = forza verso il BASSO
  nodo_fx:  val > 0 = forza verso SINISTRA
  nodo_m:   val > 0 = momento ORARIO
  distribuiti: val > 0 = carico verso il BASSO
  Reazioni: R = K·U - F  (positivo = reazione verso l'alto/destra/antiorario)

Esegui con:
    pytest tests/ -v
"""

import math
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pytest

from fem_travi import Struttura, Nodo, Trave, Vincolo, Carico
from fem_travi.assembler import k_locale, matrice_rotazione


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def mk_appoggiata(L=6.0, q=10.0, EI=10000, EA=100000):
    s = Struttura()
    s.aggiungi_nodo(Nodo(0, 0)); s.aggiungi_nodo(Nodo(L, 0))
    s.aggiungi_trave(Trave(0, 1, EI=EI, EA=EA))
    s.aggiungi_vincolo(Vincolo(0, ux_fisso=True,  uy_fisso=True, phi_fisso=False))
    s.aggiungi_vincolo(Vincolo(1, ux_fisso=False, uy_fisso=True, phi_fisso=False))
    return s


# ---------------------------------------------------------------------------
# Test matrici elemento
# ---------------------------------------------------------------------------

class TestMatriciElemento:

    def test_k_locale_simmetria(self):
        K = k_locale(L=5.0, EI=10000, EA=100000)
        assert K.shape == (6, 6)
        assert np.allclose(K, K.T), "K locale deve essere simmetrica"

    def test_k_locale_valori_noti(self):
        K = k_locale(L=1.0, EI=1.0, EA=1.0)
        assert K[1, 1] == pytest.approx(12.0)
        assert K[2, 2] == pytest.approx(4.0)
        assert K[5, 5] == pytest.approx(4.0)
        assert K[2, 5] == pytest.approx(2.0)
        assert K[1, 4] == pytest.approx(-12.0)

    def test_k_locale_semidefinita_positiva(self):
        K = k_locale(L=4.0, EI=5000, EA=80000)
        eigvals = np.linalg.eigvalsh(K)
        assert np.all(eigvals >= -1e-8), "K locale non deve avere autovalori negativi"

    def test_matrice_rotazione_theta0(self):
        T = matrice_rotazione(0.0)
        assert np.allclose(T, np.eye(6)), "θ=0 → T = I"

    def test_matrice_rotazione_ortonormale(self):
        for theta in [0.0, math.pi/6, math.pi/4, math.pi/2, math.pi]:
            T = matrice_rotazione(theta)
            assert np.allclose(T @ T.T, np.eye(6)), f"T non ortonormale per θ={theta:.4f}"


# ---------------------------------------------------------------------------
# Test solver con soluzioni analitiche
# ---------------------------------------------------------------------------

class TestSolverAnalitico:

    def test_appoggiata_uniforme_reazioni(self):
        """Ry_A = Ry_B = q*L/2."""
        L, q = 6.0, 10.0
        s = mk_appoggiata(L=L, q=q)
        s.aggiungi_carico(Carico("uniforme", obj=0, val1=q))
        r = s.risolvi()
        rv = r.reazioni_vincolari()
        assert rv["Ry_N1"][0] == pytest.approx(q * L / 2, rel=1e-3)
        assert rv["Ry_N2"][0] == pytest.approx(q * L / 2, rel=1e-3)

    def test_appoggiata_uniforme_equilibrio(self):
        s = mk_appoggiata()
        s.aggiungi_carico(Carico("uniforme", obj=0, val1=10.0))
        r = s.risolvi()
        eq = r.verifica_equilibrio(tol=1e-4)
        assert eq["ok"], f"Equilibrio fallito: {eq}"

    def test_appoggiata_uniforme_rotazione(self):
        """φ_A = q*L³/(24*EI)."""
        L, q, EI = 6.0, 10.0, 10_000.0
        s = mk_appoggiata(L=L, q=q, EI=EI)
        s.aggiungi_carico(Carico("uniforme", obj=0, val1=q))
        r = s.risolvi()
        phi_an = q * L**3 / (24 * EI)
        assert abs(r.U[2]) == pytest.approx(phi_an, rel=1e-3)

    def test_appoggiata_freccia_centro(self):
        """δ_max = 5*q*L⁴/(384*EI) verificato con nodo al centro."""
        L, q, EI = 6.0, 10.0, 10_000.0
        s = Struttura()
        s.aggiungi_nodo(Nodo(0, 0)); s.aggiungi_nodo(Nodo(3, 0)); s.aggiungi_nodo(Nodo(6, 0))
        s.aggiungi_trave(Trave(0, 1, EI=EI, EA=100_000))
        s.aggiungi_trave(Trave(1, 2, EI=EI, EA=100_000))
        s.aggiungi_vincolo(Vincolo(0, ux_fisso=True,  uy_fisso=True, phi_fisso=False))
        s.aggiungi_vincolo(Vincolo(2, ux_fisso=False, uy_fisso=True, phi_fisso=False))
        s.aggiungi_carico(Carico("uniforme", obj=0, val1=q))
        s.aggiungi_carico(Carico("uniforme", obj=1, val1=q))
        r = s.risolvi()
        delta_an = 5 * q * L**4 / (384 * EI)
        assert abs(r.U[4]) == pytest.approx(delta_an, rel=5e-3)

    def test_sbalzo_forza_puntuale(self):
        """
        Trave a sbalzo L=4m, F=10kN verso il basso all'estremo.
        Ry = F, Mz = +F*L (antiorario), δ = F*L³/(3*EI).
        """
        L, F, EI = 4.0, 10.0, 5_000.0
        s = Struttura()
        s.aggiungi_nodo(Nodo(0, 0)); s.aggiungi_nodo(Nodo(L, 0))
        s.aggiungi_trave(Trave(0, 1, EI=EI, EA=100_000))
        s.aggiungi_vincolo(Vincolo(0, ux_fisso=True, uy_fisso=True, phi_fisso=True))
        s.aggiungi_carico(Carico("nodo_fy", obj=1, val1=F))  # val>0 = verso il basso
        r = s.risolvi()
        rv = r.reazioni_vincolari()
        assert rv["Ry_N1"][0] == pytest.approx(F, rel=1e-3)
        assert rv["Mz_N1"][0] == pytest.approx(F * L, rel=1e-3)  # antiorario
        assert abs(r.U[4])    == pytest.approx(F * L**3 / (3 * EI), rel=1e-3)
        assert r.verifica_equilibrio(tol=1e-4)["ok"]

    def test_triangolare_sx_reazioni(self):
        """R_A = q*L/3, R_B = q*L/6 per carico triangolare max a sinistra."""
        L, q = 6.0, 12.0
        s = mk_appoggiata(L=L)
        s.aggiungi_carico(Carico("triangolare_sx", obj=0, val1=q))
        r = s.risolvi()
        rv = r.reazioni_vincolari()
        assert rv["Ry_N1"][0] == pytest.approx(q * L / 3, rel=1e-3)
        assert rv["Ry_N2"][0] == pytest.approx(q * L / 6, rel=1e-3)
        assert r.verifica_equilibrio(tol=1e-4)["ok"]

    def test_triangolare_dx_reazioni(self):
        """R_A = q*L/6, R_B = q*L/3 per carico triangolare max a destra."""
        L, q = 6.0, 12.0
        s = mk_appoggiata(L=L)
        s.aggiungi_carico(Carico("triangolare_dx", obj=0, val1=q))
        r = s.risolvi()
        rv = r.reazioni_vincolari()
        assert rv["Ry_N1"][0] == pytest.approx(q * L / 6, rel=1e-3)
        assert rv["Ry_N2"][0] == pytest.approx(q * L / 3, rel=1e-3)
        assert r.verifica_equilibrio(tol=1e-4)["ok"]

    def test_colonna_verticale_forza_orizzontale(self):
        """Colonna verticale con forza orizzontale: |Rx| = F, |Mz| = F*L."""
        L, F = 5.0, 20.0
        s = Struttura()
        s.aggiungi_nodo(Nodo(0, 0)); s.aggiungi_nodo(Nodo(0, L))
        s.aggiungi_trave(Trave(0, 1, EI=15_000, EA=120_000))
        s.aggiungi_vincolo(Vincolo(0, ux_fisso=True, uy_fisso=True, phi_fisso=True))
        s.aggiungi_carico(Carico("nodo_fx", obj=1, val1=F))
        r = s.risolvi()
        rv = r.reazioni_vincolari()
        assert abs(rv["Rx_N1"][0]) == pytest.approx(F,     rel=1e-3)
        assert abs(rv["Mz_N1"][0]) == pytest.approx(F * L, rel=1e-3)
        assert r.verifica_equilibrio(tol=1e-4)["ok"]


# ---------------------------------------------------------------------------
# Test strutture più complesse
# ---------------------------------------------------------------------------

class TestStruttureComplesse:

    def test_trave_continua_equilibrio(self):
        """Trave continua su 3 appoggi con carico triangolare."""
        s = Struttura()
        s.aggiungi_nodo(Nodo(0, 0)); s.aggiungi_nodo(Nodo(4, 0)); s.aggiungi_nodo(Nodo(8, 0))
        s.aggiungi_trave(Trave(0, 1, EI=8_000, EA=80_000))
        s.aggiungi_trave(Trave(1, 2, EI=8_000, EA=80_000))
        s.aggiungi_vincolo(Vincolo(0, ux_fisso=True,  uy_fisso=True, phi_fisso=False))
        s.aggiungi_vincolo(Vincolo(1, ux_fisso=False, uy_fisso=True, phi_fisso=False))
        s.aggiungi_vincolo(Vincolo(2, ux_fisso=False, uy_fisso=True, phi_fisso=False))
        s.aggiungi_carico(Carico("triangolare_sx", obj=1, val1=12.0))
        r = s.risolvi()
        assert r.verifica_equilibrio(tol=1e-4)["ok"]

    def test_telaio_L_equilibrio(self):
        """Telaio a L: colonna + trave con carichi misti."""
        s = Struttura()
        s.aggiungi_nodo(Nodo(0, 0)); s.aggiungi_nodo(Nodo(5, 4)); s.aggiungi_nodo(Nodo(0, 4))
        s.aggiungi_trave(Trave(0, 2, EI=20_000, EA=150_000))
        s.aggiungi_trave(Trave(2, 1, EI=20_000, EA=150_000))
        s.aggiungi_vincolo(Vincolo(0, ux_fisso=True, uy_fisso=True, phi_fisso=True))
        s.aggiungi_vincolo(Vincolo(1, ux_fisso=False, uy_fisso=True, phi_fisso=False))
        s.aggiungi_carico(Carico("nodo_fx",  obj=2, val1=15.0))
        s.aggiungi_carico(Carico("uniforme", obj=1, val1=8.0))
        r = s.risolvi()
        assert r.verifica_equilibrio(tol=1e-3)["ok"]

    def test_trave_inclinata_equilibrio(self):
        """Trave inclinata a 45°: verifica equilibrio con carico verticale."""
        import math
        L = 5.0
        s = Struttura()
        s.aggiungi_nodo(Nodo(0, 0))
        s.aggiungi_nodo(Nodo(L * math.cos(math.pi/4), L * math.sin(math.pi/4)))
        s.aggiungi_trave(Trave(0, 1, EI=10_000, EA=100_000))
        s.aggiungi_vincolo(Vincolo(0, ux_fisso=True, uy_fisso=True, phi_fisso=True))
        s.aggiungi_vincolo(Vincolo(1, ux_fisso=True, uy_fisso=True, phi_fisso=False))
        s.aggiungi_carico(Carico("nodo_fy", obj=1, val1=10.0))
        r = s.risolvi()
        assert r.verifica_equilibrio(tol=1e-4)["ok"]


# ---------------------------------------------------------------------------
# Test validazione
# ---------------------------------------------------------------------------

class TestValidazione:

    def test_nodi_coincidenti_raise(self):
        s = Struttura()
        s.aggiungi_nodo(Nodo(0, 0)); s.aggiungi_nodo(Nodo(0, 0))
        s.aggiungi_vincolo(Vincolo(0, ux_fisso=True, uy_fisso=True, phi_fisso=True))
        with pytest.raises((ValueError, Exception)):
            s.aggiungi_trave(Trave(0, 1, EI=10000, EA=100000))
            s.risolvi()

    def test_struttura_labile_raise(self):
        s = Struttura()
        s.aggiungi_nodo(Nodo(0, 0)); s.aggiungi_nodo(Nodo(6, 0))
        s.aggiungi_trave(Trave(0, 1, EI=10000, EA=100000))
        s.aggiungi_vincolo(Vincolo(0, ux_fisso=False, uy_fisso=True, phi_fisso=False))
        with pytest.raises(np.linalg.LinAlgError):
            s.risolvi()

    def test_EI_negativo_raise(self):
        s = Struttura()
        s.aggiungi_nodo(Nodo(0, 0)); s.aggiungi_nodo(Nodo(3, 0))
        with pytest.raises(ValueError, match="EI"):
            s.aggiungi_trave(Trave(0, 1, EI=-1, EA=100000))

    def test_carico_tipo_sconosciuto(self):
        s = Struttura()
        s.aggiungi_nodo(Nodo(0, 0)); s.aggiungi_nodo(Nodo(3, 0))
        s.aggiungi_trave(Trave(0, 1, EI=10000, EA=100000))
        with pytest.raises(ValueError, match="sconosciuto"):
            s.aggiungi_carico(Carico("tipo_xxx", obj=0, val1=10))

    def test_nessun_nodo_raise(self):
        with pytest.raises(ValueError):
            Struttura().risolvi()
