"""
gui_fem_pro_v2.py  —  FEM Solver Pro v2
========================================
Solver FEM avanzato per strutture piane 2D.

Novità v2:
  - Vincoli completi: incastro, cerniera, carrello, biella, pattino,
    molla traslazionale/rotazionale, vincolo a terra
  - Svincoli interni alle estremità delle travi (cerniere interne)
  - Elementi Truss (asta, solo N, 2 GDL/nodo)
  - Database sezioni IPE/HEA/HEB/tubolari + sezioni custom
  - Database materiali: acciaio, alluminio, titanio, CFRP
  - Analisi modale: frequenze proprie + visualizzazione modi
  - Carichi termici α·ΔT
  - Finestra risultati avanzata separata

Requisiti:
    pip install numpy matplotlib reportlab scipy

Avvio:
    python gui_fem_pro_v2.py
"""

import sys, math, json, os, copy, time
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np

try:
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    from scipy.linalg import eigh
    HAS_SCI = True
except ImportError:
    HAS_SCI = False

try:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.platypus import (SimpleDocTemplate, Table, TableStyle,
                                    Paragraph, Spacer)
    from reportlab.lib.styles import getSampleStyleSheet
    HAS_RL = True
except ImportError:
    HAS_RL = False

# ── Palette ──────────────────────────────────────────────────────────────────
BG     = "#1a1b26"; PANEL  = "#24253a"; PANEL2 = "#1e1f2e"
BORDER = "#383a5c"; GRID   = "#252640"
ACCENT = "#7c6af7"; ACC2   = "#5dcaa5"; ACC3   = "#f7c06a"
DANGER = "#e46d6d"; WARN   = "#f7c06a"
TEXT   = "#c8cce8"; TEXT2  = "#7a7fa8"; TEXT3  = "#4a4f78"
BEAM   = "#7c6af7"; TRUSS  = "#f7c06a"; NODE   = "#c8cce8"
SEL    = "#f7c06a"; REACT  = "#5dcaa5"; LOAD   = "#f38ba8"
DEFORM = "#f38ba8"; MOMENT = "#cba6f7"; SHEAR  = "#89dceb"
NORMAL = "#a6e3a1"; SPRING = "#fab387"
SNAP_G = 0.5

# ── Database Materiali ────────────────────────────────────────────────────────
MATERIALS = {
    "Acciaio S235":   {"E": 210e6, "nu": 0.3,  "rho": 7850, "fy": 235e3, "alpha": 12e-6, "color": "#7c6af7"},
    "Acciaio S355":   {"E": 210e6, "nu": 0.3,  "rho": 7850, "fy": 355e3, "alpha": 12e-6, "color": "#9a8af7"},
    "Alluminio 6061": {"E":  69e6, "nu": 0.33, "rho": 2700, "fy": 276e3, "alpha": 23e-6, "color": "#5dcaa5"},
    "Titanio Ti-6Al-4V": {"E": 114e6, "nu": 0.34, "rho": 4430, "fy": 880e3, "alpha": 8.6e-6, "color": "#f7c06a"},
    "CFRP (0°)":      {"E": 135e6, "nu": 0.3,  "rho": 1600, "fy": 1500e3,"alpha": 0.5e-6, "color": "#f38ba8"},
    "Custom":         {"E": 210e6, "nu": 0.3,  "rho": 7850, "fy": 235e3, "alpha": 12e-6, "color": "#cdd6f4"},
}

# ── Database Sezioni ─────────────────────────────────────────────────────────
def _ipe(h, b, tf, tw, r=0):
    """Calcola proprietà sezione IPE."""
    A  = 2*b*tf + (h-2*tf)*tw
    I  = b*h**3/12 - (b-tw)*(h-2*tf)**3/12
    W  = I/(h/2)
    return {"A": A, "I": I, "W": W, "h": h, "b": b}

SECTIONS = {
    # IPE
    "IPE 80":  {**_ipe(0.080,0.046,0.0038,0.0038), "type":"IPE"},
    "IPE 100": {**_ipe(0.100,0.055,0.0041,0.0041), "type":"IPE"},
    "IPE 120": {**_ipe(0.120,0.064,0.0044,0.0044), "type":"IPE"},
    "IPE 160": {**_ipe(0.160,0.082,0.0050,0.0050), "type":"IPE"},
    "IPE 200": {**_ipe(0.200,0.100,0.0085,0.0056), "type":"IPE"},
    "IPE 240": {**_ipe(0.240,0.120,0.0098,0.0062), "type":"IPE"},
    "IPE 270": {**_ipe(0.270,0.135,0.0102,0.0066), "type":"IPE"},
    "IPE 300": {**_ipe(0.300,0.150,0.0107,0.0071), "type":"IPE"},
    "IPE 360": {**_ipe(0.360,0.170,0.0127,0.0080), "type":"IPE"},
    "IPE 400": {**_ipe(0.400,0.180,0.0135,0.0086), "type":"IPE"},
    "IPE 500": {**_ipe(0.500,0.200,0.0160,0.0102), "type":"IPE"},
    # HEA
    "HEA 100": {**_ipe(0.096,0.100,0.0080,0.0050), "type":"HEA"},
    "HEA 120": {**_ipe(0.114,0.120,0.0080,0.0050), "type":"HEA"},
    "HEA 160": {**_ipe(0.152,0.160,0.0090,0.0060), "type":"HEA"},
    "HEA 200": {**_ipe(0.190,0.200,0.0100,0.0065), "type":"HEA"},
    "HEA 240": {**_ipe(0.230,0.240,0.0120,0.0075), "type":"HEA"},
    "HEA 300": {**_ipe(0.290,0.300,0.0140,0.0085), "type":"HEA"},
    # HEB
    "HEB 100": {**_ipe(0.100,0.100,0.0100,0.0060), "type":"HEB"},
    "HEB 160": {**_ipe(0.160,0.160,0.0130,0.0080), "type":"HEB"},
    "HEB 200": {**_ipe(0.200,0.200,0.0150,0.0090), "type":"HEB"},
    "HEB 300": {**_ipe(0.300,0.300,0.0190,0.0110), "type":"HEB"},
    # Tubolare quadrato
    "CHS 60x4":  {"A":6.97e-4,  "I":2.59e-7,  "W":8.63e-6,  "h":0.060,"b":0.060,"type":"CHS"},
    "CHS 100x5": {"A":15.1e-4,  "I":1.39e-6,  "W":2.78e-5,  "h":0.100,"b":0.100,"type":"CHS"},
    "RHS 100x60x5":{"A":14.4e-4,"I":9.10e-7,  "W":1.82e-5,  "h":0.100,"b":0.060,"type":"RHS"},
    "Custom":    {"A":1e-3,     "I":1e-5,     "W":2e-4,     "h":0.200,"b":0.100,"type":"Custom"},
}

# ── Libreria strutture ────────────────────────────────────────────────────────
LIBRERIA = {
    "Trave appoggiata": {
        "nodi":    [{"x":0,"y":0},{"x":6,"y":0}],
        "elementi":[{"tipo":"trave","i":0,"j":1,"mat":"Acciaio S235","sez":"IPE 300",
                     "svincolo_i":False,"svincolo_j":False}],
        "vincoli": [{"nodo":0,"tipo":"cerniera","kx":0,"ky":0,"kphi":0},
                    {"nodo":1,"tipo":"carrello","kx":0,"ky":0,"kphi":0}],
        "carichi": [{"tipo":"uniforme","elem":0,"val":10,"dT":0}],
    },
    "Trave a sbalzo": {
        "nodi":    [{"x":0,"y":0},{"x":5,"y":0}],
        "elementi":[{"tipo":"trave","i":0,"j":1,"mat":"Acciaio S355","sez":"IPE 240",
                     "svincolo_i":False,"svincolo_j":False}],
        "vincoli": [{"nodo":0,"tipo":"incastro","kx":0,"ky":0,"kphi":0}],
        "carichi": [{"tipo":"Fy","nodo":1,"val":30,"dT":0}],
    },
    "Portale": {
        "nodi":    [{"x":0,"y":0},{"x":0,"y":4},{"x":6,"y":4},{"x":6,"y":0}],
        "elementi":[
            {"tipo":"trave","i":0,"j":1,"mat":"Acciaio S235","sez":"HEB 200","svincolo_i":False,"svincolo_j":False},
            {"tipo":"trave","i":1,"j":2,"mat":"Acciaio S235","sez":"IPE 300","svincolo_i":False,"svincolo_j":False},
            {"tipo":"trave","i":2,"j":3,"mat":"Acciaio S235","sez":"HEB 200","svincolo_i":False,"svincolo_j":False},
        ],
        "vincoli": [{"nodo":0,"tipo":"incastro","kx":0,"ky":0,"kphi":0},
                    {"nodo":3,"tipo":"incastro","kx":0,"ky":0,"kphi":0}],
        "carichi": [{"tipo":"uniforme","elem":1,"val":20,"dT":0},
                    {"tipo":"Fx","nodo":1,"val":15,"dT":0}],
    },
    "Capriata Pratt": {
        "nodi": [{"x":0,"y":0},{"x":2,"y":0},{"x":4,"y":0},{"x":6,"y":0},
                 {"x":1,"y":2},{"x":3,"y":2},{"x":5,"y":2}],
        "elementi":[
            {"tipo":"truss","i":0,"j":1,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
            {"tipo":"truss","i":1,"j":2,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
            {"tipo":"truss","i":2,"j":3,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
            {"tipo":"truss","i":4,"j":5,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
            {"tipo":"truss","i":5,"j":6,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
            {"tipo":"truss","i":0,"j":4,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
            {"tipo":"truss","i":1,"j":4,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
            {"tipo":"truss","i":1,"j":5,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
            {"tipo":"truss","i":2,"j":5,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
            {"tipo":"truss","i":2,"j":6,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
            {"tipo":"truss","i":3,"j":6,"mat":"Acciaio S235","sez":"CHS 60x4","svincolo_i":False,"svincolo_j":False},
        ],
        "vincoli": [{"nodo":0,"tipo":"cerniera","kx":0,"ky":0,"kphi":0},
                    {"nodo":3,"tipo":"carrello","kx":0,"ky":0,"kphi":0}],
        "carichi": [{"tipo":"Fy","nodo":4,"val":20,"dT":0},
                    {"tipo":"Fy","nodo":5,"val":20,"dT":0},
                    {"tipo":"Fy","nodo":6,"val":20,"dT":0}],
    },
}

# ── FEM Core ─────────────────────────────────────────────────────────────────

def elem_geo(ni, nj):
    dx,dy = nj['x']-ni['x'], nj['y']-ni['y']
    L = math.hypot(dx,dy)
    theta = math.atan2(dy,dx)
    return L, theta

def k_trave_loc(L, EI, EA):
    """Matrice rigidezza locale 6x6 trave Euler-Bernoulli."""
    a=EA/L; b=12*EI/L**3; c=6*EI/L**2; d=4*EI/L; e=2*EI/L
    K=np.zeros((6,6))
    K[0,0]=a;  K[0,3]=-a; K[3,0]=-a; K[3,3]=a
    K[1,1]=b;  K[1,2]=c;  K[1,4]=-b; K[1,5]=c
    K[2,1]=c;  K[2,2]=d;  K[2,4]=-c; K[2,5]=e
    K[4,1]=-b; K[4,2]=-c; K[4,4]=b;  K[4,5]=-c
    K[5,1]=c;  K[5,2]=e;  K[5,4]=-c; K[5,5]=d
    return K

def k_trave_loc_svincolo(L, EI, EA, sv_i, sv_j):
    """Matrice rigidezza con svincoli interni (cerniere alle estremità)."""
    K = k_trave_loc(L, EI, EA)
    # svincolo al nodo i: phi_i libero → condensazione statica
    if sv_i and not sv_j:
        # condensa phi_i (dof locale 2)
        k22 = K[2,2]
        if abs(k22) > 1e-12:
            for r in range(6):
                for c in range(6):
                    if r != 2 and c != 2:
                        K[r,c] -= K[r,2]*K[2,c]/k22
            K[2,:] = 0; K[:,2] = 0
    elif sv_j and not sv_i:
        # condensa phi_j (dof locale 5)
        k55 = K[5,5]
        if abs(k55) > 1e-12:
            for r in range(6):
                for c in range(6):
                    if r != 5 and c != 5:
                        K[r,c] -= K[r,5]*K[5,c]/k55
            K[5,:] = 0; K[:,5] = 0
    elif sv_i and sv_j:
        # condensa entrambi con eliminazione simultanea
        dofs_cond = [2, 5]
        for dc in dofs_cond:
            kdd = K[dc,dc]
            if abs(kdd) > 1e-12:
                for r in range(6):
                    for c in range(6):
                        if r != dc and c != dc:
                            K[r,c] -= K[r,dc]*K[dc,c]/kdd
                K[dc,:] = 0; K[:,dc] = 0
    return K

def T_trave(theta):
    """Matrice rotazione 6x6."""
    c,s = math.cos(theta), math.sin(theta)
    T = np.zeros((6,6))
    T[0,0]=c; T[0,1]=s; T[1,0]=-s; T[1,1]=c; T[2,2]=1
    T[3,3]=c; T[3,4]=s; T[4,3]=-s; T[4,4]=c; T[5,5]=1
    return T

def k_truss_loc(L, EA):
    """Matrice rigidezza locale 4x4 asta (truss)."""
    a = EA/L
    return np.array([[ a, 0,-a, 0],
                     [ 0, 0, 0, 0],
                     [-a, 0, a, 0],
                     [ 0, 0, 0, 0]])

def T_truss(theta):
    """Matrice rotazione 4x4 per truss."""
    c,s = math.cos(theta), math.sin(theta)
    T = np.zeros((4,4))
    T[0,0]=c; T[0,1]=s; T[1,0]=-s; T[1,1]=c
    T[2,2]=c; T[2,3]=s; T[3,2]=-s; T[3,3]=c
    return T

def m_trave_loc(L, rho, A):
    """Matrice di massa consistente 6x6 (trave)."""
    m = rho*A*L/420
    return m*np.array([
        [140,  0,    0,   70,  0,    0  ],
        [  0,156,   22*L,  0, 54,  -13*L],
        [  0, 22*L, 4*L**2,0, 13*L,-3*L**2],
        [ 70,  0,    0,  140,  0,    0  ],
        [  0, 54,  13*L,  0, 156, -22*L],
        [  0,-13*L,-3*L**2,0,-22*L, 4*L**2],
    ])

def m_truss_loc(L, rho, A):
    """Matrice di massa consistente 4x4 (truss)."""
    m = rho*A*L/6
    return np.array([
        [2,0,1,0],
        [0,2,0,1],
        [1,0,2,0],
        [0,1,0,2],
    ])*m

def get_EI_EA(el, sez_db, mat_db):
    mat = mat_db.get(el.get('mat','Acciaio S235'), mat_db['Acciaio S235'])
    sez = sez_db.get(el.get('sez','Custom'), sez_db['Custom'])
    E = mat['E']; I = sez['I']; A = sez['A']
    return E*I, E*A, E, A, mat, sez

def dofs_elemento(el, idx):
    """Restituisce lista GDL globali e tipo."""
    i,j = el['i'], el['j']
    if el['tipo'] == 'truss':
        return [i*3, i*3+1, j*3, j*3+1], 4
    else:
        return [i*3, i*3+1, i*3+2, j*3, j*3+1, j*3+2], 6

def assembla_sistema(nodi, elementi, vincoli, carichi):
    """
    Assembla K, M, F globali.
    Ritorna K, M, F, is_fixed, n_dof.
    """
    n = len(nodi); ndof = 3*n
    K = np.zeros((ndof,ndof)); M = np.zeros((ndof,ndof)); F = np.zeros(ndof)

    for el in elementi:
        i,j = el['i'], el['j']
        ni,nj = nodi[i], nodi[j]
        L,theta = elem_geo(ni,nj)
        if L < 1e-9: continue
        EI,EA,E,A,mat,sez = get_EI_EA(el, SECTIONS, MATERIALS)
        rho = mat['rho']

        if el['tipo'] == 'truss':
            Kl = k_truss_loc(L, EA)
            Tv = T_truss(theta)
            Kg = Tv.T @ Kl @ Tv
            Ml = m_truss_loc(L, rho, A)
            Mg = Tv.T @ Ml @ Tv
            dofs = [i*3, i*3+1, j*3, j*3+1]
            nd = 4
        else:
            sv_i = el.get('svincolo_i', False)
            sv_j = el.get('svincolo_j', False)
            Kl = k_trave_loc_svincolo(L, EI, EA, sv_i, sv_j)
            Tv = T_trave(theta)
            Kg = Tv.T @ Kl @ Tv
            Ml = m_trave_loc(L, rho, A)
            Mg = Tv.T @ Ml @ Tv
            dofs = [i*3, i*3+1, i*3+2, j*3, j*3+1, j*3+2]
            nd = 6

        for a in range(nd):
            for b in range(nd):
                K[dofs[a], dofs[b]] += Kg[a,b]
                M[dofs[a], dofs[b]] += Mg[a,b]

        # Carico termico
        alpha = mat['alpha']
        for c in carichi:
            if c['tipo'] in ('uniforme','triang_sx','triang_dx'):
                if c.get('elem') == elementi.index(el) and abs(c.get('dT',0)) > 1e-9:
                    dT = c['dT']; fT = np.zeros(nd)
                    if el['tipo'] == 'truss':
                        N_th = E*A*alpha*dT
                        fT[0] = -N_th; fT[2] = N_th
                    else:
                        N_th = E*A*alpha*dT
                        fT[0] = -N_th; fT[3] = N_th
                    fT_g = Tv.T @ fT
                    for a,da in enumerate(dofs): F[da] += fT_g[a]

    # Carichi esterni
    for c in carichi:
        tp = c['tipo']
        if tp == 'Fy':
            if c.get('nodo') is not None and c['nodo']<len(nodi):
                F[c['nodo']*3+1] -= c['val']
        elif tp == 'Fx':
            if c.get('nodo') is not None and c['nodo']<len(nodi):
                F[c['nodo']*3] -= c['val']
        elif tp == 'M':
            if c.get('nodo') is not None and c['nodo']<len(nodi):
                F[c['nodo']*3+2] -= c['val']
        elif tp in ('uniforme','triang_sx','triang_dx'):
            ei = c.get('elem')
            if ei is None or ei >= len(elementi): continue
            el = elementi[ei]
            if el['tipo'] == 'truss': continue
            i,j = el['i'],el['j']
            L,theta = elem_geo(nodi[i],nodi[j])
            if L<1e-9: continue
            q = c['val']; fne = np.zeros(6)
            if tp=='uniforme':
                fne[1]=-q*L/2; fne[2]=-q*L**2/12
                fne[4]=-q*L/2; fne[5]= q*L**2/12
            elif tp=='triang_sx':
                fne[1]=-7*q*L/20; fne[2]=-q*L**2/20
                fne[4]=-3*q*L/20; fne[5]= q*L**2/30
            else:
                fne[1]=-3*q*L/20; fne[2]=-q*L**2/30
                fne[4]=-7*q*L/20; fne[5]= q*L**2/20
            Tv = T_trave(theta); fg = Tv.T @ fne
            dofs = [i*3,i*3+1,i*3+2,j*3,j*3+1,j*3+2]
            for a,da in enumerate(dofs): F[da] += fg[a]

    # Vincoli con molle
    is_fixed = np.zeros(ndof, dtype=bool)
    for v in vincoli:
        nd = v['nodo']; tp = v['tipo']
        kx   = v.get('kx',0)
        ky   = v.get('ky',0)
        kphi = v.get('kphi',0)

        if tp == 'incastro':
            is_fixed[nd*3]=True; is_fixed[nd*3+1]=True; is_fixed[nd*3+2]=True
        elif tp == 'cerniera':
            is_fixed[nd*3]=True; is_fixed[nd*3+1]=True
        elif tp == 'carrello':
            is_fixed[nd*3+1]=True
        elif tp == 'pattino':      # blocca solo ux
            is_fixed[nd*3]=True
        elif tp == 'biella':       # blocca solo nella direzione dell'asta
            ang = v.get('angolo',0)
            # implementato via penalizzazione direzionale
            c2,s2 = math.cos(ang)**2, math.sin(ang)**2
            cs   = math.cos(ang)*math.sin(ang)
            pen  = 1e12
            K[nd*3,   nd*3]   += pen*c2
            K[nd*3,   nd*3+1] += pen*cs
            K[nd*3+1, nd*3]   += pen*cs
            K[nd*3+1, nd*3+1] += pen*s2
        elif tp == 'incastro_scorrevole':  # blocca uy e phi, libero ux
            is_fixed[nd*3+1]=True; is_fixed[nd*3+2]=True
        elif tp == 'molla':
            if kx   > 0: K[nd*3,   nd*3]   += kx
            if ky   > 0: K[nd*3+1, nd*3+1] += ky
            if kphi > 0: K[nd*3+2, nd*3+2] += kphi
        elif tp == 'carrello_x':   # blocca ux
            is_fixed[nd*3]=True
        elif tp == 'carrello_inclinato':
            ang = v.get('angolo',0)
            pen = 1e12
            cx,cy = math.cos(ang), math.sin(ang)
            K[nd*3,   nd*3]   += pen*cx**2
            K[nd*3,   nd*3+1] += pen*cx*cy
            K[nd*3+1, nd*3]   += pen*cx*cy
            K[nd*3+1, nd*3+1] += pen*cy**2

    return K, M, F, is_fixed

def solve_statico(K, F, is_fixed):
    ndof = K.shape[0]
    free  = [i for i in range(ndof) if not is_fixed[i]]
    fixed = [i for i in range(ndof) if     is_fixed[i]]
    if not free: raise ValueError("Struttura completamente bloccata")
    Kff = K[np.ix_(free,free)]; Ff = F[free]
    try:
        Uf = np.linalg.solve(Kff, Ff)
    except np.linalg.LinAlgError:
        raise ValueError("Sistema singolare: struttura labile")
    U = np.zeros(ndof)
    for li,gi in enumerate(free): U[gi] = Uf[li]
    R = K@U - F
    return U, R, free, fixed

def solve_modale(K, M, is_fixed, n_modi=6):
    """Analisi agli autovalori: [K]{U} = ω²[M]{U}."""
    ndof = K.shape[0]
    free = [i for i in range(ndof) if not is_fixed[i]]
    Kff = K[np.ix_(free,free)]
    Mff = M[np.ix_(free,free)]
    # regolarizza M
    Mff += np.eye(len(free))*1e-12
    n_max = min(n_modi, len(free)-1)
    if n_max < 1: raise ValueError("Troppi pochi GDL liberi per analisi modale")
    if HAS_SCI:
        vals, vecs = eigh(Kff, Mff, subset_by_index=[0,n_max-1])
    else:
        vals, vecs = np.linalg.eigh(Kff)
        vals = vals[:n_max]; vecs = vecs[:,:n_max]
    omega = np.sqrt(np.maximum(vals, 0))
    freq  = omega / (2*math.pi)
    # estendi autovettori a ndof
    modes = []
    for k in range(n_max):
        phi = np.zeros(ndof)
        for li,gi in enumerate(free): phi[gi] = vecs[li,k]
        modes.append(phi)
    return freq, omega, modes

def internal_forces(el, nodi, U):
    i,j = el['i'],el['j']
    L,theta = elem_geo(nodi[i],nodi[j])
    if L<1e-9: return None
    EI,EA,E,A,mat,sez = get_EI_EA(el, SECTIONS, MATERIALS)
    T = T_trave(theta)
    u_g = np.array([U[i*3],U[i*3+1],U[i*3+2],U[j*3],U[j*3+1],U[j*3+2]])
    u_l = T @ u_g
    ul_i,vl_i,phi_i = u_l[0],u_l[1],u_l[2]
    ul_j,vl_j,phi_j = u_l[3],u_l[4],u_l[5]
    steps = 50
    xs = np.linspace(0,L,steps+1)
    M_vals,V_vals,N_vals = [],[],[]
    for xi in xs:
        s = xi/L
        d2N1=(-6+12*s)/L**2; d2N2=(-4+6*s)/L
        d2N3=(6-12*s)/L**2;  d2N4=(-2+6*s)/L
        M = EI*(d2N1*vl_i+d2N2*phi_i+d2N3*vl_j+d2N4*phi_j)
        d3 = 12/L**3*(vl_i-vl_j)+(6/L**2)*(phi_i-phi_j)
        V = -EI*d3
        N = EA*(ul_j-ul_i)/L
        M_vals.append(M); V_vals.append(V); N_vals.append(N)
    return xs, np.array(M_vals), np.array(V_vals), np.array(N_vals)

def truss_force(el, nodi, U):
    i,j = el['i'],el['j']
    L,theta = elem_geo(nodi[i],nodi[j])
    if L<1e-9: return 0
    _,EA,E,A,_,_ = get_EI_EA(el, SECTIONS, MATERIALS)
    c,s = math.cos(theta),math.sin(theta)
    u = np.array([U[i*3],U[i*3+1],U[j*3],U[j*3+1]])
    dl = c*(u[2]-u[0]) + s*(u[3]-u[1])
    return EA*dl/L

# ── Widget helpers ────────────────────────────────────────────────────────────

def mk_btn(parent, text, cmd, bg=ACCENT, fg=TEXT, **kw):
    return tk.Button(parent, text=text, command=cmd, bg=bg, fg=fg,
                     relief="flat", activebackground=ACC2, activeforeground=BG,
                     font=("Segoe UI",10,"bold"), padx=8, pady=4,
                     cursor="hand2", **kw)

def mk_lbl(parent, text, fg=TEXT2, size=9, bold=False):
    return tk.Label(parent, text=text, bg=PANEL, fg=fg,
                    font=("Segoe UI", size, "bold" if bold else "normal"))

def hsep(parent, bg=PANEL):
    tk.Frame(parent, bg=BORDER, height=1).pack(fill="x", pady=5)

def sec_lbl(parent, text):
    tk.Label(parent, text=text, bg=PANEL, fg=TEXT2,
             font=("Segoe UI",9,"bold")).pack(anchor="w", padx=12, pady=(8,2))

# ── App ───────────────────────────────────────────────────────────────────────

class FEMProV2(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FEM Solver Pro v2")
        self.configure(bg=BG)
        self.minsize(1200,720); self.geometry("1500x860")

        # struttura
        self.nodi     = []
        self.elementi = []   # ogni el ha 'tipo': 'trave' | 'truss'
        self.vincoli  = []
        self.carichi  = []

        # risultati statici
        self.U=None; self.R=None; self.fixed=[]; self.K=None; self.F=None; self.M_mat=None
        # risultati modali
        self.freq=None; self.modi=None; self.modo_idx=0

        # undo
        self._undo=[]; self._redo=[]

        # UI vars
        self.mode        = tk.StringVar(value="nodo")
        self.el_tipo     = tk.StringVar(value="trave")
        self.mat_var     = tk.StringVar(value="Acciaio S235")
        self.sez_var     = tk.StringVar(value="IPE 300")
        self.vtype       = tk.StringVar(value="incastro")
        self.v_ang       = tk.StringVar(value="0")
        self.v_kx        = tk.StringVar(value="0")
        self.v_ky        = tk.StringVar(value="0")
        self.v_kphi      = tk.StringVar(value="0")
        self.sv_i        = tk.BooleanVar(value=False)
        self.sv_j        = tk.BooleanVar(value=False)
        self.ctype       = tk.StringVar(value="Fy")
        self.cval        = tk.StringVar(value="10")
        self.cdT         = tk.StringVar(value="0")
        self.show_def    = tk.BooleanVar(value=False)
        self.show_mvn    = tk.BooleanVar(value=False)
        self.show_heat   = tk.BooleanVar(value=False)
        self.show_labels = tk.BooleanVar(value=True)
        self.snap_g      = tk.BooleanVar(value=True)
        self.snap_n      = tk.BooleanVar(value=True)
        self.def_scale   = tk.DoubleVar(value=1.0)
        self.n_modi_var  = tk.IntVar(value=6)
        self.modo_var    = tk.IntVar(value=0)
        self.show_modo   = tk.BooleanVar(value=False)

        # canvas
        self.cam_x=0.0; self.cam_y=0.0; self.cam_scale=80.0
        self.sel_nodo=None; self.sel_elem=None
        self.drag_from=None; self.drag_pos=None
        self.dragging=False; self.pan_start=None
        self.anim_run=False; self.anim_ph=0.0
        self.current_file=None

        self._build_ui()
        self._reset_view()
        self.bind_all("<Control-z>", lambda e: self._undo_fn())
        self.bind_all("<Control-y>", lambda e: self._redo_fn())
        self.bind_all("<Control-s>", lambda e: self._save())
        self.bind_all("<Control-o>", lambda e: self._load())
        self.bind_all("<Control-n>", lambda e: self._new())

    # ── UI ───────────────────────────────────────────────────────────────────

    def _build_ui(self):
        self.columnconfigure(1,weight=1)
        self.rowconfigure(1,weight=1)
        self._build_menu()
        self._build_topbar()
        self._build_left()
        self._build_center()
        self._build_right()

    def _build_menu(self):
        mb=tk.Menu(self,bg=PANEL2,fg=TEXT,activebackground=ACCENT,activeforeground=BG,relief="flat")
        self.config(menu=mb)
        def m(label, items):
            mn=tk.Menu(mb,tearoff=0,bg=PANEL2,fg=TEXT,activebackground=ACCENT,activeforeground=BG)
            mb.add_cascade(label=label,menu=mn)
            for it in items:
                if it is None: mn.add_separator()
                else: mn.add_command(label=it[0],command=it[1])
        m("File",[("🗋 Nuovo  Ctrl+N",self._new),("📂 Apri  Ctrl+O",self._load),
                  ("💾 Salva  Ctrl+S",self._save),("Salva con nome...",self._save_as),
                  None,("📄 Esporta PDF",self._export_pdf),("Esci",self.quit)])
        m("Modifica",[("↩ Undo  Ctrl+Z",self._undo_fn),("↪ Redo  Ctrl+Y",self._redo_fn),
                      None,("🗑 Cancella tutto",self._clear_all)])
        lm=tk.Menu(mb,tearoff=0,bg=PANEL2,fg=TEXT,activebackground=ACCENT,activeforeground=BG)
        mb.add_cascade(label="Libreria",menu=lm)
        for nome in LIBRERIA:
            lm.add_command(label=nome,command=lambda n=nome:self._load_preset(n))
        m("Analisi",[("▶ Analisi statica",self._solve),
                     ("🎵 Analisi modale",self._solve_modale),
                     None,("📊 Risultati avanzati",self._open_results_window)])
        m("Vista",[("Deformata",self._toggle_def),("Diagrammi M/V/N",self._toggle_mvn),
                   ("Mappa tensioni",self._toggle_heat),None,("Reset vista",self._reset_view)])

    def _build_topbar(self):
        tb=tk.Frame(self,bg=PANEL2,height=44)
        tb.grid(row=0,column=0,columnspan=3,sticky="ew"); tb.grid_propagate(False)
        def tbtn(t,cmd,col=PANEL2):
            b=tk.Button(tb,text=t,command=cmd,bg=col,fg=TEXT,relief="flat",
                        font=("Segoe UI",9),padx=8,pady=6,cursor="hand2",
                        activebackground=ACCENT,activeforeground=BG)
            b.pack(side="left",padx=1); return b
        tbtn("🗋",self._new); tbtn("📂",self._load); tbtn("💾",self._save); tbtn("📄",self._export_pdf)
        tk.Frame(tb,bg=BORDER,width=1).pack(side="left",fill="y",padx=4)
        tbtn("↩",self._undo_fn); tbtn("↪",self._redo_fn)
        tk.Frame(tb,bg=BORDER,width=1).pack(side="left",fill="y",padx=4)
        tbtn("⊕",lambda:self._zoom(1.15)); tbtn("⊖",lambda:self._zoom(0.87)); tbtn("⌂",self._reset_view)
        tk.Frame(tb,bg=BORDER,width=1).pack(side="left",fill="y",padx=4)
        self._mbns={}
        for txt,val in [("✦ Nodo","nodo"),("─ Elem","trave"),("⬛ Vincolo","vincolo"),
                        ("↓ Carico","carico"),("✥ Sposta","sposta")]:
            b=tbtn(txt,lambda v=val:self._set_mode(v))
            self._mbns[val]=b
        tk.Frame(tb,bg=BORDER,width=1).pack(side="left",fill="y",padx=4)
        tbtn("▶ RISOLVI",self._solve,col="#2a2a5a")
        tbtn("🎵 MODALE",self._solve_modale,col="#1a3a2a")
        tbtn("📊 RISULTATI",self._open_results_window,col="#3a1a2a")
        self.anim_btn=tbtn("🎬 Anima",self._toggle_anim)
        tk.Label(tb,text="Scala def:",bg=PANEL2,fg=TEXT2,font=("Segoe UI",8)).pack(side="left",padx=(8,2))
        tk.Scale(tb,variable=self.def_scale,from_=0.1,to=20,resolution=0.1,
                 orient="horizontal",bg=PANEL2,fg=TEXT,highlightthickness=0,
                 troughcolor=PANEL,length=70,showvalue=False,
                 command=lambda v:self._redraw()).pack(side="left")
        self.file_lbl=tk.Label(tb,text="senza titolo",bg=PANEL2,fg=TEXT2,font=("Segoe UI",9))
        self.file_lbl.pack(side="right",padx=12)
        self._upd_mode_btns()

    def _set_mode(self,v):
        self.mode.set(v); self.sel_nodo=None; self.sel_elem=None
        self._upd_mode_btns(); self._redraw()

    def _upd_mode_btns(self):
        m=self.mode.get()
        for v,b in self._mbns.items():
            b.config(bg=ACCENT if v==m else PANEL2,fg=BG if v==m else TEXT)

    def _build_left(self):
        frame=tk.Frame(self,bg=PANEL,width=300)
        frame.grid(row=1,column=0,sticky="nsew"); frame.grid_propagate(False)
        cv=tk.Canvas(frame,bg=PANEL,highlightthickness=0)
        sb=tk.Scrollbar(frame,orient="vertical",command=cv.yview,bg=PANEL)
        cv.configure(yscrollcommand=sb.set)
        sb.pack(side="right",fill="y"); cv.pack(side="left",fill="both",expand=True)
        inn=tk.Frame(cv,bg=PANEL)
        win=cv.create_window((0,0),window=inn,anchor="nw")
        inn.bind("<Configure>",lambda e:cv.configure(scrollregion=cv.bbox("all")))
        cv.bind("<Configure>",lambda e:cv.itemconfig(win,width=e.width))
        p=inn

        # Tipo elemento
        sec_lbl(p,"TIPO ELEMENTO")
        for txt,val in [("— Trave (M,V,N)","trave"),("⊿ Asta/Truss (N)","truss")]:
            tk.Radiobutton(p,text=txt,variable=self.el_tipo,value=val,
                           bg=PANEL,fg=TEXT,selectcolor=ACCENT,
                           activebackground=PANEL,font=("Segoe UI",10),
                           indicatoron=False,relief="flat",padx=8,pady=3,cursor="hand2"
                           ).pack(fill="x",padx=12,pady=1)
        hsep(p)

        # Materiale
        sec_lbl(p,"MATERIALE")
        self.mat_cb=ttk.Combobox(p,textvariable=self.mat_var,
                                  values=list(MATERIALS.keys()),state="readonly",
                                  font=("Segoe UI",9))
        self.mat_cb.pack(fill="x",padx=12,pady=2)
        self.mat_info=tk.Label(p,text="",bg=PANEL,fg=TEXT2,font=("Segoe UI",8),
                                justify="left",anchor="w")
        self.mat_info.pack(fill="x",padx=12)
        self.mat_cb.bind("<<ComboboxSelected>>",self._upd_mat_info)
        self._upd_mat_info()

        # Sezione
        sec_lbl(p,"SEZIONE")
        self.sez_cb=ttk.Combobox(p,textvariable=self.sez_var,
                                   values=list(SECTIONS.keys()),state="readonly",
                                   font=("Segoe UI",9))
        self.sez_cb.pack(fill="x",padx=12,pady=2)
        self.sez_info=tk.Label(p,text="",bg=PANEL,fg=TEXT2,font=("Segoe UI",8),
                                justify="left",anchor="w")
        self.sez_info.pack(fill="x",padx=12)
        self.sez_cb.bind("<<ComboboxSelected>>",self._upd_sez_info)
        self._upd_sez_info()

        # Svincoli
        sec_lbl(p,"SVINCOLI INTERNI (solo travi)")
        sf=tk.Frame(p,bg=PANEL); sf.pack(fill="x",padx=12)
        tk.Checkbutton(sf,text="Cerniera estremo i",variable=self.sv_i,
                       bg=PANEL,fg=TEXT,selectcolor=ACCENT,activebackground=PANEL,
                       font=("Segoe UI",9)).pack(anchor="w")
        tk.Checkbutton(sf,text="Cerniera estremo j",variable=self.sv_j,
                       bg=PANEL,fg=TEXT,selectcolor=ACCENT,activebackground=PANEL,
                       font=("Segoe UI",9)).pack(anchor="w")
        hsep(p)

        # Vincoli
        sec_lbl(p,"TIPO VINCOLO")
        vincoli_tipi = [
            ("⬛ Incastro",           "incastro"),
            ("◯ Cerniera",            "cerniera"),
            ("△ Carrello (vert.)",    "carrello"),
            ("▷ Carrello (oriz.)",    "carrello_x"),
            ("⟋ Carrello inclinato",  "carrello_inclinato"),
            ("⟂ Pattino (sciv. x)",   "pattino"),
            ("═ Incastro scorrevole", "incastro_scorrevole"),
            ("/ Biella",              "biella"),
            ("🌀 Molla",              "molla"),
        ]
        for txt,val in vincoli_tipi:
            tk.Radiobutton(p,text=txt,variable=self.vtype,value=val,
                           bg=PANEL,fg=TEXT,selectcolor=ACCENT,
                           activebackground=PANEL,font=("Segoe UI",9),
                           indicatoron=False,relief="flat",padx=6,pady=2,cursor="hand2"
                           ).pack(fill="x",padx=12,pady=1)
        vp=tk.Frame(p,bg=PANEL); vp.pack(fill="x",padx=12,pady=2)
        for r,(lbl,var) in enumerate([("Angolo biella/carr.°",self.v_ang),
                                       ("k_x (kN/m)",self.v_kx),
                                       ("k_y (kN/m)",self.v_ky),
                                       ("k_phi (kNm/rad)",self.v_kphi)]):
            tk.Label(vp,text=lbl,bg=PANEL,fg=TEXT2,font=("Segoe UI",8)
                     ).grid(row=r,column=0,sticky="w",pady=1)
            tk.Entry(vp,textvariable=var,bg=PANEL2,fg=TEXT,insertbackground=TEXT,
                     relief="flat",font=("Segoe UI",9),bd=3,width=8
                     ).grid(row=r,column=1,sticky="ew",padx=(4,0))
        vp.columnconfigure(1,weight=1)
        hsep(p)

        # Carico
        sec_lbl(p,"TIPO CARICO")
        for txt,val in [("↓ Fy nodale","Fy"),("→ Fx nodale","Fx"),
                        ("↺ Momento M","M"),("▬ Uniforme","uniforme"),
                        ("◤ Triang. sx","triang_sx"),("◥ Triang. dx","triang_dx")]:
            tk.Radiobutton(p,text=txt,variable=self.ctype,value=val,
                           bg=PANEL,fg=TEXT,selectcolor=ACCENT,
                           activebackground=PANEL,font=("Segoe UI",9),
                           indicatoron=False,relief="flat",padx=6,pady=2,cursor="hand2"
                           ).pack(fill="x",padx=12,pady=1)
        cp=tk.Frame(p,bg=PANEL); cp.pack(fill="x",padx=12,pady=2)
        for r,(lbl,var) in enumerate([("Valore (kN/kNm/kN/m)",self.cval),
                                       ("ΔT termico (°C)",self.cdT)]):
            tk.Label(cp,text=lbl,bg=PANEL,fg=TEXT2,font=("Segoe UI",8)
                     ).grid(row=r,column=0,sticky="w",pady=1)
            tk.Entry(cp,textvariable=var,bg=PANEL2,fg=TEXT,insertbackground=TEXT,
                     relief="flat",font=("Segoe UI",9),bd=3,width=8
                     ).grid(row=r,column=1,sticky="ew",padx=(4,0))
        cp.columnconfigure(1,weight=1)
        hsep(p)

        # Vista
        sec_lbl(p,"VISUALIZZAZIONE")
        for txt,var,cmd in [("Deformata",self.show_def,self._redraw),
                             ("Diagrammi M/V/N",self.show_mvn,self._toggle_mvn),
                             ("Mappa tensioni",self.show_heat,self._redraw),
                             ("Etichette",self.show_labels,self._redraw),
                             ("Snap griglia",self.snap_g,None),
                             ("Snap nodi",self.snap_n,None)]:
            tk.Checkbutton(p,text=txt,variable=var,command=cmd,
                           bg=PANEL,fg=TEXT,selectcolor=ACCENT,
                           activebackground=PANEL,font=("Segoe UI",9)
                           ).pack(anchor="w",padx=12,pady=1)
        hsep(p)

        # Modale
        sec_lbl(p,"ANALISI MODALE")
        mp=tk.Frame(p,bg=PANEL); mp.pack(fill="x",padx=12)
        tk.Label(mp,text="N° modi:",bg=PANEL,fg=TEXT2,font=("Segoe UI",9)).grid(row=0,column=0,sticky="w")
        tk.Spinbox(mp,from_=1,to=20,textvariable=self.n_modi_var,width=4,
                   bg=PANEL2,fg=TEXT,insertbackground=TEXT,relief="flat",
                   font=("Segoe UI",9)).grid(row=0,column=1,padx=4)
        tk.Label(mp,text="Modo:",bg=PANEL,fg=TEXT2,font=("Segoe UI",9)).grid(row=1,column=0,sticky="w",pady=2)
        tk.Spinbox(mp,from_=0,to=19,textvariable=self.modo_var,width=4,
                   bg=PANEL2,fg=TEXT,insertbackground=TEXT,relief="flat",
                   font=("Segoe UI",9),command=self._redraw).grid(row=1,column=1,padx=4)
        tk.Checkbutton(p,text="Mostra modo proprio",variable=self.show_modo,
                       command=self._redraw,bg=PANEL,fg=TEXT,selectcolor=ACCENT,
                       activebackground=PANEL,font=("Segoe UI",9)
                       ).pack(anchor="w",padx=12)
        hsep(p)

        # Libreria
        sec_lbl(p,"STRUTTURE PREDEFINITE")
        for nome in LIBRERIA:
            tk.Button(p,text=nome,command=lambda n=nome:self._load_preset(n),
                      bg=PANEL2,fg=TEXT,relief="flat",font=("Segoe UI",9),
                      padx=6,pady=2,cursor="hand2",
                      activebackground=ACCENT,activeforeground=BG
                      ).pack(fill="x",padx=12,pady=1)
        hsep(p)
        mk_btn(p,"▶ RISOLVI",self._solve).pack(fill="x",padx=12,pady=2)
        mk_btn(p,"🎵 MODALE",self._solve_modale,bg="#1a4a3a").pack(fill="x",padx=12,pady=2)
        mk_btn(p,"🗑 Cancella tutto",self._clear_all,bg=DANGER).pack(fill="x",padx=12,pady=2)

    def _upd_mat_info(self,*_):
        mat=MATERIALS.get(self.mat_var.get(),{})
        self.mat_info.config(text=f"E={mat.get('E',0)/1e6:.0f}GPa  fy={mat.get('fy',0)/1e3:.0f}MPa  ρ={mat.get('rho',0):.0f}kg/m³  α={mat.get('alpha',0)*1e6:.1f}μ/°C")

    def _upd_sez_info(self,*_):
        sez=SECTIONS.get(self.sez_var.get(),{})
        self.sez_info.config(text=f"A={sez.get('A',0)*1e4:.2f}cm²  I={sez.get('I',0)*1e8:.2f}cm⁴  W={sez.get('W',0)*1e6:.1f}cm³")

    def _build_center(self):
        cf=tk.Frame(self,bg=BG)
        cf.grid(row=1,column=1,sticky="nsew")
        cf.rowconfigure(0,weight=1); cf.columnconfigure(0,weight=1)
        self.canvas=tk.Canvas(cf,bg=BG,highlightthickness=0,cursor="crosshair")
        self.canvas.grid(row=0,column=0,sticky="nsew")
        self.mvn_frame=tk.Frame(cf,bg=PANEL2,height=230)
        self.mvn_frame.grid(row=1,column=0,sticky="ew"); self.mvn_frame.grid_remove()
        self.status=tk.StringVar(value="Pronto")
        tk.Label(cf,textvariable=self.status,bg=PANEL2,fg=TEXT2,
                 font=("Segoe UI",9),anchor="w",padx=10).grid(row=2,column=0,sticky="ew")
        self._bind_canvas()

    def _build_right(self):
        frame=tk.Frame(self,bg=PANEL,width=290)
        frame.grid(row=1,column=2,sticky="nsew"); frame.grid_propagate(False)
        tk.Label(frame,text="RISULTATI",bg=PANEL,fg=ACCENT,
                 font=("Segoe UI",11,"bold")).pack(anchor="w",padx=12,pady=(12,4))
        nb=ttk.Notebook(frame)
        nb.pack(fill="both",expand=True,padx=6,pady=4)
        style=ttk.Style()
        style.configure("TNotebook",background=PANEL,borderwidth=0)
        style.configure("TNotebook.Tab",background=PANEL2,foreground=TEXT2,
                        font=("Segoe UI",9),padding=[6,3])
        style.map("TNotebook.Tab",background=[("selected",ACCENT)],foreground=[("selected",BG)])
        self._rt={}
        for name,lbl in [("react","Reazioni"),("displ","Spostamenti"),("modal","Modale"),("elem","Elementi")]:
            t=tk.Frame(nb,bg=PANEL2); nb.add(t,text=lbl)
            txt=tk.Text(t,bg=PANEL2,fg=TEXT,font=("Consolas",8),relief="flat",state="disabled",wrap="none")
            sb=tk.Scrollbar(t,command=txt.yview,bg=PANEL2)
            txt.config(yscrollcommand=sb.set)
            sb.pack(side="right",fill="y"); txt.pack(fill="both",expand=True,padx=4,pady=4)
            self._rt[name]=txt
        self.verify_lbl=tk.Label(frame,text="—",bg=PANEL,fg=TEXT2,
                                  font=("Segoe UI",9),wraplength=270,justify="left")
        self.verify_lbl.pack(anchor="w",padx=12,pady=(0,8))

    # ── coordinate ───────────────────────────────────────────────────────────

    def w2s(self,wx,wy):
        W=self.canvas.winfo_width() or 900; H=self.canvas.winfo_height() or 600
        return W/2+self.cam_x+wx*self.cam_scale, H/2+self.cam_y-wy*self.cam_scale

    def s2w(self,sx,sy):
        W=self.canvas.winfo_width() or 900; H=self.canvas.winfo_height() or 600
        return (sx-W/2-self.cam_x)/self.cam_scale, -(sy-H/2-self.cam_y)/self.cam_scale

    def _snap(self,wx,wy,excl=None):
        if self.snap_n.get():
            bd=12/self.cam_scale
            for i,n in enumerate(self.nodi):
                if i==excl: continue
                d=math.hypot(wx-n['x'],wy-n['y'])
                if d<bd: return n['x'],n['y'],True
        if self.snap_g.get():
            g=SNAP_G; return round(wx/g)*g,round(wy/g)*g,False
        return wx,wy,False

    # ── canvas binding ───────────────────────────────────────────────────────

    def _bind_canvas(self):
        c=self.canvas
        c.bind("<Configure>",      lambda e:self._redraw())
        c.bind("<ButtonPress-1>",  self._on_click)
        c.bind("<B1-Motion>",      self._on_drag)
        c.bind("<ButtonRelease-1>",self._on_release)
        c.bind("<ButtonPress-3>",  self._pan_start)
        c.bind("<B3-Motion>",      self._pan_move)
        c.bind("<MouseWheel>",     self._on_wheel)
        c.bind("<Button-4>",       lambda e:self._zoom(1.1))
        c.bind("<Button-5>",       lambda e:self._zoom(0.9))
        c.bind("<Delete>",         self._del_sel)
        c.bind("<BackSpace>",      self._del_sel)
        c.focus_set()

    def _find_nodo(self,sx,sy,tol=14):
        best,bd=None,tol
        for i,n in enumerate(self.nodi):
            px,py=self.w2s(n['x'],n['y'])
            d=math.hypot(sx-px,sy-py)
            if d<bd: best,bd=i,d
        return best

    def _find_elem(self,sx,sy,tol=8):
        for i,el in enumerate(self.elementi):
            ni,nj=self.nodi[el['i']],self.nodi[el['j']]
            ax,ay=self.w2s(ni['x'],ni['y']); bx,by=self.w2s(nj['x'],nj['y'])
            dx,dy=bx-ax,by-ay; L2=dx*dx+dy*dy
            if L2<1: continue
            t=max(0,min(1,((sx-ax)*dx+(sy-ay)*dy)/L2))
            px,py=ax+t*dx,ay+t*dy
            if math.hypot(sx-px,sy-py)<tol: return i
        return None

    def _on_click(self,e):
        self.canvas.focus_set(); m=self.mode.get()
        wx,wy,_=self._snap(*self.s2w(e.x,e.y))
        if m=="nodo":
            ni=self._find_nodo(e.x,e.y)
            if ni is None:
                self._push_undo()
                self.nodi.append({'x':wx,'y':wy}); self.U=None
            else: self.sel_nodo=ni
            self._redraw()
        elif m=="trave":
            ni=self._find_nodo(e.x,e.y)
            if ni is not None: self.drag_from=ni; self.drag_pos=(e.x,e.y)
        elif m=="vincolo":
            ni=self._find_nodo(e.x,e.y)
            if ni is not None:
                self._push_undo(); self._apply_vincolo(ni); self.U=None; self._redraw()
        elif m=="carico":
            ct=self.ctype.get()
            try: val=float(self.cval.get()); dT=float(self.cdT.get())
            except: return messagebox.showerror("Errore","Valore non valido")
            if ct in ('Fy','Fx','M'):
                ni=self._find_nodo(e.x,e.y)
                if ni is not None:
                    self._push_undo()
                    self.carichi.append({'tipo':ct,'nodo':ni,'val':val,'dT':dT})
                    self.U=None; self._redraw()
            else:
                ei=self._find_elem(e.x,e.y)
                if ei is not None:
                    self._push_undo()
                    self.carichi.append({'tipo':ct,'elem':ei,'val':val,'dT':dT})
                    self.U=None; self._redraw()
        elif m=="sposta":
            ni=self._find_nodo(e.x,e.y)
            self.sel_nodo=ni; self.dragging=ni is not None; self._redraw()

    def _on_drag(self,e):
        m=self.mode.get()
        if m=="trave" and self.drag_from is not None:
            self.drag_pos=(e.x,e.y); self._redraw()
        elif m=="sposta" and self.dragging and self.sel_nodo is not None:
            wx,wy,_=self._snap(*self.s2w(e.x,e.y),excl=self.sel_nodo)
            self.nodi[self.sel_nodo]['x']=wx; self.nodi[self.sel_nodo]['y']=wy
            self.U=None; self._redraw()

    def _on_release(self,e):
        m=self.mode.get()
        if m=="trave" and self.drag_from is not None:
            nj=self._find_nodo(e.x,e.y)
            if nj is not None and nj!=self.drag_from:
                ex=any((t['i']==self.drag_from and t['j']==nj) or
                       (t['i']==nj and t['j']==self.drag_from)
                       for t in self.elementi)
                if not ex:
                    self._push_undo()
                    self.elementi.append({
                        'tipo': self.el_tipo.get(),
                        'i': self.drag_from, 'j': nj,
                        'mat': self.mat_var.get(),
                        'sez': self.sez_var.get(),
                        'svincolo_i': self.sv_i.get(),
                        'svincolo_j': self.sv_j.get(),
                    }); self.U=None
            self.drag_from=None; self.drag_pos=None; self._redraw()
        elif m=="sposta":
            if self.dragging: self._push_undo()
            self.dragging=False

    def _pan_start(self,e): self.pan_start=(e.x-self.cam_x,e.y-self.cam_y)
    def _pan_move(self,e):
        if self.pan_start:
            self.cam_x=e.x-self.pan_start[0]; self.cam_y=e.y-self.pan_start[1]; self._redraw()

    def _on_wheel(self,e):
        f=1.15 if e.delta>0 else 0.87
        W=self.canvas.winfo_width(); H=self.canvas.winfo_height()
        wx0,wy0=self.s2w(e.x,e.y)
        self.cam_scale=max(10,min(600,self.cam_scale*f))
        self.cam_x=e.x-W/2-wx0*self.cam_scale; self.cam_y=e.y-H/2+wy0*self.cam_scale
        self._redraw()

    def _zoom(self,f):
        W=self.canvas.winfo_width(); H=self.canvas.winfo_height()
        wx0,wy0=self.s2w(W/2,H/2)
        self.cam_scale=max(10,min(600,self.cam_scale*f))
        self.cam_x=W/2-W/2-wx0*self.cam_scale; self.cam_y=H/2-H/2+wy0*self.cam_scale
        self._redraw()

    def _reset_view(self):
        if not self.nodi: self.cam_x=0;self.cam_y=0;self.cam_scale=80;self._redraw();return
        xs=[n['x'] for n in self.nodi]; ys=[n['y'] for n in self.nodi]
        W=self.canvas.winfo_width() or 900; H=self.canvas.winfo_height() or 550
        spanX=max(max(xs)-min(xs),1); spanY=max(max(ys)-min(ys),1)
        self.cam_scale=min(W/(spanX*3),H/(spanY*3),200); self.cam_scale=max(self.cam_scale,20)
        cx=(max(xs)+min(xs))/2; cy=(max(ys)+min(ys))/2
        self.cam_x=-cx*self.cam_scale; self.cam_y=cy*self.cam_scale
        self._redraw()

    # ── vincoli ──────────────────────────────────────────────────────────────

    def _apply_vincolo(self,ni):
        self.vincoli=[v for v in self.vincoli if v['nodo']!=ni]
        try:
            ang=float(self.v_ang.get())
            kx=float(self.v_kx.get())
            ky=float(self.v_ky.get())
            kphi=float(self.v_kphi.get())
        except: ang=kx=ky=kphi=0
        self.vincoli.append({'nodo':ni,'tipo':self.vtype.get(),
                              'angolo':math.radians(ang),
                              'kx':kx,'ky':ky,'kphi':kphi})

    # ── solve ─────────────────────────────────────────────────────────────────

    def _solve(self):
        if len(self.nodi)<2 or len(self.elementi)<1 or len(self.vincoli)<1:
            return messagebox.showwarning("Attenzione","Servono almeno 2 nodi, 1 elemento e 1 vincolo")
        try:
            K,M,F,is_fixed=assembla_sistema(self.nodi,self.elementi,self.vincoli,self.carichi)
            U,R,free,fixed=solve_statico(K,F,is_fixed)
            self.U=U; self.R=R; self.fixed=fixed; self.K=K; self.F=F; self.M_mat=M
            self._show_results_static()
            self._show_results_elem()
            self._redraw()
            if self.show_mvn.get(): self._draw_mvn()
            self.status.set("✓ Analisi statica completata")
        except Exception as ex:
            messagebox.showerror("Errore solver",str(ex)); self.U=None

    def _solve_modale(self):
        if len(self.nodi)<2 or len(self.elementi)<1 or len(self.vincoli)<1:
            return messagebox.showwarning("Attenzione","Prima imposta la struttura completa")
        try:
            K,M,F,is_fixed=assembla_sistema(self.nodi,self.elementi,self.vincoli,self.carichi)
            n_m=self.n_modi_var.get()
            freq,omega,modi=solve_modale(K,M,is_fixed,n_m)
            self.freq=freq; self.modi=modi
            self._show_results_modal(freq,omega)
            self.status.set(f"✓ Analisi modale: {len(freq)} modi calcolati")
        except Exception as ex:
            messagebox.showerror("Errore modale",str(ex))

    def _show_results_static(self):
        U,R,fixed=self.U,self.R,self.fixed
        comp=['Rx','Ry','Mz']; unit=['kN','kN','kN·m']
        lines=["═"*34,"  REAZIONI VINCOLARI","═"*34]
        for vi in fixed:
            nd=vi//3; k=vi%3
            lines.append(f"  {comp[k]:3s} @ N{nd+1:2d} = {R[vi]:+12.4f}  {unit[k]}")
        KU=self.K@U
        sfx=sum(KU[i] for i in range(len(U)) if i%3==0)
        sfy=sum(KU[i] for i in range(len(U)) if i%3==1)
        ok=abs(sfx)<1e-4 and abs(sfy)<1e-4
        lines+=["","─"*34,f"  ΣFx = {sfx:+.4e} kN",f"  ΣFy = {sfy:+.4e} kN",
                f"  {'✓ Equilibrio OK' if ok else '✗ ERRORE equilibrio'}","═"*34]
        self._set_txt(self._rt['react'],"\n".join(lines))
        self.verify_lbl.config(text="✓ Equilibrio OK" if ok else "✗ Errore",
                                fg=ACC2 if ok else DANGER)
        lines2=["═"*34,"  SPOSTAMENTI","═"*34]
        for i in range(len(self.nodi)):
            lines2.append(f"  N{i+1}: ux={U[i*3]:+.4e}m  uy={U[i*3+1]:+.4e}m  φ={U[i*3+2]:+.4e}rad")
        self._set_txt(self._rt['displ'],"\n".join(lines2))

    def _show_results_elem(self):
        if self.U is None: return
        lines=["═"*38,"  AZIONI INTERNE (max per elemento)","═"*38]
        for idx,el in enumerate(self.elementi):
            if el['tipo']=='truss':
                N=truss_force(el,self.nodi,self.U)
                stato="TRAZIONE" if N>0 else "COMPRESSIONE"
                lines.append(f"  T{idx+1} (Asta): N={N:+.4f} kN  [{stato}]")
            else:
                res=internal_forces(el,self.nodi,self.U)
                if res:
                    _,M,V,N=res
                    _,_,E,A,mat,sez=get_EI_EA(el,SECTIONS,MATERIALS)
                    fy=mat['fy']; W=sez['W']
                    sig_max=abs(M).max()/W if W>0 else 0
                    eta=sig_max/fy if fy>0 else 0
                    lines.append(f"  T{idx+1}: |M|max={abs(M).max():.3f}kNm  |V|max={abs(V).max():.3f}kN  |N|max={abs(N).max():.3f}kN")
                    lines.append(f"       σmax={sig_max/1e3:.2f}MPa  η={eta:.3f}  {'✓OK' if eta<1 else '✗SNERVAMENTO'}")
        self._set_txt(self._rt['elem'],"\n".join(lines))

    def _show_results_modal(self,freq,omega):
        lines=["═"*36,"  FREQUENZE PROPRIE","═"*36]
        for i,(f,w) in enumerate(zip(freq,omega)):
            T=1/f if f>1e-9 else float('inf')
            lines.append(f"  Modo {i+1:2d}: f={f:8.4f} Hz  T={T:8.4f} s  ω={w:8.3f} rad/s")
        lines+=["","═"*36]
        self._set_txt(self._rt['modal'],"\n".join(lines))

    def _set_txt(self,w,txt):
        w.config(state="normal"); w.delete("1.0","end")
        w.insert("1.0",txt); w.config(state="disabled")

    # ── disegno ──────────────────────────────────────────────────────────────

    def _redraw(self,*_):
        c=self.canvas; c.delete("all")
        W=c.winfo_width() or 900; H=c.winfo_height() or 600
        self._draw_grid(W,H)
        self._draw_elementi()
        if self.drag_from is not None and self.drag_pos:
            self._draw_preview()
        if self.show_def.get() and self.U is not None:
            self._draw_deformed(1.0)
        if self.show_modo.get() and self.modi is not None:
            idx=min(self.modo_var.get(),len(self.modi)-1)
            if idx>=0: self._draw_deformed_mode(self.modi[idx])
        if self.show_heat.get() and self.U is not None:
            self._draw_heatmap()
        self._draw_carichi()
        self._draw_vincoli()
        if self.U is not None: self._draw_reazioni()
        self._draw_nodi()

    def _draw_grid(self,W,H):
        c=self.canvas
        s=self._grid_step()
        tl=self.s2w(0,0); br=self.s2w(W,H)
        x0=math.floor(min(tl[0],br[0])/s)*s; x1=math.ceil(max(tl[0],br[0])/s)*s
        y0=math.floor(min(tl[1],br[1])/s)*s; y1=math.ceil(max(tl[1],br[1])/s)*s
        ox,oy=self.w2s(0,0)
        def fr(a,b,st):
            out=[]; x=a
            while x<=b+1e-9: out.append(round(x/st)*st); x+=st
            return out
        for gx in fr(x0,x1,s):
            sx,_=self.w2s(gx,0)
            col="#444466" if abs(gx)<1e-9 else GRID
            c.create_line(sx,0,sx,H,fill=col,width=1 if abs(gx)<1e-9 else 0.5)
            if abs(gx)>1e-9 and self.cam_scale>25:
                c.create_text(sx,min(oy+11,H-6),text=f"{gx:.1f}",fill=TEXT3,font=("Segoe UI",7))
        for gy in fr(y0,y1,s):
            _,sy=self.w2s(0,gy)
            col="#444466" if abs(gy)<1e-9 else GRID
            c.create_line(0,sy,W,sy,fill=col,width=1 if abs(gy)<1e-9 else 0.5)
            if abs(gy)>1e-9 and self.cam_scale>25:
                c.create_text(max(ox-18,20),sy,text=f"{gy:.1f}",fill=TEXT3,font=("Segoe UI",7))

    def _grid_step(self):
        mpp=1/self.cam_scale; raw=mpp*55
        if raw<=0: return 1
        exp=10**math.floor(math.log10(raw)); n=raw/exp
        if n<2: return exp
        if n<5: return 2*exp
        return 5*exp

    def _draw_elementi(self):
        c=self.canvas
        for idx,el in enumerate(self.elementi):
            ni,nj=self.nodi[el['i']],self.nodi[el['j']]
            ax,ay=self.w2s(ni['x'],ni['y']); bx,by=self.w2s(nj['x'],nj['y'])
            sel=self.sel_elem==idx
            is_truss=el['tipo']=='truss'
            col=SEL if sel else (TRUSS if is_truss else BEAM)
            w=3+sel
            if is_truss:
                c.create_line(ax,ay,bx,by,fill=col,width=w,dash=(8,4))
            else:
                c.create_line(ax,ay,bx,by,fill=col,width=w,capstyle="round")
            # svincoli
            if not is_truss and self.show_labels.get():
                if el.get('svincolo_i'):
                    c.create_oval(ax-5,ay-5,ax+5,ay+5,fill=BG,outline=ACC3,width=2)
                if el.get('svincolo_j'):
                    c.create_oval(bx-5,by-5,bx+5,by+5,fill=BG,outline=ACC3,width=2)
            if self.show_labels.get() and self.cam_scale>35:
                mx,my=(ax+bx)/2,(ay+by)/2
                L,_=elem_geo(ni,nj)
                lbl=f"{'A' if is_truss else 'T'}{idx+1} {L:.2f}m"
                c.create_text(mx,my-10,text=lbl,fill=col,font=("Segoe UI",7))

    def _draw_preview(self):
        ni=self.nodi[self.drag_from]
        ax,ay=self.w2s(ni['x'],ni['y']); bx,by=self.drag_pos
        self.canvas.create_line(ax,ay,bx,by,fill=ACCENT,width=2,dash=(6,4))

    def _draw_nodi(self):
        c=self.canvas; R=7
        for i,n in enumerate(self.nodi):
            sx,sy=self.w2s(n['x'],n['y'])
            col=SEL if i==self.sel_nodo else NODE
            c.create_oval(sx-R,sy-R,sx+R,sy+R,fill=BG,outline=col,width=2)
            if self.show_labels.get() and self.cam_scale>25:
                c.create_text(sx+12,sy-10,text=f"N{i+1}",fill=col,font=("Segoe UI",8,"bold"))

    def _draw_vincoli(self):
        c=self.canvas
        VINCOLO_COLORS={
            'incastro':ACC2,'cerniera':ACC2,'carrello':ACC2,'carrello_x':ACC2,
            'pattino':"#89dceb",'incastro_scorrevole':"#cba6f7",
            'biella':"#f38ba8",'molla':SPRING,'carrello_inclinato':ACC2,
        }
        for v in self.vincoli:
            n=self.nodi[v['nodo']]; sx,sy=self.w2s(n['x'],n['y'])
            tp=v['tipo']; col=VINCOLO_COLORS.get(tp,ACC2)
            sz=13
            if tp=='incastro':
                c.create_rectangle(sx-sz,sy-sz,sx+sz,sy+sz,fill="",outline=col,width=2)
                for off in range(-sz,sz+1,5):
                    c.create_line(sx+off,sy+sz,sx+off-4,sy+sz+7,fill=col,width=1)
            elif tp in ('cerniera','carrello','carrello_x','carrello_inclinato'):
                ang=v.get('angolo',0)
                if tp=='carrello': ang=0
                elif tp=='carrello_x': ang=math.pi/2
                nx_,ny_=math.sin(ang),math.cos(ang)
                # triangolo orientato
                pts=[sx,sy,
                     sx-sz*math.cos(ang-math.pi/6)*math.cos(0)+sz*nx_,
                     sy+sz*math.sin(ang-math.pi/6)*math.cos(0)+sz*ny_,
                     sx+sz*math.cos(ang+math.pi/6)*math.cos(0)+sz*nx_,
                     sy-sz*math.sin(ang+math.pi/6)*math.cos(0)+sz*ny_]
                # semplificato: triangolo verso il basso (locale y)
                pts=[sx,sy, sx-sz,sy+sz*1.6, sx+sz,sy+sz*1.6]
                c.create_polygon(pts,fill="",outline=col,width=2)
                if tp=='cerniera':
                    c.create_oval(sx-4,sy-4,sx+4,sy+4,fill=col)
                else:
                    c.create_oval(sx-5,sy+sz*1.6-5,sx+5,sy+sz*1.6+5,fill="",outline=col,width=2)
                    c.create_line(sx-sz-4,sy+sz*1.6+9,sx+sz+4,sy+sz*1.6+9,fill=col,width=2)
            elif tp=='pattino':
                c.create_rectangle(sx-sz,sy-sz//2,sx-1,sy+sz//2,fill="",outline=col,width=2)
                for off in [-sz//2,0,sz//2]:
                    c.create_line(sx-sz-1,sy+off,sx-sz-7,sy+off,fill=col,width=1)
            elif tp=='incastro_scorrevole':
                c.create_rectangle(sx-sz,sy-sz,sx+sz,sy+sz,fill="",outline=col,width=2)
                c.create_line(sx-sz,sy,sx+sz,sy,fill=col,width=1,dash=(3,2))
            elif tp=='biella':
                ang=v.get('angolo',0)
                ex=sx+sz*2*math.cos(ang); ey=sy-sz*2*math.sin(ang)
                c.create_line(sx,sy,ex,ey,fill=col,width=2)
                c.create_oval(sx-3,sy-3,sx+3,sy+3,fill=col)
                c.create_oval(ex-4,ey-4,ex+4,ey+4,fill="",outline=col,width=2)
            elif tp=='molla':
                kx=v.get('kx',0); ky=v.get('ky',0)
                if ky>0:
                    for k2 in range(5):
                        yy=sy+sz/2+k2*4
                        xoff=4 if k2%2==0 else -4
                        c.create_line(sx,yy,sx+xoff,yy+2,fill=col,width=1)
                    c.create_line(sx,sy+sz/2+20,sx,sy+sz/2+26,fill=col,width=2)
                if kx>0:
                    for k2 in range(5):
                        xx=sx-sz/2-k2*4
                        yoff=4 if k2%2==0 else -4
                        c.create_line(xx,sy,xx-2,sy+yoff,fill=col,width=1)
            c.create_text(sx+sz+4,sy-sz-2,text=tp[:3].upper(),fill=col,font=("Segoe UI",7))

    def _draw_carichi(self):
        c=self.canvas; aL=max(22,self.cam_scale*0.32)
        for load in self.carichi:
            tp=load['tipo']; val=load['val']
            if tp in ('Fy','Fx','M'):
                ni_idx=load.get('nodo')
                if ni_idx is None or ni_idx>=len(self.nodi): continue
                n=self.nodi[ni_idx]; sx,sy=self.w2s(n['x'],n['y'])
                if tp=='Fy':
                    dy=-aL if val>0 else aL
                    c.create_line(sx,sy+dy,sx,sy,fill=LOAD,width=2,arrow="last",arrowshape=(8,10,3))
                    c.create_text(sx+8,sy+dy/2,text=f"{abs(val):.1f}kN",fill=LOAD,font=("Segoe UI",8))
                elif tp=='Fx':
                    dx=aL if val>0 else -aL
                    c.create_line(sx-dx,sy,sx,sy,fill=LOAD,width=2,arrow="last",arrowshape=(8,10,3))
                    c.create_text(sx-dx/2,sy-12,text=f"{abs(val):.1f}kN",fill=LOAD,font=("Segoe UI",8))
                else:
                    c.create_text(sx,sy-22,text=f"↺{abs(val):.1f}kNm",fill=LOAD,font=("Segoe UI",9,"bold"))
            else:
                ei=load.get('elem')
                if ei is None or ei>=len(self.elementi): continue
                el=self.elementi[ei]
                if el['tipo']=='truss': continue
                ni,nj=self.nodi[el['i']],self.nodi[el['j']]
                ax,ay=self.w2s(ni['x'],ni['y']); bx,by=self.w2s(nj['x'],nj['y'])
                dx,dy=bx-ax,by-ay; L=math.hypot(dx,dy)
                if L<1: continue
                nx_,ny_=-dy/L,dx/L; q=abs(val); sc=min(26,aL*0.5)
                for k in range(9):
                    xi=k/8
                    qk=(q*(1-xi) if tp=='triang_sx' else q*xi if tp=='triang_dx' else q)
                    if q<1e-9: continue
                    px,py=ax+xi*dx,ay+xi*dy
                    tx,ty=px+nx_*qk/q*sc,py+ny_*qk/q*sc
                    c.create_line(tx,ty,px,py,fill=WARN,width=1,arrow="last",arrowshape=(5,6,2))
                mx,my=(ax+bx)/2,(ay+by)/2
                c.create_text(mx+nx_*sc,my+ny_*sc,text=f"{q:.1f}kN/m",fill=WARN,font=("Segoe UI",8))
                # carico termico
                dT=load.get('dT',0)
                if abs(dT)>1e-9:
                    c.create_text(mx,my+ny_*sc+12,text=f"ΔT={dT:+.1f}°C",fill="#fab387",font=("Segoe UI",8))

    def _draw_reazioni(self):
        c=self.canvas; aL=max(28,self.cam_scale*0.45)
        comp=['Rx','Ry','Mz']
        for vi in self.fixed:
            nd=vi//3; k=vi%3; val=self.R[vi]
            if abs(val)<1e-6 or nd>=len(self.nodi): continue
            n=self.nodi[nd]; sx,sy=self.w2s(n['x'],n['y'])
            if k==1:
                dy=aL if val>0 else -aL
                c.create_line(sx,sy,sx,sy-dy,fill=REACT,width=2,arrow="last",arrowshape=(10,12,4))
                c.create_text(sx+6,sy-dy/2,text=f"Ry={val:+.2f}kN",fill=REACT,font=("Consolas",8))
            elif k==0:
                dx=aL if val>0 else -aL
                c.create_line(sx,sy,sx+dx,sy,fill=REACT,width=2,arrow="last",arrowshape=(10,12,4))
                c.create_text(sx+dx/2,sy-12,text=f"Rx={val:+.2f}kN",fill=REACT,font=("Consolas",8))
            else:
                c.create_text(sx+32,sy+18,text=f"Mz={val:+.2f}kNm",fill=REACT,font=("Consolas",8))

    def _draw_deformed(self,amp_factor=1.0):
        if self.U is None: return
        U=self.U
        absU=np.abs(U); maxU=float(absU[absU>0].max()) if (absU>0).any() else 1.0
        if not self.nodi: return
        xs=[n['x'] for n in self.nodi]; ys=[n['y'] for n in self.nodi]
        span=max(max(xs)-min(xs),max(ys)-min(ys),0.1)
        amp=span*0.08/maxU*float(self.def_scale.get())*amp_factor
        phase=math.sin(self.anim_ph) if self.anim_run else 1.0
        c=self.canvas
        for el in self.elementi:
            if el['tipo']=='truss': continue
            ni,nj=self.nodi[el['i']],self.nodi[el['j']]
            L,theta=elem_geo(ni,nj)
            if L<1e-9: continue
            cs,sn=math.cos(theta),math.sin(theta)
            i,j=el['i'],el['j']
            u=[U[i*3]*amp*phase,U[i*3+1]*amp*phase,U[i*3+2]*amp*phase,
               U[j*3]*amp*phase,U[j*3+1]*amp*phase,U[j*3+2]*amp*phase]
            steps=28; pts=[]
            for k in range(steps+1):
                xi=k/steps
                N1=1-3*xi**2+2*xi**3; N2=xi*(1-xi)**2*L
                N3=3*xi**2-2*xi**3;   N4=xi**2*(xi-1)*L
                ul_i=cs*u[0]+sn*u[1]; vl_i=-sn*u[0]+cs*u[1]
                ul_j=cs*u[3]+sn*u[4]; vl_j=-sn*u[3]+cs*u[4]
                ul=ul_i*(1-xi)+ul_j*xi
                vl=N1*vl_i+N2*u[2]+N3*vl_j+N4*u[5]
                xg=ni['x']+(nj['x']-ni['x'])*xi+cs*ul-sn*vl
                yg=ni['y']+(nj['y']-ni['y'])*xi+sn*ul+cs*vl
                sx,sy=self.w2s(xg,yg); pts.extend([sx,sy])
            if len(pts)>=4:
                c.create_line(pts,fill=DEFORM,width=1.5,dash=(5,3),smooth=True)

    def _draw_deformed_mode(self,phi):
        absP=np.abs(phi); maxP=float(absP[absP>0].max()) if (absP>0).any() else 1.0
        if not self.nodi: return
        xs=[n['x'] for n in self.nodi]; ys=[n['y'] for n in self.nodi]
        span=max(max(xs)-min(xs),max(ys)-min(ys),0.1)
        amp=span*0.12/maxP*float(self.def_scale.get())
        phase=math.sin(self.anim_ph) if self.anim_run else 1.0
        c=self.canvas
        for el in self.elementi:
            if el['tipo']=='truss': continue
            ni,nj=self.nodi[el['i']],self.nodi[el['j']]
            L,theta=elem_geo(ni,nj)
            if L<1e-9: continue
            cs,sn=math.cos(theta),math.sin(theta)
            i,j=el['i'],el['j']
            u=[phi[i*3]*amp*phase,phi[i*3+1]*amp*phase,phi[i*3+2]*amp*phase,
               phi[j*3]*amp*phase,phi[j*3+1]*amp*phase,phi[j*3+2]*amp*phase]
            steps=28; pts=[]
            for k in range(steps+1):
                xi=k/steps
                N1=1-3*xi**2+2*xi**3; N2=xi*(1-xi)**2*L
                N3=3*xi**2-2*xi**3;   N4=xi**2*(xi-1)*L
                ul_i=cs*u[0]+sn*u[1]; vl_i=-sn*u[0]+cs*u[1]
                ul_j=cs*u[3]+sn*u[4]; vl_j=-sn*u[3]+cs*u[4]
                ul=ul_i*(1-xi)+ul_j*xi
                vl=N1*vl_i+N2*u[2]+N3*vl_j+N4*u[5]
                xg=ni['x']+(nj['x']-ni['x'])*xi+cs*ul-sn*vl
                yg=ni['y']+(nj['y']-ni['y'])*xi+sn*ul+cs*vl
                sx,sy=self.w2s(xg,yg); pts.extend([sx,sy])
            if len(pts)>=4:
                c.create_line(pts,fill=MOMENT,width=2,dash=(4,3),smooth=True)

    def _draw_heatmap(self):
        if self.U is None: return
        c=self.canvas
        for el in self.elementi:
            if el['tipo']=='truss':
                N=truss_force(el,self.nodi,self.U)
                ni,nj=self.nodi[el['i']],self.nodi[el['j']]
                ax,ay=self.w2s(ni['x'],ni['y']); bx,by=self.w2s(nj['x'],nj['y'])
                col="#e46d6d" if N<0 else "#a6e3a1"
                c.create_line(ax,ay,bx,by,fill=col,width=5,capstyle="round")
                mx,my=(ax+bx)/2,(ay+by)/2
                c.create_text(mx,my-10,text=f"N={N:+.1f}kN",fill=col,font=("Consolas",8))
                continue
            res=internal_forces(el,self.nodi,self.U)
            if res is None: continue
            _,M,V,N=res; Mabs=np.abs(M)
            sig_max=Mabs.max()
            if sig_max<1e-9: continue
            ni,nj=self.nodi[el['i']],self.nodi[el['j']]
            L,theta=elem_geo(ni,nj)
            cs,sn=math.cos(theta),math.sin(theta)
            steps=len(M)-1
            xs_=np.linspace(0,L,steps+1)
            for k in range(steps):
                xi=xs_[k]/L; ratio=Mabs[k]/sig_max
                if ratio>0.66: r,g=255,int((1-ratio)*3*255)
                elif ratio>0.33: r,g=int((ratio-0.33)*3*255),255
                else: r,g=0,int(ratio*3*255)
                col=f"#{min(r,255):02x}{min(g,255):02x}00"
                xi2=xs_[k+1]/L
                ax2=ni['x']+xi*(nj['x']-ni['x']); ay2=ni['y']+xi*(nj['y']-ni['y'])
                bx2=ni['x']+xi2*(nj['x']-ni['x']); by2=ni['y']+xi2*(nj['y']-ni['y'])
                px1,py1=self.w2s(ax2,ay2); px2,py2=self.w2s(bx2,by2)
                c.create_line(px1,py1,px2,py2,fill=col,width=5,capstyle="round")

    # ── diagrammi M/V/N ──────────────────────────────────────────────────────

    def _toggle_def(self):
        self.show_def.set(not self.show_def.get()); self._redraw()

    def _toggle_mvn(self):
        if self.show_mvn.get():
            self.mvn_frame.grid()
            if self.U is not None: self._draw_mvn()
        else:
            self.mvn_frame.grid_remove()

    def _toggle_heat(self):
        self.show_heat.set(not self.show_heat.get()); self._redraw()

    def _draw_mvn(self):
        if not HAS_MPL:
            messagebox.showwarning("matplotlib","pip install matplotlib")
            self.show_mvn.set(False); return
        if self.U is None: return
        for w in self.mvn_frame.winfo_children(): w.destroy()
        fig,axes=plt.subplots(1,3,figsize=(13,2.4),facecolor=PANEL2)
        plt.subplots_adjust(left=0.04,right=0.98,top=0.82,bottom=0.22,wspace=0.35)
        titles=["Momento M [kNm]","Taglio V [kN]","Normale N [kN]"]
        cols=[MOMENT,SHEAR,NORMAL]
        x_off=0.0; xticks=[0.0]; xlbls=["0"]
        for el in self.elementi:
            if el['tipo']=='truss': continue
            res=internal_forces(el,self.nodi,self.U)
            if res is None: continue
            xs,M,V,N=res; L=xs[-1]; xp=xs+x_off
            for ax,vals,col in zip(axes,[M,V,N],cols):
                ax.fill_between(xp,vals,0,alpha=0.25,color=col)
                ax.plot(xp,vals,color=col,linewidth=1.5)
                ax.axvline(x_off+L,color=BORDER,linewidth=0.5,linestyle=":")
            x_off+=L; xticks.append(x_off); xlbls.append(f"{x_off:.1f}")
        for ax,t,col in zip(axes,titles,cols):
            ax.set_facecolor(PANEL2); ax.tick_params(colors=TEXT2,labelsize=7)
            for sp in ax.spines.values(): sp.set_color(BORDER)
            ax.set_title(t,color=TEXT2,fontsize=8,pad=3)
            ax.axhline(0,color=BORDER,linewidth=0.7)
            ax.set_xticks(xticks); ax.set_xticklabels(xlbls,fontsize=7,color=TEXT2)
            ax.grid(True,alpha=0.12,linestyle="--",color=BORDER)
        fc=FigureCanvasTkAgg(fig,master=self.mvn_frame)
        fc.draw(); fc.get_tk_widget().pack(fill="both",expand=True); plt.close(fig)

    # ── animazione ────────────────────────────────────────────────────────────

    def _toggle_anim(self):
        if self.U is None and self.modi is None:
            return messagebox.showwarning("Attenzione","Esegui prima il calcolo")
        if self.anim_run:
            self.anim_run=False; self.anim_btn.config(text="🎬 Anima")
        else:
            self.anim_run=True; self.show_def.set(True)
            self.anim_btn.config(text="⏹ Stop"); self._anim_loop()

    def _anim_loop(self):
        if not self.anim_run: return
        self.anim_ph+=0.07; self._redraw(); self.after(28,self._anim_loop)

    # ── finestra risultati avanzata ───────────────────────────────────────────

    def _open_results_window(self):
        if self.U is None:
            return messagebox.showwarning("Attenzione","Esegui prima l'analisi statica")
        win=tk.Toplevel(self); win.title("Risultati Avanzati — FEM Solver Pro v2")
        win.configure(bg=BG); win.geometry("1100x750")

        nb=ttk.Notebook(win); nb.pack(fill="both",expand=True,padx=8,pady=8)
        style=ttk.Style(); style.configure("TNotebook",background=BG)
        style.configure("TNotebook.Tab",background=PANEL2,foreground=TEXT2,
                        font=("Segoe UI",10),padding=[10,5])
        style.map("TNotebook.Tab",background=[("selected",ACCENT)],foreground=[("selected",BG)])

        # Tab 1: Verifica resistenza
        t1=tk.Frame(nb,bg=BG); nb.add(t1,text="Verifica Resistenza")
        self._build_verify_tab(t1)

        # Tab 2: Diagrammi interattivi
        t2=tk.Frame(nb,bg=BG); nb.add(t2,text="Diagrammi M/V/N")
        self._build_mvn_tab(t2)

        # Tab 3: Analisi modale
        t3=tk.Frame(nb,bg=BG); nb.add(t3,text="Analisi Modale")
        self._build_modal_tab(t3)

        # Tab 4: Spostamenti nodali
        t4=tk.Frame(nb,bg=BG); nb.add(t4,text="Spostamenti")
        self._build_displ_tab(t4)

        # Tab 5: Tabella elementi
        t5=tk.Frame(nb,bg=BG); nb.add(t5,text="Tabella Completa")
        self._build_full_table_tab(t5)

    def _build_verify_tab(self,parent):
        if not HAS_MPL: tk.Label(parent,text="Installa matplotlib",bg=BG,fg=DANGER).pack(); return
        fig,axes=plt.subplots(1,2,figsize=(11,5),facecolor=BG)
        plt.subplots_adjust(left=0.08,right=0.96,top=0.88,bottom=0.12,wspace=0.35)

        nomi=[]; eta_vals=[]; M_maxs=[]; fy_vals=[]
        for idx,el in enumerate(self.elementi):
            if el['tipo']=='truss':
                N=truss_force(el,self.nodi,self.U)
                _,EA,E,A,mat,sez=get_EI_EA(el,SECTIONS,MATERIALS)
                fy=mat['fy']; sig_N=abs(N)/A if A>0 else 0
                eta=sig_N/fy if fy>0 else 0
                nomi.append(f"A{idx+1}"); eta_vals.append(eta); M_maxs.append(abs(N)); fy_vals.append(fy*A/1000)
            else:
                res=internal_forces(el,self.nodi,self.U)
                if res is None: continue
                _,M,V,N=res
                _,_,E,A,mat,sez=get_EI_EA(el,SECTIONS,MATERIALS)
                fy=mat['fy']; W=sez['W']
                sig=abs(M).max()/W if W>0 else 0
                eta=sig/fy if fy>0 else 0
                nomi.append(f"T{idx+1}"); eta_vals.append(eta); M_maxs.append(abs(M).max()); fy_vals.append(fy*W/1000)

        # bar chart coefficienti
        ax1=axes[0]; ax1.set_facecolor(PANEL2)
        colors_=[DANGER if e>1 else WARN if e>0.8 else NORMAL for e in eta_vals]
        bars=ax1.bar(nomi,eta_vals,color=colors_,edgecolor=BORDER,linewidth=0.5)
        ax1.axhline(1.0,color=DANGER,linewidth=1.5,linestyle="--",label="Limite η=1")
        ax1.axhline(0.8,color=WARN,linewidth=1,linestyle=":",alpha=0.7,label="Attenzione η=0.8")
        ax1.set_title("Coefficiente di sfruttamento η = σ/fy",color=TEXT,fontsize=11,pad=8)
        ax1.set_ylabel("η [ ]",color=TEXT2); ax1.tick_params(colors=TEXT2)
        for sp in ax1.spines.values(): sp.set_color(BORDER)
        ax1.legend(facecolor=PANEL,edgecolor=BORDER,labelcolor=TEXT2,fontsize=9)
        ax1.set_ylim(0,max(max(eta_vals,default=0)*1.2,1.2))
        for bar,eta in zip(bars,eta_vals):
            ax1.text(bar.get_x()+bar.get_width()/2,bar.get_height()+0.02,
                     f"{eta:.3f}",ha='center',fontsize=8,color=TEXT2)

        # M vs Mr
        ax2=axes[1]; ax2.set_facecolor(PANEL2)
        x_pos=range(len(nomi))
        ax2.bar([x-0.2 for x in x_pos],M_maxs,width=0.35,label="|M|max / |N|max act.",color=MOMENT,alpha=0.8)
        ax2.bar([x+0.2 for x in x_pos],fy_vals,width=0.35,label="Resistenza Mr [kNm/kN]",color=NORMAL,alpha=0.8)
        ax2.set_xticks(list(x_pos)); ax2.set_xticklabels(nomi)
        ax2.set_title("Azione vs Resistenza",color=TEXT,fontsize=11,pad=8)
        ax2.set_ylabel("kNm / kN",color=TEXT2); ax2.tick_params(colors=TEXT2)
        for sp in ax2.spines.values(): sp.set_color(BORDER)
        ax2.legend(facecolor=PANEL,edgecolor=BORDER,labelcolor=TEXT2,fontsize=9)

        fc=FigureCanvasTkAgg(fig,master=parent)
        fc.draw(); fc.get_tk_widget().pack(fill="both",expand=True); plt.close(fig)

    def _build_mvn_tab(self,parent):
        if not HAS_MPL or self.U is None:
            tk.Label(parent,text="Nessun risultato disponibile",bg=BG,fg=TEXT2).pack(pady=20); return
        fig,axes=plt.subplots(3,1,figsize=(11,6),facecolor=BG,sharex=True)
        plt.subplots_adjust(left=0.08,right=0.97,top=0.93,bottom=0.1,hspace=0.3)
        titles=["Momento M [kNm]","Taglio V [kN]","Sforzo normale N [kN]"]
        cols=[MOMENT,SHEAR,NORMAL]
        x_off=0.0; xticks=[0.0]; xlbls=["0"]
        for el in self.elementi:
            if el['tipo']=='truss': continue
            res=internal_forces(el,self.nodi,self.U)
            if res is None: continue
            xs,M,V,N=res; L=xs[-1]; xp=xs+x_off
            for ax,vals,col in zip(axes,[M,V,N],cols):
                ax.fill_between(xp,vals,0,alpha=0.3,color=col)
                ax.plot(xp,vals,color=col,linewidth=2)
                ax.axvline(x_off+L,color=BORDER,linewidth=0.7,linestyle=":")
            x_off+=L; xticks.append(x_off); xlbls.append(f"{x_off:.2f}")
        for ax,t,col in zip(axes,titles,cols):
            ax.set_facecolor(PANEL2); ax.tick_params(colors=TEXT2)
            for sp in ax.spines.values(): sp.set_color(BORDER)
            ax.set_title(t,color=TEXT2,fontsize=10,pad=4)
            ax.axhline(0,color=BORDER,linewidth=0.8)
            ax.grid(True,alpha=0.15,linestyle="--",color=BORDER)
            ax.set_ylabel(t.split()[-1],color=TEXT2)
        axes[-1].set_xticks(xticks); axes[-1].set_xticklabels(xlbls,color=TEXT2,fontsize=8)
        axes[-1].set_xlabel("Posizione [m]",color=TEXT2)
        fc=FigureCanvasTkAgg(fig,master=parent)
        fc.draw(); fc.get_tk_widget().pack(fill="both",expand=True); plt.close(fig)

    def _build_modal_tab(self,parent):
        if not HAS_MPL:
            tk.Label(parent,text="Installa matplotlib",bg=BG,fg=DANGER).pack(); return
        ctrl=tk.Frame(parent,bg=PANEL2,height=40); ctrl.pack(fill="x")
        ctrl.pack_propagate(False)
        tk.Label(ctrl,text="Seleziona modo:",bg=PANEL2,fg=TEXT2,font=("Segoe UI",10)).pack(side="left",padx=8,pady=8)
        modo_cb_var=tk.IntVar(value=0)
        cb=ttk.Combobox(ctrl,textvariable=modo_cb_var,width=12,
                         values=list(range(len(self.modi) if self.modi else 0)),state="readonly")
        cb.pack(side="left",padx=4,pady=8)

        plot_frame=tk.Frame(parent,bg=BG); plot_frame.pack(fill="both",expand=True)

        def draw_mode(idx=None):
            if idx is None: idx=modo_cb_var.get()
            for w in plot_frame.winfo_children(): w.destroy()
            if self.modi is None or idx>=len(self.modi):
                tk.Label(plot_frame,text="Esegui prima l'analisi modale",bg=BG,fg=TEXT2).pack(pady=30); return
            phi=self.modi[idx]
            fig,ax=plt.subplots(1,1,figsize=(10,4.5),facecolor=BG)
            ax.set_facecolor(PANEL2); ax.set_aspect("equal")
            for sp in ax.spines.values(): sp.set_color(BORDER)
            ax.tick_params(colors=TEXT2)
            # struttura non deformata
            for el in self.elementi:
                ni,nj=self.nodi[el['i']],self.nodi[el['j']]
                ax.plot([ni['x'],nj['x']],[ni['y'],nj['y']],color=BORDER,linewidth=1,linestyle="--")
            # modo
            absP=np.abs(phi); maxP=float(absP.max()) if absP.max()>0 else 1.0
            xs_=[n['x'] for n in self.nodi]; ys_=[n['y'] for n in self.nodi]
            span=max(max(xs_)-min(xs_),max(ys_)-min(ys_),0.1)
            amp=span*0.15/maxP
            for el in self.elementi:
                if el['tipo']=='truss': continue
                ni,nj=self.nodi[el['i']],self.nodi[el['j']]
                L,theta=elem_geo(ni,nj)
                if L<1e-9: continue
                cs,sn=math.cos(theta),math.sin(theta)
                i,j=el['i'],el['j']
                u=[phi[i*3]*amp,phi[i*3+1]*amp,phi[i*3+2]*amp,
                   phi[j*3]*amp,phi[j*3+1]*amp,phi[j*3+2]*amp]
                steps=30; pts_x=[]; pts_y=[]
                for k in range(steps+1):
                    xi=k/steps
                    N1=1-3*xi**2+2*xi**3; N2=xi*(1-xi)**2*L
                    N3=3*xi**2-2*xi**3;   N4=xi**2*(xi-1)*L
                    ul_i=cs*u[0]+sn*u[1]; vl_i=-sn*u[0]+cs*u[1]
                    ul_j=cs*u[3]+sn*u[4]; vl_j=-sn*u[3]+cs*u[4]
                    ul=ul_i*(1-xi)+ul_j*xi
                    vl=N1*vl_i+N2*u[2]+N3*vl_j+N4*u[5]
                    pts_x.append(ni['x']+(nj['x']-ni['x'])*xi+cs*ul-sn*vl)
                    pts_y.append(ni['y']+(nj['y']-ni['y'])*xi+sn*ul+cs*vl)
                ax.plot(pts_x,pts_y,color=MOMENT,linewidth=2.5)
            f=self.freq[idx] if self.freq is not None and idx<len(self.freq) else 0
            T=1/f if f>1e-9 else float('inf')
            ax.set_title(f"Modo {idx+1}  —  f = {f:.4f} Hz   T = {T:.4f} s",
                         color=TEXT,fontsize=11,pad=8)
            ax.grid(True,alpha=0.1,color=BORDER)
            fc=FigureCanvasTkAgg(fig,master=plot_frame)
            fc.draw(); fc.get_tk_widget().pack(fill="both",expand=True); plt.close(fig)

        cb.bind("<<ComboboxSelected>>",lambda e:draw_mode())
        draw_mode()

    def _build_displ_tab(self,parent):
        if not HAS_MPL or self.U is None:
            tk.Label(parent,text="Nessun risultato",bg=BG,fg=TEXT2).pack(pady=20); return
        fig,axes=plt.subplots(1,3,figsize=(11,4.5),facecolor=BG)
        plt.subplots_adjust(left=0.07,right=0.97,top=0.88,bottom=0.15,wspace=0.4)
        labels=[f"N{i+1}" for i in range(len(self.nodi))]
        ux=[self.U[i*3]*1000 for i in range(len(self.nodi))]
        uy=[self.U[i*3+1]*1000 for i in range(len(self.nodi))]
        ph=[self.U[i*3+2]*1000 for i in range(len(self.nodi))]
        for ax,vals,t,col in zip(axes,[ux,uy,ph],
                                  ["ux [mm]","uy [mm]","φ [mrad]"],
                                  [SHEAR,NORMAL,MOMENT]):
            ax.set_facecolor(PANEL2)
            c2=[DANGER if v<0 else col for v in vals]
            ax.bar(labels,vals,color=c2,edgecolor=BORDER,linewidth=0.5)
            ax.set_title(t,color=TEXT2,fontsize=10,pad=5)
            ax.tick_params(colors=TEXT2); ax.axhline(0,color=BORDER,linewidth=0.8)
            for sp in ax.spines.values(): sp.set_color(BORDER)
            ax.grid(True,alpha=0.12,axis='y',linestyle="--",color=BORDER)
        fc=FigureCanvasTkAgg(fig,master=parent)
        fc.draw(); fc.get_tk_widget().pack(fill="both",expand=True); plt.close(fig)

    def _build_full_table_tab(self,parent):
        cols=("ID","Tipo","Nodi","Materiale","Sezione","L [m]","|M|max [kNm]","|N|max [kN]","η","Stato")
        tv=ttk.Treeview(parent,columns=cols,show="headings",height=20)
        style=ttk.Style()
        style.configure("Treeview",background=PANEL2,foreground=TEXT,fieldbackground=PANEL2,
                        rowheight=24,font=("Segoe UI",9))
        style.configure("Treeview.Heading",background=PANEL,foreground=TEXT2,
                        font=("Segoe UI",9,"bold"))
        style.map("Treeview",background=[("selected",ACCENT)],foreground=[("selected",BG)])
        widths=[50,60,80,130,100,70,120,120,60,80]
        for col,w in zip(cols,widths):
            tv.heading(col,text=col); tv.column(col,width=w,anchor="center")
        sb=tk.Scrollbar(parent,command=tv.yview,bg=PANEL2)
        tv.config(yscrollcommand=sb.set)
        sb.pack(side="right",fill="y"); tv.pack(fill="both",expand=True)
        for idx,el in enumerate(self.elementi):
            ni,nj=self.nodi[el['i']],self.nodi[el['j']]
            L,_=elem_geo(ni,nj)
            _,_,E,A,mat,sez=get_EI_EA(el,SECTIONS,MATERIALS)
            fy=mat['fy']; W=sez['W']
            if el['tipo']=='truss':
                N=truss_force(el,self.nodi,self.U) if self.U is not None else 0
                eta=abs(N)/(A*fy/1000) if A*fy>0 else 0
                stato="✓" if eta<1 else "✗SNERV."
                tv.insert("","end",values=(f"A{idx+1}","Asta",
                          f"N{el['i']+1}→N{el['j']+1}",el['mat'],el['sez'],
                          f"{L:.3f}","—",f"{abs(N):.3f}",f"{eta:.3f}",stato))
            else:
                if self.U is not None:
                    res=internal_forces(el,self.nodi,self.U)
                    Mmax=abs(res[1]).max() if res else 0
                    Nmax=abs(res[3]).max() if res else 0
                    sig=Mmax/W if W>0 else 0; eta=sig/fy if fy>0 else 0
                    stato="✓" if eta<1 else "✗SNERV."
                else:
                    Mmax=Nmax=eta=0; stato="—"
                tv.insert("","end",values=(f"T{idx+1}","Trave",
                          f"N{el['i']+1}→N{el['j']+1}",el['mat'],el['sez'],
                          f"{L:.3f}",f"{Mmax:.3f}",f"{Nmax:.3f}",f"{eta:.3f}",stato))

    # ── undo/redo ─────────────────────────────────────────────────────────────

    def _state(self):
        return copy.deepcopy({'nodi':self.nodi,'elementi':self.elementi,
                               'vincoli':self.vincoli,'carichi':self.carichi})
    def _restore(self,s):
        self.nodi=s['nodi']; self.elementi=s['elementi']
        self.vincoli=s['vincoli']; self.carichi=s['carichi']
        self.U=None; self._redraw()
    def _push_undo(self):
        self._undo.append(self._state())
        if len(self._undo)>60: self._undo.pop(0)
        self._redo.clear()
    def _undo_fn(self):
        if not self._undo: self.status.set("Nulla da annullare"); return
        self._redo.append(self._state()); self._restore(self._undo.pop()); self.status.set("Annullato")
    def _redo_fn(self):
        if not self._redo: self.status.set("Nulla da ripetere"); return
        self._undo.append(self._state()); self._restore(self._redo.pop()); self.status.set("Ripetuto")

    # ── file I/O ──────────────────────────────────────────────────────────────

    def _new(self):
        if messagebox.askyesno("Nuovo","Creare nuova struttura?"):
            self._push_undo(); self._clear_all(False)
            self.current_file=None; self.file_lbl.config(text="senza titolo")

    def _save(self):
        if self.current_file: self._write(self.current_file)
        else: self._save_as()

    def _save_as(self):
        p=filedialog.asksaveasfilename(defaultextension=".json",
              filetypes=[("FEM JSON","*.json")],title="Salva")
        if p: self._write(p)

    def _write(self,path):
        data={'nodi':self.nodi,'elementi':self.elementi,
              'vincoli':self.vincoli,'carichi':self.carichi,'version':'2.0'}
        with open(path,'w') as f: json.dump(data,f,indent=2)
        self.current_file=path; self.file_lbl.config(text=os.path.basename(path))
        self.status.set(f"Salvato: {path}")

    def _load(self):
        p=filedialog.askopenfilename(filetypes=[("FEM JSON","*.json")],title="Apri")
        if not p: return
        try:
            with open(p) as f: d=json.load(f)
            self._push_undo()
            self.nodi=d.get('nodi',[]); self.elementi=d.get('elementi',d.get('travi',[]))
            self.vincoli=d.get('vincoli',[]); self.carichi=d.get('carichi',[])
            # migrazione da v1
            for el in self.elementi:
                if 'EI' in el:
                    el['mat']=el.get('mat','Acciaio S235')
                    el['sez']=el.get('sez','Custom')
                    el['tipo']=el.get('tipo','trave')
            self.U=None; self.current_file=p
            self.file_lbl.config(text=os.path.basename(p))
            self._reset_view(); self.status.set(f"Aperto: {p}")
        except Exception as ex:
            messagebox.showerror("Errore",str(ex))

    def _load_preset(self,nome):
        self._push_undo()
        p=LIBRERIA[nome]
        self.nodi=copy.deepcopy(p['nodi']); self.elementi=copy.deepcopy(p['elementi'])
        self.vincoli=copy.deepcopy(p['vincoli']); self.carichi=copy.deepcopy(p['carichi'])
        self.U=None; self.current_file=None
        self.file_lbl.config(text=nome); self._reset_view()
        self.status.set(f"Caricato: {nome}")

    # ── export PDF ────────────────────────────────────────────────────────────

    def _export_pdf(self):
        if not HAS_RL:
            return messagebox.showwarning("ReportLab","pip install reportlab")
        if self.U is None:
            return messagebox.showwarning("Attenzione","Esegui prima il calcolo")
        p=filedialog.asksaveasfilename(defaultextension=".pdf",
              filetypes=[("PDF","*.pdf")],title="Esporta PDF")
        if not p: return
        try:
            self._write_pdf(p)
            messagebox.showinfo("OK",f"PDF salvato:\n{p}")
        except Exception as ex:
            messagebox.showerror("Errore PDF",str(ex))

    def _write_pdf(self,path):
        doc=SimpleDocTemplate(path,pagesize=A4,leftMargin=36,rightMargin=36,
                               topMargin=36,bottomMargin=36)
        s=getSampleStyleSheet(); story=[]
        story.append(Paragraph("FEM Solver Pro v2 — Report",s['Title']))
        story.append(Spacer(1,10))

        def mk_table(data,widths=None):
            t=Table(data,colWidths=widths)
            t.setStyle(TableStyle([
                ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#3a3a6a')),
                ('TEXTCOLOR',(0,0),(-1,0),colors.white),
                ('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
                ('FONTSIZE',(0,0),(-1,-1),8),
                ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#f0f0f8'),colors.white]),
                ('GRID',(0,0),(-1,-1),0.4,colors.HexColor('#ccccdd')),
                ('LEFTPADDING',(0,0),(-1,-1),6),
                ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
            ])); return t

        story.append(Paragraph("Nodi",s['Heading2']))
        nd=[["N°","x [m]","y [m]"]]
        for i,n in enumerate(self.nodi): nd.append([f"N{i+1}",f"{n['x']:.3f}",f"{n['y']:.3f}"])
        story.append(mk_table(nd,[40,80,80])); story.append(Spacer(1,8))

        story.append(Paragraph("Elementi",s['Heading2']))
        el_d=[["N°","Tipo","Nodi","Materiale","Sezione","L [m]"]]
        for i,el in enumerate(self.elementi):
            ni,nj=self.nodi[el['i']],self.nodi[el['j']]; L,_=elem_geo(ni,nj)
            el_d.append([f"{'A' if el['tipo']=='truss' else 'T'}{i+1}",
                         el['tipo'],f"N{el['i']+1}→N{el['j']+1}",
                         el.get('mat','—'),el.get('sez','—'),f"{L:.3f}"])
        story.append(mk_table(el_d,[30,50,70,120,100,60])); story.append(Spacer(1,8))

        story.append(Paragraph("Reazioni vincolari",s['Heading2']))
        comp=['Rx','Ry','Mz']; unit=['kN','kN','kN·m']
        rd=[["GDL","Componente","Valore","Unità"]]
        for vi in self.fixed:
            nd_=vi//3; k=vi%3
            rd.append([f"N{nd_+1}",comp[k],f"{self.R[vi]:+.6f}",unit[k]])
        story.append(mk_table(rd,[50,80,150,60])); story.append(Spacer(1,8))

        story.append(Paragraph("Spostamenti nodali",s['Heading2']))
        dd=[["Nodo","ux [m]","uy [m]","φ [rad]"]]
        for i in range(len(self.nodi)):
            dd.append([f"N{i+1}",f"{self.U[i*3]:+.6e}",
                       f"{self.U[i*3+1]:+.6e}",f"{self.U[i*3+2]:+.6e}"])
        story.append(mk_table(dd,[50,130,130,130])); story.append(Spacer(1,8))

        story.append(Paragraph("Verifica resistenza",s['Heading2']))
        vd=[["Elem.","Tipo","|M|max [kNm]","|N|max [kN]","σmax [MPa]","η","Stato"]]
        for idx,el in enumerate(self.elementi):
            _,_,E,A,mat,sez=get_EI_EA(el,SECTIONS,MATERIALS)
            fy=mat['fy']; W=sez['W']
            if el['tipo']=='truss':
                N=truss_force(el,self.nodi,self.U)
                sig=abs(N)/A if A>0 else 0; eta=sig/fy if fy>0 else 0
                vd.append([f"A{idx+1}","Asta","—",f"{abs(N):.3f}",
                           f"{sig/1e3:.2f}",f"{eta:.3f}","✓" if eta<1 else "✗"])
            else:
                res=internal_forces(el,self.nodi,self.U)
                if res:
                    _,M,V,N=res; Mmax=abs(M).max(); Nmax=abs(N).max()
                    sig=Mmax/W if W>0 else 0; eta=sig/fy if fy>0 else 0
                    vd.append([f"T{idx+1}","Trave",f"{Mmax:.3f}",f"{Nmax:.3f}",
                               f"{sig/1e3:.2f}",f"{eta:.3f}","✓" if eta<1 else "✗"])
        story.append(mk_table(vd,[40,50,100,100,100,60,60]))

        if self.freq is not None:
            story.append(Spacer(1,8))
            story.append(Paragraph("Frequenze proprie",s['Heading2']))
            md=[["Modo","f [Hz]","T [s]","ω [rad/s]"]]
            for i,(f,w) in enumerate(zip(self.freq,
                                          [math.sqrt(max(0,v)) for v in
                                           [f*(2*math.pi) for f in self.freq]])):
                T=1/f if f>1e-9 else float('inf')
                md.append([f"{i+1}",f"{f:.4f}",f"{T:.4f}",f"{f*2*math.pi:.3f}"])
            story.append(mk_table(md,[50,100,100,100]))

        doc.build(story)

    # ── utils ─────────────────────────────────────────────────────────────────

    def _del_sel(self,e=None):
        if self.sel_nodo is not None:
            self._push_undo(); i=self.sel_nodo
            self.elementi=[el for el in self.elementi if el['i']!=i and el['j']!=i]
            self.vincoli=[v for v in self.vincoli if v['nodo']!=i]
            self.carichi=[c for c in self.carichi
                          if not(c['tipo'] in('Fy','Fx','M') and c.get('nodo')==i)]
            self.nodi.pop(i)
            for el in self.elementi:
                if el['i']>i: el['i']-=1
                if el['j']>i: el['j']-=1
            for v in self.vincoli:
                if v['nodo']>i: v['nodo']-=1
            for c in self.carichi:
                if c['tipo'] in('Fy','Fx','M') and c.get('nodo',0)>i: c['nodo']-=1
            self.sel_nodo=None; self.U=None
            self._redraw()

    def _clear_all(self,confirm=True):
        if confirm and not messagebox.askyesno("Conferma","Cancellare tutto?"): return
        self.nodi.clear(); self.elementi.clear()
        self.vincoli.clear(); self.carichi.clear()
        self.U=None; self.R=None; self.freq=None; self.modi=None
        self.sel_nodo=None; self.sel_elem=None
        for t in self._rt.values(): self._set_txt(t,"")
        self.verify_lbl.config(text="—",fg=TEXT2)
        self._redraw()


# ── main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    app = FEMProV2()
    app.mainloop()
