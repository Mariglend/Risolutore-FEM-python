"""
gui_fem.py
==========
Interfaccia grafica interattiva per il solver FEM travi 2D.

Utilizzo:
    python gui_fem.py

Requisiti:
    pip install numpy matplotlib
    (tkinter è incluso in Python)

Interazione:
    - Click sinistro sul canvas → aggiungi nodo
    - Click su nodo esistente → seleziona/trascina
    - Drag tra due nodi → collega con trave
    - Pannello laterale → aggiungi vincoli e carichi
    - Tasto CANC → elimina elemento selezionato
"""

import sys, math, tkinter as tk
from tkinter import ttk, messagebox, font as tkfont
import numpy as np

# ── colori ──────────────────────────────────────────────────────────────────
BG        = "#1e1e2e"   # sfondo scuro
PANEL     = "#2a2a3e"
ACCENT    = "#7c6af7"   # viola
ACCENT2   = "#5DCAA5"   # teal
WARN      = "#EF9F27"   # amber
DANGER    = "#E24B4A"
TEXT      = "#cdd6f4"
TEXT2     = "#a6adc8"
GRID      = "#313244"
BEAM_COL  = "#7c6af7"
NODE_COL  = "#cdd6f4"
SEL_COL   = "#EF9F27"
REACT_COL = "#5DCAA5"
LOAD_COL  = "#f38ba8"
DEFORM    = "#f38ba8"

SNAP_GRID = 0.5   # m


# ── FEM core (standalone, non dipende dal pacchetto) ─────────────────────────

def k_locale(L, EI, EA):
    a = EA / L
    b = 12*EI/L**3; c = 6*EI/L**2; d = 4*EI/L; e = 2*EI/L
    K = np.zeros((6,6))
    K[0,0]=a;  K[0,3]=-a; K[3,0]=-a; K[3,3]=a
    K[1,1]=b;  K[1,2]=c;  K[1,4]=-b; K[1,5]=c
    K[2,1]=c;  K[2,2]=d;  K[2,4]=-c; K[2,5]=e
    K[4,1]=-b; K[4,2]=-c; K[4,4]=b;  K[4,5]=-c
    K[5,1]=c;  K[5,2]=e;  K[5,4]=-c; K[5,5]=d
    return K

def mat_rot(theta):
    c,s = math.cos(theta), math.sin(theta)
    T = np.zeros((6,6))
    T[0,0]=c; T[0,1]=s; T[1,0]=-s; T[1,1]=c; T[2,2]=1
    T[3,3]=c; T[3,4]=s; T[4,3]=-s; T[4,4]=c; T[5,5]=1
    return T

def solve_fem(nodi, travi, vincoli, carichi):
    n = len(nodi); ndof = 3*n
    K = np.zeros((ndof,ndof))
    F = np.zeros(ndof)

    for el in travi:
        i,j = el['i'], el['j']
        ni,nj = nodi[i], nodi[j]
        dx,dy = nj['x']-ni['x'], nj['y']-ni['y']
        L = math.hypot(dx,dy)
        if L < 1e-9: continue
        theta = math.atan2(dy,dx)
        Kl = k_locale(L, el['EI'], el['EA'])
        T  = mat_rot(theta)
        Kg = T.T @ Kl @ T
        dofs = [i*3,i*3+1,i*3+2, j*3,j*3+1,j*3+2]
        for a,da in enumerate(dofs):
            for b,db in enumerate(dofs):
                K[da,db] += Kg[a,b]

    for c in carichi:
        t = c['type']
        if t == 'Fy':   F[c['nodo']*3+1] -= c['val']
        elif t == 'Fx': F[c['nodo']*3]   -= c['val']
        elif t == 'M':  F[c['nodo']*3+2] -= c['val']
        else:
            el = travi[c['trave']]
            i,j = el['i'],el['j']
            ni,nj = nodi[i],nodi[j]
            dx,dy = nj['x']-ni['x'],nj['y']-ni['y']
            L = math.hypot(dx,dy)
            theta = math.atan2(dy,dx)
            q = c['val']; fne = np.zeros(6)
            if t=='uniforme':
                fne[1]=-q*L/2; fne[2]=-q*L**2/12
                fne[4]=-q*L/2; fne[5]=q*L**2/12
            elif t=='triang_sx':
                fne[1]=-7*q*L/20; fne[2]=-q*L**2/20
                fne[4]=-3*q*L/20; fne[5]=q*L**2/30
            elif t=='triang_dx':
                fne[1]=-3*q*L/20; fne[2]=-q*L**2/30
                fne[4]=-7*q*L/20; fne[5]=q*L**2/20
            Tv = mat_rot(theta)
            fg = Tv.T @ fne
            dofs=[i*3,i*3+1,i*3+2,j*3,j*3+1,j*3+2]
            for a,da in enumerate(dofs): F[da]+=fg[a]

    is_fixed = np.zeros(ndof, dtype=bool)
    for v in vincoli:
        nd = v['nodo']
        if v['ux']: is_fixed[nd*3]   = True
        if v['uy']: is_fixed[nd*3+1] = True
        if v['phi']:is_fixed[nd*3+2] = True

    free  = [i for i in range(ndof) if not is_fixed[i]]
    fixed = [i for i in range(ndof) if     is_fixed[i]]

    if not free:
        raise ValueError("Struttura completamente bloccata")

    Kff = K[np.ix_(free,free)]
    Ff  = F[free]

    try:
        Uf = np.linalg.solve(Kff, Ff)
    except np.linalg.LinAlgError:
        raise ValueError("Sistema singolare: struttura labile o iperstatica non supportata")

    U = np.zeros(ndof)
    for li,gi in enumerate(free): U[gi] = Uf[li]
    R = K@U - F
    return U, R, fixed


# ── App principale ────────────────────────────────────────────────────────────

class FEMApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FEM Travi 2D — Solver Interattivo")
        self.configure(bg=BG)
        self.minsize(1100, 680)

        # stato struttura
        self.nodi    = []   # {'x','y','id'}
        self.travi   = []   # {'i','j','EI','EA'}
        self.vincoli = []   # {'nodo','ux','uy','phi'}
        self.carichi = []   # {'type','nodo'/'trave','val'}

        # risultati
        self.U = None; self.R = None; self.fixed = []

        # interazione canvas
        self.mode       = tk.StringVar(value="nodo")
        self.sel_nodo   = None   # indice nodo selezionato
        self.drag_from  = None   # indice nodo da cui parte la trave
        self.drag_pos   = None
        self.dragging   = False
        self.show_deformed = tk.BooleanVar(value=False)
        self.snap_to_grid  = tk.BooleanVar(value=True)

        # vista canvas
        self.cam_x = 0.0; self.cam_y = 0.0
        self.cam_scale = 80.0  # px per metro
        self.pan_start = None

        self._build_ui()
        self._update_canvas()

    # ── costruzione UI ────────────────────────────────────────────────────────

    def _build_ui(self):
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        # pannello sinistro
        left = tk.Frame(self, bg=PANEL, width=280)
        left.grid(row=0, column=0, sticky="nsew")
        left.grid_propagate(False)
        self._build_left_panel(left)

        # canvas centrale
        right = tk.Frame(self, bg=BG)
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)
        self._build_toolbar(right)
        self._build_canvas(right)
        self._build_statusbar(right)

    def _lbl(self, parent, text, size=11, color=TEXT2, bold=False):
        f = ("Segoe UI", size, "bold" if bold else "normal")
        return tk.Label(parent, text=text, bg=PANEL, fg=color, font=f)

    def _btn(self, parent, text, cmd, color=ACCENT, fg=TEXT):
        b = tk.Button(parent, text=text, command=cmd,
                      bg=color, fg=fg, relief="flat",
                      activebackground=ACCENT2, activeforeground=BG,
                      font=("Segoe UI", 10, "bold"), padx=8, pady=4,
                      cursor="hand2")
        return b

    def _entry(self, parent, default=""):
        e = tk.Entry(parent, bg="#313244", fg=TEXT, insertbackground=TEXT,
                     relief="flat", font=("Segoe UI", 10), bd=4)
        e.insert(0, str(default))
        return e

    def _build_left_panel(self, parent):
        # titolo
        tk.Label(parent, text="⬡  FEM Solver", bg=PANEL, fg=ACCENT,
                 font=("Segoe UI", 14, "bold")).pack(pady=(14,2), padx=12, anchor="w")
        tk.Label(parent, text="Travi 2D · Euler-Bernoulli", bg=PANEL,
                 fg=TEXT2, font=("Segoe UI", 9)).pack(padx=12, anchor="w")

        sep = lambda: tk.Frame(parent, bg=GRID, height=1).pack(fill="x", pady=6)
        sep()

        # modalità
        self._lbl(parent,"MODALITÀ CURSORE",10,TEXT2,True).pack(padx=12,anchor="w")
        modes = tk.Frame(parent, bg=PANEL)
        modes.pack(fill="x", padx=12, pady=(4,0))
        for txt, val in [("+ Nodo","nodo"),("— Trave","trave"),("↑ Vincolo","vincolo"),("↓ Carico","carico"),("✥ Sposta","sposta")]:
            rb = tk.Radiobutton(modes, text=txt, variable=self.mode, value=val,
                                bg=PANEL, fg=TEXT, selectcolor=ACCENT,
                                activebackground=PANEL, activeforeground=ACCENT,
                                font=("Segoe UI",10), indicatoron=False,
                                relief="flat", padx=6, pady=3, cursor="hand2",
                                command=self._on_mode_change)
            rb.pack(fill="x", pady=1)
        sep()

        # Proprietà trave
        self._lbl(parent,"PROPRIETÀ TRAVE",10,TEXT2,True).pack(padx=12,anchor="w")
        pf = tk.Frame(parent, bg=PANEL)
        pf.pack(fill="x", padx=12, pady=4)
        self._lbl(pf,"EI (kN·m²)").grid(row=0,column=0,sticky="w",pady=2)
        self.ei_var = self._entry(pf,"10000"); self.ei_var.grid(row=0,column=1,sticky="ew",padx=(6,0))
        self._lbl(pf,"EA (kN)").grid(row=1,column=0,sticky="w",pady=2)
        self.ea_var = self._entry(pf,"100000"); self.ea_var.grid(row=1,column=1,sticky="ew",padx=(6,0))
        pf.columnconfigure(1,weight=1)
        sep()

        # Vincolo
        self._lbl(parent,"TIPO VINCOLO",10,TEXT2,True).pack(padx=12,anchor="w")
        vf = tk.Frame(parent, bg=PANEL); vf.pack(fill="x", padx=12, pady=4)
        self.vtype = tk.StringVar(value="incastro")
        for txt,val in [("⬛ Incastro","incastro"),("◯ Cerniera","cerniera"),("△ Carrello","carrello")]:
            tk.Radiobutton(vf,text=txt,variable=self.vtype,value=val,
                           bg=PANEL,fg=TEXT,selectcolor=ACCENT,
                           activebackground=PANEL,font=("Segoe UI",10),
                           indicatoron=False,relief="flat",padx=6,pady=2,
                           cursor="hand2").pack(fill="x",pady=1)
        sep()

        # Carico
        self._lbl(parent,"TIPO CARICO",10,TEXT2,True).pack(padx=12,anchor="w")
        cf = tk.Frame(parent, bg=PANEL); cf.pack(fill="x", padx=12, pady=4)
        self.ctype = tk.StringVar(value="Fy")
        for txt,val in [("↓ Forza Fy (nodo)","Fy"),("→ Forza Fx (nodo)","Fx"),
                        ("↺ Momento (nodo)","M"),("▬ Uniforme (trave)","uniforme"),
                        ("◤ Triang. sx (trave)","triang_sx"),("◥ Triang. dx (trave)","triang_dx")]:
            tk.Radiobutton(cf,text=txt,variable=self.ctype,value=val,
                           bg=PANEL,fg=TEXT,selectcolor=ACCENT,
                           activebackground=PANEL,font=("Segoe UI",9),
                           indicatoron=False,relief="flat",padx=6,pady=2,
                           cursor="hand2").pack(fill="x",pady=1)
        lf = tk.Frame(parent, bg=PANEL); lf.pack(fill="x", padx=12, pady=(0,4))
        self._lbl(lf,"Valore (kN / kNm / kN/m)").pack(anchor="w")
        self.cval = self._entry(lf,"10"); self.cval.pack(fill="x")
        sep()

        # Azioni
        self._btn(parent,"▶  RISOLVI", self._solve, ACCENT).pack(fill="x",padx=12,pady=2)
        self._btn(parent,"🗑  CANCELLA TUTTO", self._clear_all, DANGER).pack(fill="x",padx=12,pady=2)
        sep()

        # Risultati
        self._lbl(parent,"RISULTATI",10,TEXT2,True).pack(padx=12,anchor="w")
        self.result_text = tk.Text(parent, bg="#1e1e2e", fg=REACT_COL,
                                   font=("Consolas",9), relief="flat",
                                   height=10, state="disabled",
                                   insertbackground=TEXT)
        self.result_text.pack(fill="both", expand=True, padx=12, pady=(4,12))

    def _build_toolbar(self, parent):
        tb = tk.Frame(parent, bg=PANEL, height=36)
        tb.grid(row=0, column=0, sticky="ew")

        def tbtn(text, cmd, tip=""):
            b = tk.Button(tb, text=text, command=cmd,
                          bg=PANEL, fg=TEXT, relief="flat",
                          activebackground=GRID, font=("Segoe UI",9),
                          padx=8, pady=4, cursor="hand2")
            b.pack(side="left", padx=1)
            return b

        tbtn("⊕ Zoom +",  lambda: self._zoom(1.2))
        tbtn("⊖ Zoom −",  lambda: self._zoom(0.8))
        tbtn("⌂ Reset vista", self._reset_view)

        tk.Checkbutton(tb, text="Mostra deformata", variable=self.show_deformed,
                       bg=PANEL, fg=TEXT, selectcolor=ACCENT,
                       activebackground=PANEL, font=("Segoe UI",9),
                       command=self._update_canvas).pack(side="left", padx=8)
        tk.Checkbutton(tb, text="Snap griglia", variable=self.snap_to_grid,
                       bg=PANEL, fg=TEXT, selectcolor=ACCENT,
                       activebackground=PANEL, font=("Segoe UI",9)).pack(side="left")

        # legenda
        for col,txt in [(BEAM_COL,"Trave"),(REACT_COL,"Reazione"),(LOAD_COL,"Carico"),(DEFORM,"Deformata")]:
            tk.Label(tb, text="●", fg=col, bg=PANEL, font=("Segoe UI",12)).pack(side="right")
            tk.Label(tb, text=txt, fg=TEXT2, bg=PANEL, font=("Segoe UI",9)).pack(side="right")

    def _build_canvas(self, parent):
        self.canvas = tk.Canvas(parent, bg=BG, highlightthickness=0, cursor="crosshair")
        self.canvas.grid(row=1, column=0, sticky="nsew")
        self.canvas.bind("<Configure>",      lambda e: self._on_resize(e))
        self.canvas.bind("<ButtonPress-1>",  self._on_click)
        self.canvas.bind("<B1-Motion>",      self._on_drag)
        self.canvas.bind("<ButtonRelease-1>",self._on_release)
        self.canvas.bind("<ButtonPress-2>",  self._on_pan_start)
        self.canvas.bind("<B2-Motion>",      self._on_pan)
        self.canvas.bind("<ButtonPress-3>",  self._on_pan_start)
        self.canvas.bind("<B3-Motion>",      self._on_pan)
        self.canvas.bind("<MouseWheel>",     self._on_wheel)
        self.canvas.bind("<Delete>",         self._delete_selected)
        self.canvas.bind("<BackSpace>",      self._delete_selected)
        self.canvas.focus_set()

    def _build_statusbar(self, parent):
        self.status_var = tk.StringVar(value="Modalità: aggiungi nodo  |  Click sinistro per aggiungere")
        sb = tk.Label(parent, textvariable=self.status_var,
                      bg=PANEL, fg=TEXT2, font=("Segoe UI",9),
                      anchor="w", padx=10)
        sb.grid(row=2, column=0, sticky="ew")

    # ── coordinate ───────────────────────────────────────────────────────────

    def w2s(self, wx, wy):
        """World → screen."""
        W = self.canvas.winfo_width() or 800
        H = self.canvas.winfo_height() or 600
        return (W/2 + self.cam_x + wx*self.cam_scale,
                H/2 + self.cam_y - wy*self.cam_scale)

    def s2w(self, sx, sy):
        """Screen → world."""
        W = self.canvas.winfo_width() or 800
        H = self.canvas.winfo_height() or 600
        return ((sx - W/2 - self.cam_x)/self.cam_scale,
                -(sy - H/2 - self.cam_y)/self.cam_scale)

    def _snap(self, wx, wy):
        if self.snap_to_grid.get():
            g = SNAP_GRID
            return round(wx/g)*g, round(wy/g)*g
        return wx, wy

    # ── eventi canvas ─────────────────────────────────────────────────────────

    def _on_resize(self, e):
        self._update_canvas()

    def _on_mode_change(self):
        tips = {
            "nodo":    "Click → aggiungi nodo",
            "trave":   "Click su nodo A → trascina su nodo B → rilascia",
            "vincolo": "Click su nodo → applica vincolo selezionato",
            "carico":  "Click su nodo o trave → applica carico selezionato",
            "sposta":  "Trascina i nodi per riposizionarli",
        }
        m = self.mode.get()
        self.status_var.set(f"Modalità: {m}  |  {tips.get(m,'')}")
        self.sel_nodo = None
        self.drag_from = None
        self._update_canvas()

    def _find_nodo_at(self, sx, sy, tol=14):
        best, best_d = None, tol
        for i, n in enumerate(self.nodi):
            px, py = self.w2s(n['x'], n['y'])
            d = math.hypot(sx-px, sy-py)
            if d < best_d:
                best, best_d = i, d
        return best

    def _find_trave_at(self, sx, sy, tol=8):
        for i, el in enumerate(self.travi):
            ni, nj = self.nodi[el['i']], self.nodi[el['j']]
            ax,ay = self.w2s(ni['x'],ni['y'])
            bx,by = self.w2s(nj['x'],nj['y'])
            dx,dy = bx-ax, by-ay
            L2 = dx*dx+dy*dy
            if L2 < 1: continue
            t = max(0,min(1,((sx-ax)*dx+(sy-ay)*dy)/L2))
            px,py = ax+t*dx, ay+t*dy
            if math.hypot(sx-px,sy-py) < tol:
                return i
        return None

    def _on_click(self, e):
        self.canvas.focus_set()
        m = self.mode.get()
        wx, wy = self._snap(*self.s2w(e.x, e.y))

        if m == "nodo":
            # non sovrapporre nodi
            if self._find_nodo_at(e.x,e.y) is None:
                self.nodi.append({'x':wx,'y':wy,'id':len(self.nodi)})
                self.U = None
            self._update_canvas()

        elif m == "trave":
            ni = self._find_nodo_at(e.x,e.y)
            if ni is not None:
                self.drag_from = ni
                self.drag_pos  = (e.x, e.y)

        elif m == "vincolo":
            ni = self._find_nodo_at(e.x,e.y)
            if ni is not None:
                self._apply_vincolo(ni)

        elif m == "carico":
            ctype = self.ctype.get()
            if ctype in ('Fy','Fx','M'):
                ni = self._find_nodo_at(e.x,e.y)
                if ni is not None:
                    try:
                        val = float(self.cval.get())
                    except ValueError:
                        return messagebox.showerror("Errore","Valore non valido")
                    self.carichi.append({'type':ctype,'nodo':ni,'val':val})
                    self.U = None
                    self._update_canvas()
            else:
                ti = self._find_trave_at(e.x,e.y)
                if ti is not None:
                    try:
                        val = float(self.cval.get())
                    except ValueError:
                        return messagebox.showerror("Errore","Valore non valido")
                    self.carichi.append({'type':ctype,'trave':ti,'val':val})
                    self.U = None
                    self._update_canvas()

        elif m == "sposta":
            ni = self._find_nodo_at(e.x,e.y)
            self.sel_nodo = ni
            self.dragging = ni is not None
            self._update_canvas()

    def _on_drag(self, e):
        m = self.mode.get()
        if m == "trave" and self.drag_from is not None:
            self.drag_pos = (e.x, e.y)
            self._update_canvas()
        elif m == "sposta" and self.dragging and self.sel_nodo is not None:
            wx, wy = self._snap(*self.s2w(e.x, e.y))
            self.nodi[self.sel_nodo]['x'] = wx
            self.nodi[self.sel_nodo]['y'] = wy
            # ricalcola lunghezze travi
            self.U = None
            self._update_canvas()

    def _on_release(self, e):
        m = self.mode.get()
        if m == "trave" and self.drag_from is not None:
            nj = self._find_nodo_at(e.x,e.y)
            if nj is not None and nj != self.drag_from:
                # evita doppi
                exists = any(
                    (t['i']==self.drag_from and t['j']==nj) or
                    (t['i']==nj and t['j']==self.drag_from)
                    for t in self.travi
                )
                if not exists:
                    try:
                        EI = float(self.ei_var.get())
                        EA = float(self.ea_var.get())
                    except ValueError:
                        EI, EA = 10000, 100000
                    self.travi.append({'i':self.drag_from,'j':nj,'EI':EI,'EA':EA})
                    self.U = None
            self.drag_from = None
            self.drag_pos  = None
            self._update_canvas()
        elif m == "sposta":
            self.dragging = False

    def _on_pan_start(self, e):
        self.pan_start = (e.x - self.cam_x, e.y - self.cam_y)

    def _on_pan(self, e):
        if self.pan_start:
            self.cam_x = e.x - self.pan_start[0]
            self.cam_y = e.y - self.pan_start[1]
            self._update_canvas()

    def _on_wheel(self, e):
        factor = 1.15 if e.delta > 0 else 0.87
        wx0, wy0 = self.s2w(e.x, e.y)
        self.cam_scale *= factor
        self.cam_scale = max(10, min(400, self.cam_scale))
        # zoom centrato sul cursore
        W = self.canvas.winfo_width(); H = self.canvas.winfo_height()
        self.cam_x = e.x - W/2 - wx0*self.cam_scale
        self.cam_y = e.y - H/2 + wy0*self.cam_scale
        self._update_canvas()

    def _zoom(self, factor):
        W = self.canvas.winfo_width(); H = self.canvas.winfo_height()
        cx, cy = W/2, H/2
        self._on_wheel(type('E',(),{'x':cx,'y':cy,'delta':1 if factor>1 else -1})())

    def _reset_view(self):
        self.cam_x = 0; self.cam_y = 0; self.cam_scale = 80
        self._update_canvas()

    # ── vincoli ───────────────────────────────────────────────────────────────

    def _apply_vincolo(self, ni):
        # rimuovi vincolo esistente sullo stesso nodo
        self.vincoli = [v for v in self.vincoli if v['nodo'] != ni]
        vt = self.vtype.get()
        cfg = {
            'incastro':  {'ux':True,'uy':True,'phi':True},
            'cerniera':  {'ux':True,'uy':True,'phi':False},
            'carrello':  {'ux':False,'uy':True,'phi':False},
        }[vt]
        self.vincoli.append({'nodo':ni,**cfg})
        self.U = None
        self._update_canvas()

    # ── solve ─────────────────────────────────────────────────────────────────

    def _solve(self):
        if len(self.nodi) < 2:
            return messagebox.showwarning("Attenzione","Aggiungi almeno 2 nodi e 1 trave")
        if len(self.travi) < 1:
            return messagebox.showwarning("Attenzione","Aggiungi almeno 1 trave")
        if len(self.vincoli) < 1:
            return messagebox.showwarning("Attenzione","Aggiungi almeno 1 vincolo")
        try:
            U, R, fixed = solve_fem(self.nodi, self.travi, self.vincoli, self.carichi)
            self.U = U; self.R = R; self.fixed = fixed
            self._show_results(U, R, fixed)
            self._update_canvas()
        except Exception as ex:
            messagebox.showerror("Errore solver", str(ex))
            self.U = None

    def _show_results(self, U, R, fixed):
        comp = ['Rx','Ry','Mz']
        unit = ['kN','kN','kN·m']
        lines = ["═"*34, "  REAZIONI VINCOLARI", "═"*34]
        seen = set()
        for vi in fixed:
            nd = vi//3; k = vi%3
            if nd not in seen:
                lines.append(f"  Nodo {nd+1}:")
                seen.add(nd)
            val = R[vi]
            lines.append(f"    {comp[k]:3s} = {val:+10.3f}  {unit[k]}")

        lines += ["", "  SPOSTAMENTI", "─"*34]
        for i in range(len(self.nodi)):
            ux = U[i*3]; uy = U[i*3+1]; ph = U[i*3+2]
            if abs(ux)>1e-9 or abs(uy)>1e-9 or abs(ph)>1e-9:
                lines.append(f"  N{i+1}: ux={ux:+.4e}m  uy={uy:+.4e}m")

        # verifica equilibrio
        KU = np.zeros(len(U))
        for i in range(len(self.nodi)):
            for j in range(len(self.nodi)):
                pass  # già calcolato in solve
        sumFy = sum(R[i] for i in fixed if i%3==1)
        lines += ["", f"  Σ reazioni Fy = {sumFy:+.4f} kN", "═"*34]

        txt = "\n".join(lines)
        self.result_text.config(state="normal")
        self.result_text.delete("1.0","end")
        self.result_text.insert("1.0", txt)
        self.result_text.config(state="disabled")

    # ── disegno ───────────────────────────────────────────────────────────────

    def _update_canvas(self):
        c = self.canvas
        c.delete("all")
        W = c.winfo_width() or 800
        H = c.winfo_height() or 600

        self._draw_grid(W, H)
        self._draw_travi()
        if self.drag_from is not None and self.drag_pos:
            self._draw_drag_preview()
        if self.show_deformed.get() and self.U is not None:
            self._draw_deformed()
        self._draw_carichi()
        self._draw_vincoli()
        self._draw_reazioni()
        self._draw_nodi()

    def _draw_grid(self, W, H):
        c = self.canvas
        step = self._grid_step()
        tl = self.s2w(0,0); br = self.s2w(W,H)
        x0 = math.floor(min(tl[0],br[0])/step)*step
        x1 = math.ceil( max(tl[0],br[0])/step)*step
        y0 = math.floor(min(tl[1],br[1])/step)*step
        y1 = math.ceil( max(tl[1],br[1])/step)*step
        for gx in self._frange(x0,x1,step):
            sx,_ = self.w2s(gx,0)
            c.create_line(sx,0,sx,H, fill=GRID, width=0.5)
            if abs(gx) < 1e-9:
                c.create_line(sx,0,sx,H, fill="#444466", width=1)
            else:
                _,oy = self.w2s(0,0)
                c.create_text(sx, min(oy+12,H-10), text=f"{gx:.1f}",
                              fill="#555577", font=("Segoe UI",7))
        for gy in self._frange(y0,y1,step):
            _,sy = self.w2s(0,gy)
            c.create_line(0,sy,W,sy, fill=GRID, width=0.5)
            if abs(gy) < 1e-9:
                c.create_line(0,sy,W,sy, fill="#444466", width=1)
            else:
                ox,_ = self.w2s(0,0)
                c.create_text(max(ox-18,18),sy, text=f"{gy:.1f}",
                              fill="#555577", font=("Segoe UI",7))

    def _grid_step(self):
        mpp = 1/self.cam_scale
        raw = mpp*60
        exp = 10**math.floor(math.log10(raw)) if raw>0 else 1
        n = raw/exp
        if n < 2: return exp
        if n < 5: return 2*exp
        return 5*exp

    def _frange(self, start, stop, step):
        vals = []
        x = start
        while x <= stop+1e-9:
            vals.append(round(x/step)*step)
            x += step
        return vals

    def _draw_travi(self):
        c = self.canvas
        for idx, el in enumerate(self.travi):
            ni, nj = self.nodi[el['i']], self.nodi[el['j']]
            ax,ay = self.w2s(ni['x'],ni['y'])
            bx,by = self.w2s(nj['x'],nj['y'])
            c.create_line(ax,ay,bx,by, fill=BEAM_COL, width=3, capstyle="round")
            mx,my = (ax+bx)/2, (ay+by)/2
            L = math.hypot(nj['x']-ni['x'],nj['y']-ni['y'])
            c.create_text(mx,my-10, text=f"T{idx+1}  L={L:.2f}m",
                          fill=BEAM_COL, font=("Segoe UI",8))

    def _draw_drag_preview(self):
        c = self.canvas
        ni = self.nodi[self.drag_from]
        ax,ay = self.w2s(ni['x'],ni['y'])
        bx,by = self.drag_pos
        c.create_line(ax,ay,bx,by, fill=ACCENT, width=2, dash=(6,4))

    def _draw_nodi(self):
        c = self.canvas
        R = 7
        for i, n in enumerate(self.nodi):
            sx,sy = self.w2s(n['x'],n['y'])
            col = SEL_COL if i == self.sel_nodo else NODE_COL
            c.create_oval(sx-R,sy-R,sx+R,sy+R, fill=BG, outline=col, width=2)
            c.create_text(sx+12,sy-10, text=f"N{i+1}",
                          fill=col, font=("Segoe UI",8,"bold"))
            c.create_text(sx+12,sy+2, text=f"({n['x']:.1f},{n['y']:.1f})",
                          fill=TEXT2, font=("Segoe UI",7))

    def _draw_vincoli(self):
        c = self.canvas
        for v in self.vincoli:
            n = self.nodi[v['nodo']]
            sx,sy = self.w2s(n['x'],n['y'])
            vt = ('incastro' if (v['ux'] and v['uy'] and v['phi'])
                  else 'cerniera' if (v['ux'] and v['uy'])
                  else 'carrello')
            sz = 12
            if vt == 'incastro':
                c.create_rectangle(sx-sz,sy-sz,sx+sz,sy+sz,
                                   fill="", outline=ACCENT2, width=2)
                for off in range(-sz,sz+1,5):
                    c.create_line(sx+off,sy+sz,sx+off-4,sy+sz+6,
                                  fill=ACCENT2,width=1)
            elif vt == 'cerniera':
                pts = [sx,sy, sx-sz,sy+sz*1.5, sx+sz,sy+sz*1.5]
                c.create_polygon(pts, fill="", outline=ACCENT2, width=2)
                c.create_oval(sx-4,sy-4,sx+4,sy+4, fill=ACCENT2)
                for off in range(-sz,sz+1,5):
                    c.create_line(sx-sz+off,sy+sz*1.5,sx-sz+off-4,sy+sz*1.5+6,
                                  fill=ACCENT2,width=1)
            else:  # carrello
                pts = [sx,sy, sx-sz,sy+sz*1.5, sx+sz,sy+sz*1.5]
                c.create_polygon(pts, fill="", outline=ACCENT2, width=2)
                c.create_oval(sx-4,sy+sz*1.5-4,sx+4,sy+sz*1.5+4,
                              fill="", outline=ACCENT2,width=2)
                c.create_line(sx-sz-4,sy+sz*1.5+8,sx+sz+4,sy+sz*1.5+8,
                              fill=ACCENT2,width=2)

    def _draw_carichi(self):
        c = self.canvas
        arrow_len = max(25, self.cam_scale*0.4)
        for load in self.carichi:
            t = load['type']
            val = load['val']
            if t in ('Fy','Fx','M'):
                n = self.nodi[load['nodo']]
                sx,sy = self.w2s(n['x'],n['y'])
                if t == 'Fy':
                    dy = -arrow_len if val > 0 else arrow_len
                    c.create_line(sx,sy+dy,sx,sy, fill=LOAD_COL,width=2,
                                  arrow="last",arrowshape=(8,10,3))
                    c.create_text(sx+8,sy+dy/2, text=f"{abs(val):.1f}kN",
                                  fill=LOAD_COL,font=("Segoe UI",8))
                elif t == 'Fx':
                    dx = arrow_len if val > 0 else -arrow_len
                    c.create_line(sx-dx,sy,sx,sy, fill=LOAD_COL,width=2,
                                  arrow="last",arrowshape=(8,10,3))
                    c.create_text(sx-dx/2,sy-10, text=f"{abs(val):.1f}kN",
                                  fill=LOAD_COL,font=("Segoe UI",8))
                else:
                    c.create_text(sx,sy-20, text=f"↺{abs(val):.1f}kNm",
                                  fill=LOAD_COL,font=("Segoe UI",9,"bold"))
            else:
                if load.get('trave') is None or load['trave'] >= len(self.travi):
                    continue
                el = self.travi[load['trave']]
                ni,nj = self.nodi[el['i']],self.nodi[el['j']]
                ax,ay = self.w2s(ni['x'],ni['y'])
                bx,by = self.w2s(nj['x'],nj['y'])
                dx,dy = bx-ax,by-ay
                L = math.hypot(dx,dy)
                if L < 1: continue
                nx,ny = -dy/L, dx/L   # perpendicolare (verso "alto" locale)
                q = abs(val); steps = 7
                scale = min(30, arrow_len*0.5)
                c.create_line(ax,ay,bx,by, fill=WARN, width=1, dash=(3,3))
                for k in range(steps+1):
                    xi = k/steps
                    if t=='triang_sx': qk=q*(1-xi)
                    elif t=='triang_dx': qk=q*xi
                    else: qk=q
                    px,py = ax+xi*dx, ay+xi*dy
                    tx,ty = px+nx*qk/q*scale if q>0 else px, py+ny*qk/q*scale if q>0 else py
                    c.create_line(tx,ty,px,py, fill=WARN,width=1,
                                  arrow="last",arrowshape=(5,7,2))
                mx,my = (ax+bx)/2,(ay+by)/2
                c.create_text(mx+nx*scale*0.8,my+ny*scale*0.8,
                              text=f"{q:.1f}kN/m", fill=WARN,font=("Segoe UI",8))

    def _draw_reazioni(self):
        if self.U is None or self.R is None: return
        c = self.canvas
        comp = ['Rx','Ry','Mz']
        for vi in self.fixed:
            nd = vi//3; k = vi%3
            val = self.R[vi]
            if abs(val) < 1e-6: continue
            n = self.nodi[nd]
            sx,sy = self.w2s(n['x'],n['y'])
            arrow_len = max(30, self.cam_scale*0.5)
            if k == 1:  # Ry
                dy = arrow_len if val > 0 else -arrow_len
                c.create_line(sx,sy,sx,sy-dy, fill=REACT_COL,width=2,
                              arrow="last",arrowshape=(10,12,4))
                c.create_text(sx+8,sy-dy/2,
                              text=f"Ry={val:+.2f}kN",
                              fill=REACT_COL,font=("Consolas",8))
            elif k == 0:  # Rx
                dx = arrow_len if val > 0 else -arrow_len
                c.create_line(sx,sy,sx+dx,sy, fill=REACT_COL,width=2,
                              arrow="last",arrowshape=(10,12,4))
                c.create_text(sx+dx/2,sy-12,
                              text=f"Rx={val:+.2f}kN",
                              fill=REACT_COL,font=("Consolas",8))
            else:  # Mz
                c.create_text(sx+30,sy+16,
                              text=f"Mz={val:+.2f}kN·m",
                              fill=REACT_COL,font=("Consolas",8))

    def _draw_deformed(self):
        if self.U is None: return
        c = self.canvas
        U = self.U
        absU = np.abs(U)
        max_u = float(absU[absU>0].max()) if (absU>0).any() else 1.0
        span = max(
            max((n['x'] for n in self.nodi),default=1) - min((n['x'] for n in self.nodi),default=0),
            max((n['y'] for n in self.nodi),default=1) - min((n['y'] for n in self.nodi),default=0),
            0.1
        )
        amp = span*0.08/max_u

        for el in self.travi:
            ni,nj = self.nodi[el['i']],self.nodi[el['j']]
            dx,dy = nj['x']-ni['x'],nj['y']-ni['y']
            L = math.hypot(dx,dy)
            if L < 1e-9: continue
            theta = math.atan2(dy,dx)
            cs,sn = math.cos(theta),math.sin(theta)
            i,j = el['i'],el['j']
            u = [U[i*3]*amp,U[i*3+1]*amp,U[i*3+2]*amp,
                 U[j*3]*amp,U[j*3+1]*amp,U[j*3+2]*amp]
            steps = 20
            pts = []
            for k in range(steps+1):
                xi = k/steps
                N1=1-3*xi**2+2*xi**3; N2=xi*(1-xi)**2*L
                N3=3*xi**2-2*xi**3;   N4=xi**2*(xi-1)*L
                ul_i= cs*u[0]+sn*u[1]; vl_i=-sn*u[0]+cs*u[1]
                ul_j= cs*u[3]+sn*u[4]; vl_j=-sn*u[3]+cs*u[4]
                ul = ul_i*(1-xi)+ul_j*xi
                vl = N1*vl_i+N2*u[2]+N3*vl_j+N4*u[5]
                xg = ni['x']+(nj['x']-ni['x'])*xi+cs*ul-sn*vl
                yg = ni['y']+(nj['y']-ni['y'])*xi+sn*ul+cs*vl
                sx,sy = self.w2s(xg,yg)
                pts.extend([sx,sy])
            if len(pts)>=4:
                c.create_line(pts, fill=DEFORM, width=1.5,
                              dash=(5,3), smooth=True)

    # ── utilità ───────────────────────────────────────────────────────────────

    def _delete_selected(self, e=None):
        if self.sel_nodo is not None:
            i = self.sel_nodo
            self.travi   = [t for t in self.travi if t['i']!=i and t['j']!=i]
            self.vincoli = [v for v in self.vincoli if v['nodo']!=i]
            self.carichi = [c for c in self.carichi
                            if not (c['type'] in ('Fy','Fx','M') and c['nodo']==i)]
            self.nodi.pop(i)
            # aggiorna indici
            for t in self.travi:
                if t['i']>i: t['i']-=1
                if t['j']>i: t['j']-=1
            for v in self.vincoli:
                if v['nodo']>i: v['nodo']-=1
            for cc in self.carichi:
                if cc['type'] in ('Fy','Fx','M') and cc['nodo']>i:
                    cc['nodo']-=1
            self.sel_nodo=None; self.U=None
            self._update_canvas()

    def _clear_all(self):
        if messagebox.askyesno("Conferma","Cancellare tutta la struttura?"):
            self.nodi.clear(); self.travi.clear()
            self.vincoli.clear(); self.carichi.clear()
            self.U=None; self.R=None; self.sel_nodo=None
            self.result_text.config(state="normal")
            self.result_text.delete("1.0","end")
            self.result_text.config(state="disabled")
            self._update_canvas()


# ── entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    app = FEMApp()
    app.mainloop()
