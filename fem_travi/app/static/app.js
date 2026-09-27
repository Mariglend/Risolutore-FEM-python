/* FEM Travi — interfaccia grafica.
 * Il modello è nello stesso formato JSON della libreria Python (fem_travi.io):
 * convenzione y verso l'alto, forze positive nel verso degli assi, coppie antiorarie positive.
 */
'use strict';

const $ = (s) => document.querySelector(s);
const $$ = (s) => Array.from(document.querySelectorAll(s));
const svg = $('#canvas');
const NS = 'http://www.w3.org/2000/svg';

// ---------------------------------------------------------------- stato
const vuoto = () => ({ nodi: [], elementi: [], vincoli: [], carichi: [] });
let model = vuoto();
let sel = null;                 // {k: 'nodo'|'elem'|'vinc'|'car', i}
let tool = 'seleziona';
let vista = 'modello';
let risultato = null;
let view = { cx: 5, cy: 2, s: 60 };
let griglia = 0.5;
let chainStart = null;          // punto iniziale della trave in disegno {x,y}
let mouse = null;               // {sx, sy, x, y, snap}
let hover = null;
let drag = null;
let flashMsg = '', flashTimer = null;
const undoStack = [], redoStack = [];

const def = {
  trave: { EI: 1e4, EA: 1e6, inestensibile: true },
  asta: { EA: 1e6 },
  forza: { F: 10, ang: -90 },
  coppia: { M: 10, verso: 1 },
  distribuito: { qi: 10, qj: 10, direzione: 'y', verso: -1, proiezione: false },
  forza_campata: { F: 10, ang: -90 },
  coppia_campata: { M: 10, verso: 1 },
  termico: { alpha: 1.2e-5, dT: 0, dT_farfalla: 20, h: 0.5 },
};

const VINCOLI = {
  incastro: { nome: 'Incastro', blocca: 'spostamento orizzontale, verticale e rotazione', reazioni: 'Rx, Ry, M' },
  cerniera: { nome: 'Cerniera', blocca: 'spostamento orizzontale e verticale; la rotazione è libera', reazioni: 'Rx, Ry' },
  carrello: { nome: 'Carrello', blocca: 'solo lo spostamento perpendicolare al piano di scorrimento', reazioni: 'una forza perpendicolare al piano' },
  doppio_pendolo: { nome: 'Doppio pendolo', blocca: 'rotazione e spostamento perpendicolare al piano; scorre lungo il piano', reazioni: 'una forza e una coppia' },
  molla: { nome: 'Molla', blocca: 'niente in modo rigido: reagisce in proporzione allo spostamento', reazioni: 'forze elastiche k·u' },
  bloccarotazione: { nome: 'Blocco rotazione', blocca: 'solo la rotazione', reazioni: 'M' },
  personalizzato: { nome: 'Personalizzato', blocca: 'i GDL scelti', reazioni: '—' },
};
const TOOL_TIPO = {
  seleziona: 'sel', trave: 'draw', asta: 'draw',
  incastro: 'vinc', cerniera: 'vinc', carrello: 'vinc', doppio_pendolo: 'vinc', molla: 'vinc',
  cerniera_interna: 'hinge',
  forza: 'nodeload', coppia: 'nodeload',
  distribuito: 'beamload', forza_campata: 'beamload', coppia_campata: 'beamload', termico: 'beamload',
};

// ---------------------------------------------------------------- utilità
const fmt = (v, d = 3) => {
  if (!isFinite(v)) return '—';
  if (Math.abs(v) < 1e-9) return '0';
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) return v.toExponential(2).replace('.', ',');
  return (+v.toPrecision(d + (a >= 100 ? 1 : 0))).toLocaleString('it-IT', { maximumFractionDigits: 6 });
};
const num = (s, dflt = 0) => {
  const v = parseFloat(String(s).replace(',', '.'));
  return isFinite(v) ? v : dflt;
};
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const clone = (o) => JSON.parse(JSON.stringify(o));
const W = () => svg.clientWidth, H = () => svg.clientHeight;
const X = (x) => (x - view.cx) * view.s + W() / 2;
const Y = (y) => H() / 2 - (y - view.cy) * view.s;
const WX = (sx) => (sx - W() / 2) / view.s + view.cx;
const WY = (sy) => view.cy - (sy - H() / 2) / view.s;
const r2 = (v) => Math.round(v * 1e6) / 1e6;

function elemGeo(e) {
  const a = model.nodi[e.i], b = model.nodi[e.j];
  const dx = b.x - a.x, dy = b.y - a.y, L = Math.hypot(dx, dy);
  return { a, b, L, c: dx / L, s: dy / L, th: Math.atan2(dy, dx) };
}
function span() {
  if (!model.nodi.length) return 10;
  const xs = model.nodi.map((n) => n.x), ys = model.nodi.map((n) => n.y);
  return Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys), 1);
}

// ---------------------------------------------------------------- undo
function pushUndo() {
  undoStack.push(JSON.stringify(model));
  if (undoStack.length > 300) undoStack.shift();
  redoStack.length = 0;
}
function undo() {
  if (!undoStack.length) return;
  redoStack.push(JSON.stringify(model));
  model = JSON.parse(undoStack.pop());
  sel = null; chainStart = null; changed();
}
function redo() {
  if (!redoStack.length) return;
  undoStack.push(JSON.stringify(model));
  model = JSON.parse(redoStack.pop());
  sel = null; chainStart = null; changed();
}

// ---------------------------------------------------------------- modifiche al modello
function nodeAtWorld(x, y, tol = 1e-6) {
  return model.nodi.findIndex((n) => Math.abs(n.x - x) < tol && Math.abs(n.y - y) < tol);
}
function addNode(x, y) {
  const k = nodeAtWorld(x, y);
  if (k >= 0) return k;
  model.nodi.push({ x: r2(x), y: r2(y), cerniera: false });
  return model.nodi.length - 1;
}
function elemContaining(x, y) {
  for (let e = 0; e < model.elementi.length; e++) {
    const g = elemGeo(model.elementi[e]);
    const t = ((x - g.a.x) * g.c + (y - g.a.y) * g.s) / g.L;
    const d = Math.abs(-(x - g.a.x) * g.s + (y - g.a.y) * g.c);
    if (t > 1e-6 && t < 1 - 1e-6 && d < 1e-6 * Math.max(1, g.L)) return { e, t };
  }
  return null;
}
function splitElem(e, t) {
  const el = model.elementi[e];
  const g = elemGeo(el);
  const k = addNode(g.a.x + t * (g.b.x - g.a.x), g.a.y + t * (g.b.y - g.a.y));
  const nuovo = { ...clone(el), i: k, j: el.j, cerniera_i: false };
  el.j = k; el.cerniera_j = false;
  model.elementi.push(nuovo);
  const ne = model.elementi.length - 1, La = t * g.L;
  const extra = [];
  for (const c of model.carichi) {
    if (c.elem !== e) continue;
    if (c.tipo === 'distribuito') {
      const qm = c.qi + (c.qj - c.qi) * t;
      extra.push({ ...clone(c), elem: ne, qi: qm, qj: c.qj });
      c.qj = qm;
    } else if (c.tipo === 'termico') {
      extra.push({ ...clone(c), elem: ne });
    } else if ('a' in c && c.a > La + 1e-9) {
      c.elem = ne; c.a = r2(c.a - La);
    }
  }
  model.carichi.push(...extra);
  return k;
}
function realize(p) {
  let k = nodeAtWorld(p.x, p.y);
  if (k >= 0) return k;
  const hit = elemContaining(p.x, p.y);
  if (hit) return splitElem(hit.e, hit.t);
  return addNode(p.x, p.y);
}
function addElem(i, j, asta) {
  if (i === j) return -1;
  if (model.elementi.some((e) => (e.i === i && e.j === j) || (e.i === j && e.j === i))) return -1;
  const el = asta
    ? { i, j, tipo: 'asta', EA: def.asta.EA, EI: 1e4 }
    : { i, j, tipo: 'trave', EI: def.trave.EI, EA: def.trave.EA, inestensibile: def.trave.inestensibile };
  el.cerniera_i = false; el.cerniera_j = false;
  model.elementi.push(el);
  return model.elementi.length - 1;
}
function removeElem(e) {
  model.carichi = model.carichi.filter((c) => c.elem !== e);
  model.carichi.forEach((c) => { if (c.elem > e) c.elem--; });
  model.elementi.splice(e, 1);
}
function removeNode(k) {
  for (let e = model.elementi.length - 1; e >= 0; e--) {
    const el = model.elementi[e];
    if (el.i === k || el.j === k) removeElem(e);
  }
  model.vincoli = model.vincoli.filter((v) => v.nodo !== k);
  model.carichi = model.carichi.filter((c) => c.nodo !== k);
  model.nodi.splice(k, 1);
  const fix = (o, key) => { if (o[key] > k) o[key]--; };
  model.elementi.forEach((el) => { fix(el, 'i'); fix(el, 'j'); });
  model.vincoli.forEach((v) => fix(v, 'nodo'));
  model.carichi.forEach((c) => { if ('nodo' in c) fix(c, 'nodo'); });
}
function pruneNodes() {
  for (let k = model.nodi.length - 1; k >= 0; k--) {
    const used = model.elementi.some((e) => e.i === k || e.j === k);
    if (!used) removeNode(k);
  }
}
function deleteSel() {
  if (!sel) return;
  pushUndo();
  if (sel.k === 'nodo') removeNode(sel.i);
  else if (sel.k === 'elem') { removeElem(sel.i); pruneNodes(); }
  else if (sel.k === 'vinc') model.vincoli.splice(sel.i, 1);
  else if (sel.k === 'car') model.carichi.splice(sel.i, 1);
  sel = null;
  changed();
}

function autoAngolo(k, tipo) {
  if (tipo !== 'incastro') return 0;
  // parete dalla parte opposta alle travi che arrivano nel nodo
  let dx = 0, dy = 0;
  for (const e of model.elementi) {
    if (e.i !== k && e.j !== k) continue;
    const o = model.nodi[e.i === k ? e.j : e.i], n = model.nodi[k];
    const L = Math.hypot(o.x - n.x, o.y - n.y);
    dx += (o.x - n.x) / L; dy += (o.y - n.y) / L;
  }
  if (Math.hypot(dx, dy) < 1e-9) return 0;
  const a = Math.atan2(dy, dx) * 180 / Math.PI - 90;
  return Math.round(a / 15) * 15;
}
function setVincolo(k, tipo) {
  pushUndo();
  const iv = model.vincoli.findIndex((v) => v.nodo === k);
  if (iv >= 0 && model.vincoli[iv].tipo === tipo) {
    const v = model.vincoli[iv];
    v.angolo = ((v.angolo || 0) + 90 + 180) % 360 - 180;
    flash(`${VINCOLI[tipo].nome} ruotato a ${v.angolo}°`);
    sel = { k: 'vinc', i: iv };
  } else {
    const v = { nodo: k, tipo, angolo: autoAngolo(k, tipo), kx: 0, ky: 0, kphi: 0, cedimento: { x: 0, y: 0, phi: 0 } };
    if (tipo === 'molla') v.ky = 1e3;
    if (iv >= 0) model.vincoli[iv] = v; else model.vincoli.push(v);
    sel = { k: 'vinc', i: iv >= 0 ? iv : model.vincoli.length - 1 };
  }
  changed();
}
function addCarico(c) {
  pushUndo();
  model.carichi.push(c);
  sel = { k: 'car', i: model.carichi.length - 1 };
  changed();
}
function caricoDaTool(t, target) {
  const d = def[t];
  if (t === 'forza') {
    const a = d.ang * Math.PI / 180;
    return { tipo: 'forza', nodo: target, Fx: r2(d.F * Math.cos(a)), Fy: r2(d.F * Math.sin(a)) };
  }
  if (t === 'coppia') return { tipo: 'coppia', nodo: target, M: d.M * d.verso };
  if (t === 'distribuito') return { tipo: 'distribuito', elem: target.e, qi: d.qi * d.verso, qj: d.qj * d.verso, direzione: d.direzione, proiezione: d.proiezione };
  if (t === 'forza_campata') {
    const a = d.ang * Math.PI / 180;
    return { tipo: 'forza_campata', elem: target.e, a: target.a, Fx: r2(d.F * Math.cos(a)), Fy: r2(d.F * Math.sin(a)) };
  }
  if (t === 'coppia_campata') return { tipo: 'coppia_campata', elem: target.e, a: target.a, M: d.M * d.verso };
  if (t === 'termico') return { tipo: 'termico', elem: target.e, ...clone(d) };
}

// ---------------------------------------------------------------- snapping e hit test
function nodeNear(sx, sy, tol = 12) {
  let best = -1, bd = tol;
  model.nodi.forEach((n, k) => {
    const d = Math.hypot(X(n.x) - sx, Y(n.y) - sy);
    if (d < bd) { bd = d; best = k; }
  });
  return best;
}
function elemNear(sx, sy, tol = 9) {
  let best = null, bd = tol;
  model.elementi.forEach((el, e) => {
    const g = elemGeo(el);
    const ax = X(g.a.x), ay = Y(g.a.y), bx = X(g.b.x), by = Y(g.b.y);
    const vx = bx - ax, vy = by - ay, l2 = vx * vx + vy * vy;
    if (l2 < 1) return;
    const t = Math.max(0, Math.min(1, ((sx - ax) * vx + (sy - ay) * vy) / l2));
    const d = Math.hypot(ax + t * vx - sx, ay + t * vy - sy);
    if (d < bd) { bd = d; best = { e, t }; }
  });
  return best;
}
function snapPoint(sx, sy, ev) {
  const k = nodeNear(sx, sy);
  if (k >= 0) return { x: model.nodi[k].x, y: model.nodi[k].y, nodo: k };
  let x = WX(sx), y = WY(sy);
  const free = ev && ev.altKey;
  if (chainStart && ev && ev.shiftKey) {          // vincola a 0/45/90°
    const dx = x - chainStart.x, dy = y - chainStart.y;
    const L = Math.hypot(dx, dy), a = Math.round(Math.atan2(dy, dx) / (Math.PI / 4)) * Math.PI / 4;
    const Ls = free ? L : Math.max(griglia, Math.round(L / griglia) * griglia);
    return { x: r2(chainStart.x + Ls * Math.cos(a)), y: r2(chainStart.y + Ls * Math.sin(a)) };
  }
  const h = elemNear(sx, sy, 8);
  if (h && h.t > 0.001 && h.t < 0.999) {
    const g = elemGeo(model.elementi[h.e]);
    let d = h.t * g.L;
    if (!free) d = Math.round(d / griglia) * griglia;
    if (d > 1e-6 && d < g.L - 1e-6) return { x: r2(g.a.x + g.c * d), y: r2(g.a.y + g.s * d), elem: h.e };
  }
  if (!free) { x = Math.round(x / griglia) * griglia; y = Math.round(y / griglia) * griglia; }
  return { x: r2(x), y: r2(y) };
}
function beamTarget(sx, sy) {
  const h = elemNear(sx, sy, 12);
  if (!h) return null;
  const g = elemGeo(model.elementi[h.e]);
  let a = Math.round((h.t * g.L) / griglia) * griglia;
  a = Math.max(0, Math.min(g.L, a));
  if (Math.abs(a - g.L) < 1e-9) a = g.L;
  return { e: h.e, a: r2(a) };
}

// ---------------------------------------------------------------- calcolo
let solveTimer = null, solveSeq = 0;
function schedule() { clearTimeout(solveTimer); solveTimer = setTimeout(solve, 110); }
async function solve() {
  const seq = ++solveSeq;
  if (!model.elementi.length) { risultato = { vuoto: true, msg: 'Disegna almeno una trave.' }; renderAll(); return; }
  if (!model.vincoli.length) { risultato = { vuoto: true, msg: 'Aggiungi i vincoli: scegli un vincolo a sinistra e clicca un nodo.' }; renderAll(); return; }
  try {
    const res = await fetch('api/risolvi', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(model) });
    const out = await res.json();
    if (seq !== solveSeq) return;
    risultato = out;
  } catch (e) {
    if (seq !== solveSeq) return;
    risultato = { ok: false, errore: 'Il server di calcolo non risponde. È ancora aperta la finestra di "fem-travi"?' };
  }
  renderAll();
}
function changed() {
  save();
  renderAll();
  schedule();
}
function save() { try { localStorage.setItem('fem-travi-modello', JSON.stringify(model)); } catch (e) { /* storage non disponibile */ } }

// ---------------------------------------------------------------- disegno SVG
const el = (tag, attrs, inner = '') => {
  let s = `<${tag}`;
  for (const k in attrs) if (attrs[k] !== undefined && attrs[k] !== null) s += ` ${k}="${attrs[k]}"`;
  return s + (inner === null ? '/>' : `>${inner}</${tag}>`);
};
function arrow(x1, y1, x2, y2, cls = '', head = 9) {
  const a = Math.atan2(y2 - y1, x2 - x1);
  const hx = x2 - head * Math.cos(a), hy = y2 - head * Math.sin(a);
  const p1 = `${hx + head * 0.45 * Math.sin(a)},${hy - head * 0.45 * Math.cos(a)}`;
  const p2 = `${hx - head * 0.45 * Math.sin(a)},${hy + head * 0.45 * Math.cos(a)}`;
  return `<path class="l ${cls}" d="M${x1},${y1}L${hx},${hy}"/><path d="M${x2},${y2}L${p1}L${p2}Z" stroke-width="1"/>`;
}
function arcArrow(cx, cy, r, ccwScreen) {
  // arco di 270°; ccwScreen = antiorario come appare sullo schermo
  const a0 = ccwScreen ? 0.35 : Math.PI - 0.35, sweep = 1.5 * Math.PI;
  const a1 = ccwScreen ? a0 - sweep : a0 + sweep;
  const p = (a) => [cx + r * Math.cos(a), cy + r * Math.sin(a)];
  const [x0, y0] = p(a0), [x1, y1] = p(a1);
  const d = `M${x0},${y0}A${r},${r} 0 1 ${ccwScreen ? 0 : 1} ${x1},${y1}`;
  const tang = ccwScreen ? a1 - Math.PI / 2 : a1 + Math.PI / 2;
  const hx = x1 + 8 * Math.cos(tang + Math.PI), hy = y1 + 8 * Math.sin(tang + Math.PI);
  const n = [Math.cos(tang + Math.PI / 2) * 4, Math.sin(tang + Math.PI / 2) * 4];
  return `<path class="l" d="${d}"/><path d="M${x1 + 3 * Math.cos(tang)},${y1 + 3 * Math.sin(tang)}L${hx + n[0]},${hy + n[1]}L${hx - n[0]},${hy - n[1]}Z" stroke-width="1"/>`;
}
function hatch(y, w) {
  let s = `<path d="M${-w},${y}H${w}"/>`;
  for (let x = -w + 4; x <= w; x += 6) s += `<path d="M${x},${y}l-5,6" stroke-width="1"/>`;
  return s;
}
function simboloVincolo(tipo) {
  switch (tipo) {
    case 'incastro': return `<path d="M-18,0H18" stroke-width="4"/>${hatch(2, 18).replace(/<path d="M-18,2H18"\/>/, '')}`;
    case 'cerniera': return `<path class="fillc" d="M0,0L-11,17H11Z"/>${hatch(19, 16)}<circle class="fillc" r="3.5"/>`;
    case 'carrello': return `<path class="fillc" d="M0,0L-11,15H11Z"/><circle class="fillc" cx="-6" cy="19" r="3"/><circle class="fillc" cx="6" cy="19" r="3"/>${hatch(23, 16)}<circle class="fillc" r="3.5"/>`;
    case 'doppio_pendolo': return `<path d="M-12,0H12" stroke-width="3"/><path d="M-8,0V18M8,0V18"/><circle class="fillc" cx="-8" cy="18" r="2.5"/><circle class="fillc" cx="8" cy="18" r="2.5"/>${hatch(21, 16)}`;
    case 'molla': return `<path d="M0,0V5L-7,8L7,12L-7,16L7,20L0,23V27"/>${hatch(27, 12)}`;
    case 'bloccarotazione': return `<rect class="fillc" x="-7" y="-7" width="14" height="14"/>`;
    default: return `<rect class="fillc" x="-8" y="0" width="16" height="16"/>${hatch(18, 12)}`;
  }
}

function drawGrid(w, h) {
  let step = griglia;
  while (step * view.s < 9) step *= step < 1 ? (String(step).startsWith('0.25') ? 2 : 2.5) : 2;
  const major = step >= 1 ? step * 5 : 1;
  const x0 = WX(0), x1 = WX(w), y0 = WY(h), y1 = WY(0);
  let s = '';
  let tick = '';
  for (let x = Math.ceil(x0 / step) * step; x <= x1; x += step) {
    const isM = Math.abs(x / major - Math.round(x / major)) < 1e-6;
    s += `<path class="${isM ? 'grid2' : 'grid1'}" d="M${X(x)},0V${h}"/>`;
    if (isM) tick += `<text class="ticks" x="${X(x) + 2}" y="${h - 22}">${fmt(r2(x))}</text>`;
  }
  for (let y = Math.ceil(y0 / step) * step; y <= y1; y += step) {
    const isM = Math.abs(y / major - Math.round(y / major)) < 1e-6;
    s += `<path class="${isM ? 'grid2' : 'grid1'}" d="M0,${Y(y)}H${w}"/>`;
    if (isM) tick += `<text class="ticks" x="4" y="${Y(y) - 2}">${fmt(r2(y))}</text>`;
  }
  s += `<path class="axis" d="M${X(0)},0V${h}M0,${Y(0)}H${w}"/>`;
  return s + tick;
}

function drawCarichi() {
  let s = '';
  const qmax = Math.max(1e-12, ...model.carichi.filter((c) => c.tipo === 'distribuito').map((c) => Math.max(Math.abs(c.qi), Math.abs(c.qj))));
  const stack = {};
  model.carichi.forEach((c, k) => {
    const isSel = sel && sel.k === 'car' && sel.i === k;
    let g = '';
    if (c.tipo === 'forza' || c.tipo === 'forza_campata') {
      let px, py;
      if (c.tipo === 'forza') { const n = model.nodi[c.nodo]; px = X(n.x); py = Y(n.y); }
      else { const gg = elemGeo(model.elementi[c.elem]); px = X(gg.a.x + gg.c * c.a); py = Y(gg.a.y + gg.s * c.a); }
      const F = Math.hypot(c.Fx, c.Fy);
      if (F < 1e-12) return;
      const ux = c.Fx / F, uy = -c.Fy / F;
      g += arrow(px - 52 * ux, py - 52 * uy, px - 5 * ux, py - 5 * uy);
      g += `<text x="${px - 58 * ux + (Math.abs(uy) > 0.5 ? 6 : -14)}" y="${py - 58 * uy + (Math.abs(ux) > 0.5 ? -6 : 4)}">${fmt(F)} kN</text>`;
    } else if (c.tipo === 'coppia' || c.tipo === 'coppia_campata') {
      let px, py;
      if (c.tipo === 'coppia') { const n = model.nodi[c.nodo]; px = X(n.x); py = Y(n.y); }
      else { const gg = elemGeo(model.elementi[c.elem]); px = X(gg.a.x + gg.c * c.a); py = Y(gg.a.y + gg.s * c.a); }
      g += arcArrow(px, py, 17, c.M > 0);
      g += `<text x="${px + 20}" y="${py - 18}">${fmt(Math.abs(c.M))} kN·m</text>`;
    } else if (c.tipo === 'distribuito') {
      const gg = elemGeo(model.elementi[c.elem]);
      const sg = Math.sign(c.qi || c.qj || -1);
      let dx, dy;
      if (c.direzione === 'y') { dx = 0; dy = sg; }
      else if (c.direzione === 'x') { dx = sg; dy = 0; }
      else if (c.direzione === 'perp') { dx = -gg.s * sg; dy = gg.c * sg; }
      else { dx = gg.c * sg; dy = gg.s * sg; }
      const sdx = dx, sdy = -dy;
      const lev = stack[c.elem] = (stack[c.elem] || 0) + 1;
      const off = (lev - 1) * 48;
      const ax = X(gg.a.x), ay = Y(gg.a.y), bx = X(gg.b.x), by = Y(gg.b.y);
      const lenpx = Math.hypot(bx - ax, by - ay);
      const n = Math.max(3, Math.min(16, Math.round(lenpx / 28)));
      const hq = (t) => 8 + 30 * Math.abs(c.qi + (c.qj - c.qi) * t) / qmax;
      let env = '';
      for (let m = 0; m <= n; m++) {
        const t = m / n, x = ax + t * (bx - ax), y = ay + t * (by - ay);
        const h = hq(t);
        const bx0 = x - sdx * off, by0 = y - sdy * off;
        const tx = bx0 - sdx * h, ty = by0 - sdy * h;
        if (Math.abs(c.qi + (c.qj - c.qi) * t) > 1e-9) g += arrow(tx, ty, bx0 - sdx * 3, by0 - sdy * 3, '', 7);
        env += `${m ? 'L' : 'M'}${tx},${ty}`;
      }
      g += `<path class="l" d="${env}"/>`;
      const tm = 0.5, xm = (ax + bx) / 2 - sdx * (off + hq(tm) + 12), ym = (ay + by) / 2 - sdy * (off + hq(tm) + 12);
      const lab = Math.abs(c.qi - c.qj) < 1e-12 ? `q = ${fmt(Math.abs(c.qi))} kN/m` : `q = ${fmt(Math.abs(c.qi))} → ${fmt(Math.abs(c.qj))} kN/m`;
      g += `<text x="${xm}" y="${ym + 4}" text-anchor="middle">${lab}${c.proiezione ? ' (proiez.)' : ''}</text>`;
    } else if (c.tipo === 'termico') {
      const gg = elemGeo(model.elementi[c.elem]);
      const xm = X((gg.a.x + gg.b.x) / 2), ym = Y((gg.a.y + gg.b.y) / 2);
      const parts = [];
      if (c.dT) parts.push(`ΔT=${fmt(c.dT)}°`);
      if (c.dT_farfalla) parts.push(`ΔT↕=${fmt(c.dT_farfalla)}°`);
      g += `<text x="${xm + 12 * gg.s}" y="${ym + 16 * gg.c + 4}" text-anchor="middle">🌡 ${parts.join(' ') || 'ΔT=0'}</text>`;
    }
    s += `<g class="load${isSel ? ' sel' : ''}" data-k="car" data-i="${k}">${g}</g>`;
  });
  return s;
}

function drawVincoli() {
  let s = '';
  model.vincoli.forEach((v, k) => {
    const n = model.nodi[v.nodo];
    const isSel = sel && sel.k === 'vinc' && sel.i === k;
    s += `<g class="vinc${isSel ? ' sel' : ''}" data-k="vinc" data-i="${k}" transform="translate(${X(n.x)},${Y(n.y)}) rotate(${-(v.angolo || 0)})">${simboloVincolo(v.tipo)}<rect x="-20" y="-6" width="40" height="36" fill="transparent" stroke="none"/></g>`;
  });
  return s;
}

function drawStruttura(faint) {
  let s = '';
  model.elementi.forEach((e, k) => {
    const g = elemGeo(e);
    const isSel = sel && sel.k === 'elem' && sel.i === k;
    const isHov = hover && hover.k === 'elem' && hover.i === k;
    const cls = `beam${e.tipo === 'asta' ? ' asta' : ''}${isSel ? ' sel' : ''}${isHov ? ' hov' : ''}`;
    s += `<g data-k="elem" data-i="${k}"${faint ? ' opacity="0.35"' : ''}>`;
    s += `<path class="${cls}" d="M${X(g.a.x)},${Y(g.a.y)}L${X(g.b.x)},${Y(g.b.y)}"/>`;
    s += `<path class="beam-hit" d="M${X(g.a.x)},${Y(g.a.y)}L${X(g.b.x)},${Y(g.b.y)}"/>`;
    for (const [flag, node, dir] of [[e.cerniera_i, g.a, 1], [e.cerniera_j, g.b, -1]]) {
      if (flag && !node.cerniera) s += `<circle class="node hinge" cx="${X(node.x) + dir * 9 * g.c}" cy="${Y(node.y) - dir * 9 * g.s}" r="4"/>`;
    }
    s += `</g>`;
  });
  model.nodi.forEach((n, k) => {
    const soloAste = model.elementi.some((e) => e.i === k || e.j === k) && model.elementi.filter((e) => e.i === k || e.j === k).every((e) => e.tipo === 'asta');
    const hinge = n.cerniera || soloAste;
    const isSel = sel && sel.k === 'nodo' && sel.i === k;
    const isHov = hover && hover.k === 'nodo' && hover.i === k;
    s += `<circle class="node${hinge ? ' hinge' : ''}${isSel ? ' sel' : ''}${isHov ? ' hov' : ''}" data-k="nodo" data-i="${k}" cx="${X(n.x)}" cy="${Y(n.y)}" r="${hinge ? 5 : 4}" ${isSel ? 'stroke-width="3"' : ''}/>`;
    s += `<text class="lbl" x="${X(n.x) + 7}" y="${Y(n.y) - 7}">${k + 1}</text>`;
  });
  return s;
}

function diagrammi(q) {
  const R = risultato;
  let vmax = 1e-12;
  R.elementi.forEach((e) => e[q].forEach((v) => { vmax = Math.max(vmax, Math.abs(v)); }));
  const hpx = 70;                                    // ampiezza massima in pixel
  const sgn = q === 'M' ? -1 : 1;                    // M dal lato delle fibre tese
  let s = '', t = '';
  const labels = [];
  R.elementi.forEach((e, k) => {
    if (!model.elementi[k]) return;
    const g = elemGeo(model.elementi[k]);
    const nx = -g.s, ny = g.c;
    if (model.elementi[k].tipo === 'asta') {
      // aste: N costante, si colora l'asta (blu = tirante, rosso = puntone) invece della fascia
      if (q !== 'N') return;
      const N = e.N[0];
      const cls = Math.abs(N) < 1e-6 * vmax ? 'zero' : (N > 0 ? 'tira' : 'punt');
      s += `<path class="asta-N ${cls}" d="M${X(g.a.x)},${Y(g.a.y)}L${X(g.b.x)},${Y(g.b.y)}"/>`;
      if (cls !== 'zero') {
        const mx = X((g.a.x + g.b.x) / 2), my = Y((g.a.y + g.b.y) / 2);
        labels.push([mx, my - 3, fmt(N), `asta ${cls}`]);
      }
      return;
    }
    let poly = '', line = '';
    const pts = e.x.map((x, m) => {
      const bxw = g.a.x + g.c * x, byw = g.a.y + g.s * x;
      const off = sgn * e[q][m] / vmax * hpx;
      return [X(bxw), Y(byw), X(bxw) + nx * off, Y(byw) - ny * off];
    });
    poly = `M${pts[0][0]},${pts[0][1]}` + pts.map((p) => `L${p[2]},${p[3]}`).join('') + `L${pts[pts.length - 1][0]},${pts[pts.length - 1][1]}Z`;
    line = pts.map((p, m) => `${m ? 'L' : 'M'}${p[2]},${p[3]}`).join('');
    s += `<path class="dg ${q}" d="${poly}" stroke="none"/><path class="dg ${q}" d="${line}" fill="none"/>`;
    // etichette: estremi della trave e massimo interno
    const vals = e[q];
    const idx = new Set([0, vals.length - 1]);
    // massimi/minimi locali interni
    for (let m = 1; m < vals.length - 1; m++) {
      const d1 = vals[m] - vals[m - 1], d2 = vals[m + 1] - vals[m];
      if (d1 * d2 < 0 && Math.abs(d1) > 1e-6 * vmax && Math.abs(d2) > 1e-6 * vmax) idx.add(m);
    }
    // salti (forze/coppie in campata): valore a sinistra e a destra solo se diversi
    for (let m = 1; m < e.x.length; m++) {
      if (Math.abs(e.x[m] - e.x[m - 1]) < 1e-9) {
        idx.add(m - 1);
        if (Math.abs(vals[m] - vals[m - 1]) > 1e-6 * vmax) idx.add(m);
      }
    }
    const scritti = [];
    [...idx].sort((a, b) => a - b).forEach((m) => {
      const v = vals[m];
      if (Math.abs(v) < 1e-6 * vmax) return;
      if (scritti.some((k) => Math.abs(e.x[k] - e.x[m]) < 1e-9 && Math.abs(vals[k] - v) < 1e-6 * vmax)) return;
      scritti.push(m);
      const p = pts[m];
      const d = Math.sign(sgn * v) || 1;
      labels.push([p[2] + nx * 12 * d, p[3] - ny * 12 * d, fmt(v)]);
    });
  });
  const placed = [];
  labels.forEach(([x, y, txt, extra]) => {
    if (placed.some(([px, py]) => Math.abs(px - x) < 34 && Math.abs(py - y) < 13)) y += 14;
    placed.push([x, y]);
    t += `<text class="dgt ${q} ${extra || ''}" x="${x}" y="${y + 4}" text-anchor="middle">${txt}</text>`;
  });
  return s + t;
}

function scalaDeformata() {
  const R = risultato;
  let umax = 1e-300;
  R.elementi.forEach((e) => {
    e.def_u.forEach((u, m) => { umax = Math.max(umax, Math.hypot(u, e.def_v[m])); });
  });
  const base = 0.1 * span() / umax;
  return base * Math.pow(10, num($('#scala').value));
}
function deformata() {
  const R = risultato, k = scalaDeformata();
  let s = '';
  R.elementi.forEach((e, t) => {
    if (!model.elementi[t]) return;
    const g = elemGeo(model.elementi[t]);
    const d = e.def_x.map((x, m) => {
      const u = e.def_u[m], v = e.def_v[m];
      return `${m ? 'L' : 'M'}${X(g.a.x + g.c * x + k * (g.c * u - g.s * v))},${Y(g.a.y + g.s * x + k * (g.s * u + g.c * v))}`;
    }).join('');
    s += `<path class="def" d="${d}"/>`;
  });
  return s;
}
function drawReazioni() {
  const R = risultato;
  let s = '';
  let rmax = 1e-12;
  R.reazioni.forEach((r) => { rmax = Math.max(rmax, Math.abs(r.Rx), Math.abs(r.Ry)); });
  R.reazioni.forEach((r) => {
    const n = model.nodi[r.nodo];
    const v = model.vincoli.find((vv) => vv.nodo === r.nodo) || { angolo: 0 };
    if (!n) return;
    const px = X(n.x), py = Y(n.y);
    const a = (v.angolo || 0) * Math.PI / 180;
    const comps = Math.abs(v.angolo || 0) > 1e-9
      ? [[r.R_parallela, Math.cos(a), Math.sin(a)], [r.R_normale, -Math.sin(a), Math.cos(a)]]
      : [[r.Rx, 1, 0], [r.Ry, 0, 1]];
    let g = '';
    comps.forEach(([val, ex, ey]) => {
      if (Math.abs(val) < 1e-7 * Math.max(1, rmax)) return;
      const ux = Math.sign(val) * ex, uy = -Math.sign(val) * ey;     // direzione su schermo
      const t0 = 78, t1 = 34;
      g += arrow(px - t0 * ux, py - t0 * uy, px - t1 * ux, py - t1 * uy, '', 10);
      g += `<text x="${px - (t0 + 6) * ux + (Math.abs(uy) > 0.5 ? 8 : -12)}" y="${py - (t0 + 6) * uy + 4}">${fmt(Math.abs(val))}</text>`;
    });
    if (Math.abs(r.M) > 1e-7 * Math.max(1, rmax)) {
      g += arcArrow(px, py, 28, r.M > 0);
      g += `<text x="${px + 30}" y="${py + 34}">${fmt(Math.abs(r.M))}</text>`;
    }
    s += `<g class="reac">${g}</g>`;
  });
  return s;
}

function drawMeccanismo(fase) {
  const R = risultato;
  if (!R || !R.meccanismo) return '';
  const amp = 0.12 * span() * Math.sin(fase);
  let s = '';
  model.elementi.forEach((e) => {
    const a = model.nodi[e.i], b = model.nodi[e.j];
    const ua = R.meccanismo.slice(3 * e.i, 3 * e.i + 2), ub = R.meccanismo.slice(3 * e.j, 3 * e.j + 2);
    s += `<path class="mech" d="M${X(a.x + amp * ua[0])},${Y(a.y + amp * ua[1])}L${X(b.x + amp * ub[0])},${Y(b.y + amp * ub[1])}"/>`;
  });
  return s;
}

function drawPreview() {
  if (!mouse) return '';
  let s = '';
  const TT = TOOL_TIPO[tool];
  if (TT === 'draw' && mouse.snap) {
    const p = mouse.snap;
    s += `<circle class="snap" cx="${X(p.x)}" cy="${Y(p.y)}" r="${p.nodo !== undefined ? 9 : 6}"/>`;
    if (chainStart) {
      const L = Math.hypot(p.x - chainStart.x, p.y - chainStart.y);
      const ang = Math.atan2(p.y - chainStart.y, p.x - chainStart.x) * 180 / Math.PI;
      s += `<path class="preview" d="M${X(chainStart.x)},${Y(chainStart.y)}L${X(p.x)},${Y(p.y)}"/>`;
      s += `<text class="plabel" x="${(X(chainStart.x) + X(p.x)) / 2 + 10}" y="${(Y(chainStart.y) + Y(p.y)) / 2 - 10}">L = ${fmt(L)} m · ${fmt(ang)}°</text>`;
    }
  }
  if (TT === 'beamload' && mouse.beam && (tool === 'forza_campata' || tool === 'coppia_campata')) {
    const g = elemGeo(model.elementi[mouse.beam.e]);
    const x = X(g.a.x + g.c * mouse.beam.a), y = Y(g.a.y + g.s * mouse.beam.a);
    s += `<circle class="snap" cx="${x}" cy="${y}" r="6"/><text class="plabel" x="${x + 10}" y="${y - 12}">a = ${fmt(mouse.beam.a)} m</text>`;
  }
  return s;
}

let fase = 0, animId = null;
function renderAll() { render(); renderPanel(); renderRisultati(); renderHint(); }
function render() {
  const w = W(), h = H();
  svg.setAttribute('viewBox', `0 0 ${w} ${h}`);
  const R = risultato && risultato.ok ? risultato : null;
  const risultatiVisti = R && vista !== 'modello';
  let s = `<g>${drawGrid(w, h)}</g>`;
  if (risultatiVisti && ['N', 'V', 'M'].includes(vista)) s += `<g>${diagrammi(vista)}</g>`;
  s += `<g>${drawStruttura(risultatiVisti && vista === 'deformata')}</g>`;
  s += `<g>${drawVincoli()}</g>`;
  if (!risultatiVisti || vista === 'deformata') s += `<g>${drawCarichi()}</g>`;
  if (R && vista === 'deformata') s += `<g>${deformata()}</g>`;
  if (R && $('#chk-reazioni').checked && vista !== 'modello') s += `<g>${drawReazioni()}</g>`;
  if (R && $('#chk-reazioni').checked && vista === 'modello') s += `<g opacity="0.9">${drawReazioni()}</g>`;
  s += `<g id="mech">${drawMeccanismo(fase)}</g>`;
  s += `<g>${drawPreview()}</g>`;
  svg.innerHTML = s;
  $('#empty').hidden = model.elementi.length > 0 || TOOL_TIPO[tool] === 'draw';
  $('#scala-wrap').hidden = vista !== 'deformata';
  $('#btn-undo').disabled = !undoStack.length;
  $('#btn-redo').disabled = !redoStack.length;
  const needAnim = risultato && risultato.meccanismo;
  if (needAnim && !animId) {
    const loop = () => { fase += 0.06; const g = $('#mech'); if (g) g.innerHTML = drawMeccanismo(fase); animId = risultato && risultato.meccanismo ? requestAnimationFrame(loop) : null; };
    animId = requestAnimationFrame(loop);
  }
}

// ---------------------------------------------------------------- pannello proprietà
function field(label, value, onset, opts = {}) {
  const id = 'f' + Math.random().toString(36).slice(2, 8);
  const t = opts.type || 'number';
  setTimeout(() => {
    const inp = document.getElementById(id);
    if (!inp) return;
    inp.addEventListener('change', () => onset(t === 'number' ? num(inp.value, value) : inp.value));
  });
  const step = opts.step ? ` step="${opts.step}"` : ' step="any"';
  return `<div class="f"><label for="${id}">${label}</label><input id="${id}" type="${t === 'number' ? 'text' : t}" inputmode="decimal" value="${esc(t === 'number' ? fmtInput(value) : value)}"${step}${opts.ro ? ' readonly' : ''}></div>`;
}
function fmtInput(v) {
  if (!isFinite(v)) return '';
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e6 || a < 1e-4)) return v.toExponential().replace('e+', 'e');
  return String(+v.toPrecision(10));
}
function check(label, value, onset) {
  const id = 'c' + Math.random().toString(36).slice(2, 8);
  setTimeout(() => { const i = document.getElementById(id); if (i) i.addEventListener('change', () => onset(i.checked)); });
  return `<label class="check"><input type="checkbox" id="${id}" ${value ? 'checked' : ''}> ${label}</label>`;
}
function buttons(items, current, onset) {
  const id = 'b' + Math.random().toString(36).slice(2, 8);
  setTimeout(() => {
    const box = document.getElementById(id);
    if (box) box.querySelectorAll('button').forEach((b, k) => b.addEventListener('click', () => onset(items[k][0])));
  });
  return `<div class="dirs" id="${id}">${items.map(([v, lab, tip]) => `<button type="button" class="${v === current ? 'on' : ''}" title="${esc(tip || '')}">${lab}</button>`).join('')}</div>`;
}
function select(label, items, current, onset) {
  const id = 's' + Math.random().toString(36).slice(2, 8);
  setTimeout(() => { const s = document.getElementById(id); if (s) s.addEventListener('change', () => onset(s.value)); });
  return `<div class="f"><label for="${id}">${label}</label><select id="${id}">${items.map(([v, l]) => `<option value="${v}" ${v === current ? 'selected' : ''}>${l}</option>`).join('')}</select></div>`;
}
const DIR_FORZA = [[-90, '↓', 'verso il basso'], [90, '↑', "verso l'alto"], [180, '←', 'verso sinistra'], [0, '→', 'verso destra']];
// il ridisegno è rimandato di un tick: così, se l'utente preme Tab, il focus è già sul campo successivo
function edit(fn) { return (v) => { pushUndo(); fn(v); setTimeout(changed, 0); }; }

function formForza(c, target) {
  // target: oggetto con Fx, Fy (carico esistente) o default {F, ang}
  const F = target ? Math.hypot(c.Fx, c.Fy) : c.F;
  const ang = target ? Math.round(Math.atan2(c.Fy, c.Fx) * 180 / Math.PI * 1e6) / 1e6 : c.ang;
  const setFA = (F2, a2) => {
    if (target) { const r = a2 * Math.PI / 180; c.Fx = r2(F2 * Math.cos(r)); c.Fy = r2(F2 * Math.sin(r)); } else { c.F = F2; c.ang = a2; }
  };
  const wrap = target ? edit : (f) => (v) => { f(v); setTimeout(renderPanel, 0); };
  return `<div class="row">${field('Intensità [kN]', F, wrap((v) => setFA(Math.abs(v), ang)))}${field('Direzione [°]', ang, wrap((v) => setFA(F, v)))}</div>`
    + buttons(DIR_FORZA, ang, wrap((v) => setFA(F, v)));
}
function formCoppia(c, target) {
  const M = target ? Math.abs(c.M) : c.M, verso = target ? Math.sign(c.M) || 1 : c.verso;
  const wrap = target ? edit : (f) => (v) => { f(v); setTimeout(renderPanel, 0); };
  const set = (M2, v2) => { if (target) c.M = M2 * v2; else { c.M = M2; c.verso = v2; } };
  return `<div class="row one">${field('Intensità [kN·m]', M, wrap((v) => set(Math.abs(v), verso)))}</div>`
    + buttons([[1, '↺ antioraria', 'antioraria (positiva)'], [-1, '↻ oraria', 'oraria']], verso, wrap((v) => set(M, v)));
}
function formDistribuito(c, target) {
  const qi = target ? Math.abs(c.qi) : c.qi, qj = target ? Math.abs(c.qj) : c.qj;
  const verso = target ? (Math.sign(c.qi || c.qj) || -1) : c.verso;
  const wrap = target ? edit : (f) => (v) => { f(v); setTimeout(renderPanel, 0); };
  const set = (a, b, vs) => { if (target) { c.qi = a * vs; c.qj = b * vs; } else { c.qi = a; c.qj = b; c.verso = vs; } };
  const vers = {
    y: [[-1, '↓ giù'], [1, '↑ su']],
    x: [[1, '→ destra'], [-1, '← sinistra']],
    perp: [[-1, '⟂ lato destro', 'guardando la trave dal nodo i al nodo j'], [1, '⟂ lato sinistro', 'guardando la trave dal nodo i al nodo j']],
    assiale: [[1, 'i → j'], [-1, 'j → i']],
  }[c.direzione] || [];
  return `<div class="row">${field('q iniziale [kN/m]', qi, wrap((v) => set(Math.abs(v), qj, verso)))}${field('q finale [kN/m]', qj, wrap((v) => set(qi, Math.abs(v), verso)))}</div>`
    + `<p class="muted">Valori uguali: uniforme. Uno a zero: triangolare. Diversi: trapezio.</p>`
    + select('Direzione', [['y', 'Verticale'], ['x', 'Orizzontale'], ['perp', 'Perpendicolare alla trave'], ['assiale', "Lungo l'asse"]], c.direzione, wrap((v) => { c.direzione = v; }))
    + `<div style="height:6px"></div>` + buttons(vers, verso, wrap((v) => set(qi, qj, v)))
    + (c.direzione === 'y' || c.direzione === 'x' ? check('Riferito alla proiezione (es. neve su falda)', c.proiezione, wrap((v) => { c.proiezione = v; })) : '');
}
function formTermico(c, target) {
  const wrap = target ? edit : (f) => (v) => { f(v); setTimeout(renderPanel, 0); };
  return `<div class="row">${field('ΔT uniforme [°C]', c.dT, wrap((v) => { c.dT = v; }))}${field('ΔT a farfalla [°C]', c.dT_farfalla, wrap((v) => { c.dT_farfalla = v; }))}</div>`
    + `<p class="muted">Farfalla = T lato inferiore − T lato superiore (inferiore = a destra guardando da i verso j).</p>`
    + `<div class="row">${field('α [1/°C]', c.alpha, wrap((v) => { c.alpha = v; }))}${field('Altezza sezione h [m]', c.h, wrap((v) => { c.h = v > 0 ? v : c.h; }))}</div>`;
}

function renderPanel() {
  const box = $('#props');
  const focusSel = 'input,select,button';
  const act = document.activeElement;
  const focusIdx = box.contains(act) ? Array.from(box.querySelectorAll(focusSel)).indexOf(act) : -1;
  let h = '';
  if (sel && sel.k === 'nodo' && model.nodi[sel.i]) {
    const n = model.nodi[sel.i];
    h += `<h3>Nodo ${sel.i + 1}<button class="del" data-act="del">Elimina</button></h3>`;
    h += `<div class="row">${field('x [m]', n.x, edit((v) => { n.x = v; }))}${field('y [m]', n.y, edit((v) => { n.y = v; }))}</div>`;
    h += check('Cerniera interna (le travi ruotano liberamente nel nodo)', n.cerniera, edit((v) => { n.cerniera = v; }));
    const iv = model.vincoli.findIndex((v) => v.nodo === sel.i);
    h += iv >= 0 ? `<div class="list-item"><span>Vincolo: <a href="#" data-sel="vinc:${iv}">${VINCOLI[model.vincoli[iv].tipo]?.nome || model.vincoli[iv].tipo}</a></span><button data-delvinc="${iv}" title="Elimina il vincolo, non il nodo">Togli vincolo</button></div>` : '<p class="muted">Nessun vincolo. Scegli un vincolo a sinistra e clicca il nodo.</p>';
    const cn = model.carichi.map((c, k) => [c, k]).filter(([c]) => c.nodo === sel.i);
    if (cn.length) h += `<div class="sub">Carichi nel nodo</div>` + cn.map(([c, k]) => `<div class="list-item"><a href="#" data-sel="car:${k}">${descrCarico(c)}</a><button data-delcar="${k}">×</button></div>`).join('');
    if (risultato && risultato.ok) {
      const U = risultato.U.slice(3 * sel.i, 3 * sel.i + 3);
      h += `<div class="sub">Spostamenti</div><table><tr><th>ux [m]</th><th>uy [m]</th><th>φ [rad]</th></tr><tr><td>${fmt(U[0])}</td><td>${fmt(U[1])}</td><td>${fmt(U[2])}</td></tr></table>`;
    }
  } else if (sel && sel.k === 'elem' && model.elementi[sel.i]) {
    const e = model.elementi[sel.i], g = elemGeo(e);
    h += `<h3>${e.tipo === 'asta' ? 'Asta' : 'Trave'} ${sel.i + 1} <small class="muted">nodi ${e.i + 1} → ${e.j + 1}</small><button class="del" data-act="del">Elimina</button></h3>`;
    h += `<p class="muted">L = ${fmt(g.L)} m · inclinazione ${fmt(g.th * 180 / Math.PI)}°</p>`;
    h += select('Tipo', [['trave', 'Trave (N, V, M)'], ['asta', 'Asta reticolare (solo N)']], e.tipo, edit((v) => { e.tipo = v; if (v === 'asta') model.carichi = model.carichi.filter((c) => c.elem !== sel.i || (c.tipo === 'termico' && !c.dT_farfalla) || (c.tipo === 'distribuito' && c.direzione === 'assiale')); }));
    h += `<div style="height:8px"></div>`;
    if (e.tipo !== 'asta') {
      h += `<div class="row">${field('EI [kN·m²]', e.EI, edit((v) => { if (v > 0) e.EI = v; }))}${field('EA [kN]', e.EA, edit((v) => { if (v > 0) e.EA = v; }), { ro: e.inestensibile })}</div>`;
      h += check('Inestensibile (EA → ∞, come negli esercizi)', e.inestensibile, edit((v) => { e.inestensibile = v; }));
      h += check(`Cerniera all'estremo nodo ${e.i + 1}`, e.cerniera_i, edit((v) => { e.cerniera_i = v; }));
      h += check(`Cerniera all'estremo nodo ${e.j + 1}`, e.cerniera_j, edit((v) => { e.cerniera_j = v; }));
    } else {
      h += `<div class="row one">${field('EA [kN]', e.EA, edit((v) => { if (v > 0) e.EA = v; }))}</div>`;
    }
    const car = model.carichi.map((c, k) => [c, k]).filter(([c]) => c.elem === sel.i);
    h += `<div class="sub">Carichi su questa trave</div>`;
    h += car.length ? car.map(([c, k]) => `<div class="list-item"><a href="#" data-sel="car:${k}">${descrCarico(c)}</a><button data-delcar="${k}">×</button></div>`).join('') : '<p class="muted">Nessuno.</p>';
    if (risultato && risultato.ok && risultato.elementi[sel.i]) {
      const es = risultato.elementi[sel.i].estremi;
      h += `<div class="sub">Sollecitazioni</div><table><tr><th></th><th>min</th><th>max</th></tr>`
        + ['N', 'V', 'M'].map((q) => `<tr><td>${q}</td><td>${fmt(es[q].min)} <small class="muted">x=${fmt(es[q].x_min)}</small></td><td>${fmt(es[q].max)} <small class="muted">x=${fmt(es[q].x_max)}</small></td></tr>`).join('') + '</table>';
    }
  } else if (sel && sel.k === 'vinc' && model.vincoli[sel.i]) {
    const v = model.vincoli[sel.i], info = VINCOLI[v.tipo] || VINCOLI.personalizzato;
    h += `<h3>${info.nome} · nodo ${v.nodo + 1}<button class="del" data-act="del">Elimina</button></h3>`;
    h += `<div class="desc">Blocca: <b>${info.blocca}</b>.<br>Reazioni: ${info.reazioni}.</div>`;
    h += select('Tipo', Object.entries(VINCOLI).filter(([k]) => k !== 'personalizzato' || v.tipo === 'personalizzato').map(([k, o]) => [k, o.nome]), v.tipo, edit((t) => { v.tipo = t; if (t === 'molla' && !v.kx && !v.ky && !v.kphi) v.ky = 1e3; }));
    h += `<div style="height:8px"></div>`;
    h += `<div class="row">${field('Rotazione [°]', v.angolo || 0, edit((a) => { v.angolo = a; }))}<div></div></div>`;
    h += buttons([[0, '0°'], [90, '90°'], [180, '180°'], [-90, '−90°']], v.angolo || 0, edit((a) => { v.angolo = a; }));
    if (v.tipo === 'carrello' || v.tipo === 'doppio_pendolo') h += `<p class="muted">L'angolo è l'inclinazione del piano di scorrimento. Clicca di nuovo sul nodo con lo stesso strumento per ruotare di 90°.</p>`;
    if (v.tipo === 'molla') {
      h += `<div class="row three">${field('kx [kN/m]', v.kx, edit((a) => { v.kx = Math.max(0, a); }))}${field('ky [kN/m]', v.ky, edit((a) => { v.ky = Math.max(0, a); }))}${field('kφ [kN·m/rad]', v.kphi, edit((a) => { v.kphi = Math.max(0, a); }))}</div>`;
    } else {
      v.cedimento = v.cedimento || { x: 0, y: 0, phi: 0 };
      h += `<div class="sub">Cedimento vincolare (spostamento imposto)</div><div class="row three">${field('lungo x\' [m]', v.cedimento.x, edit((a) => { v.cedimento.x = a; }))}${field('lungo y\' [m]', v.cedimento.y, edit((a) => { v.cedimento.y = a; }))}${field('rotazione [rad]', v.cedimento.phi, edit((a) => { v.cedimento.phi = a; }))}</div><p class="muted">Vale solo per le direzioni bloccate. Con rotazione 0°, x' e y' coincidono con x e y.</p>`;
    }
    if (risultato && risultato.ok) {
      const r = risultato.reazioni.find((rr) => rr.nodo === v.nodo);
      if (r) h += `<div class="sub">Reazioni</div><table><tr><th>Rx [kN]</th><th>Ry [kN]</th><th>M [kN·m]</th></tr><tr><td>${fmt(r.Rx)}</td><td>${fmt(r.Ry)}</td><td>${fmt(r.M)}</td></tr></table>`;
    }
  } else if (sel && sel.k === 'car' && model.carichi[sel.i]) {
    const c = model.carichi[sel.i];
    h += `<h3>${descrCarico(c, true)}<button class="del" data-act="del">Elimina</button></h3>`;
    if (c.tipo === 'forza' || c.tipo === 'forza_campata') h += formForza(c, true);
    if (c.tipo === 'coppia' || c.tipo === 'coppia_campata') h += formCoppia(c, true);
    if (c.tipo === 'distribuito') h += formDistribuito(c, true);
    if (c.tipo === 'termico') h += formTermico(c, true);
    if ('a' in c) {
      const L = elemGeo(model.elementi[c.elem]).L;
      h += `<div class="row">${field(`Distanza dal nodo ${model.elementi[c.elem].i + 1} [m]`, c.a, edit((v) => { c.a = Math.max(0, Math.min(L, v)); }))}<div class="f"><label>Lunghezza trave</label><input value="${fmt(L)} m" readonly></div></div>`;
    }
  } else {
    h += panelTool();
  }
  box.innerHTML = h;
  if (focusIdx >= 0) {
    const f = box.querySelectorAll(focusSel)[focusIdx];
    if (f) { f.focus(); if (f.select && f.tagName === 'INPUT') f.select(); }
  }
  box.querySelectorAll('[data-act=del]').forEach((b) => b.addEventListener('click', deleteSel));
  box.querySelectorAll('[data-sel]').forEach((a) => a.addEventListener('click', (ev) => {
    ev.preventDefault(); const [k, i] = a.dataset.sel.split(':'); sel = { k, i: +i }; renderAll();
  }));
  box.querySelectorAll('[data-delvinc]').forEach((b) => b.addEventListener('click', () => {
    pushUndo(); model.vincoli.splice(+b.dataset.delvinc, 1); changed();
  }));
  box.querySelectorAll('[data-delcar]').forEach((b) => b.addEventListener('click', () => {
    pushUndo(); model.carichi.splice(+b.dataset.delcar, 1); changed();
  }));
}
function descrCarico(c, titolo) {
  switch (c.tipo) {
    case 'forza': return `Forza ${fmt(Math.hypot(c.Fx, c.Fy))} kN · nodo ${c.nodo + 1}`;
    case 'coppia': return `Coppia ${fmt(Math.abs(c.M))} kN·m ${c.M > 0 ? '↺' : '↻'} · nodo ${c.nodo + 1}`;
    case 'distribuito': return `Distribuito ${fmt(Math.abs(c.qi))}${c.qi !== c.qj ? '→' + fmt(Math.abs(c.qj)) : ''} kN/m${titolo ? ' · trave ' + (c.elem + 1) : ''}`;
    case 'forza_campata': return `Forza ${fmt(Math.hypot(c.Fx, c.Fy))} kN a ${fmt(c.a)} m${titolo ? ' · trave ' + (c.elem + 1) : ''}`;
    case 'coppia_campata': return `Coppia ${fmt(Math.abs(c.M))} kN·m a ${fmt(c.a)} m${titolo ? ' · trave ' + (c.elem + 1) : ''}`;
    case 'termico': return `Termico${titolo ? ' · trave ' + (c.elem + 1) : ''}`;
  }
  return c.tipo;
}
function panelTool() {
  const nomi = {
    seleziona: 'Seleziona', trave: 'Trave', asta: 'Asta reticolare', cerniera_interna: 'Cerniera interna',
    forza: 'Forza nel nodo', coppia: 'Coppia nel nodo', distribuito: 'Carico distribuito',
    forza_campata: 'Forza in campata', coppia_campata: 'Coppia in campata', termico: 'Variazione termica',
  };
  let h = `<h3>${nomi[tool] || VINCOLI[tool]?.nome || tool}</h3>`;
  const tt = TOOL_TIPO[tool];
  if (tool === 'seleziona') {
    h += `<p class="muted">Clicca un nodo, una trave, un vincolo o un carico per modificarlo. Trascina i nodi per spostarli. Canc elimina.</p>`;
    h += `<div class="sub">Modello</div><p class="muted">${model.nodi.length} nodi · ${model.elementi.length} travi · ${model.vincoli.length} vincoli · ${model.carichi.length} carichi</p>`;
  } else if (tool === 'trave') {
    h += `<p class="muted">Clicca il punto iniziale e poi i successivi. <b>Esc</b> o tasto destro per finire. <b>Shift</b> blocca l'angolo a multipli di 45°, <b>Alt</b> disattiva l'aggancio alla griglia.</p>`;
    h += `<div class="sub">Proprietà delle nuove travi</div>`;
    h += `<div class="row">${field('EI [kN·m²]', def.trave.EI, (v) => { if (v > 0) def.trave.EI = v; })}${field('EA [kN]', def.trave.EA, (v) => { if (v > 0) def.trave.EA = v; })}</div>`;
    h += check('Inestensibili (EA → ∞)', def.trave.inestensibile, (v) => { def.trave.inestensibile = v; });
    h += coordEntry();
  } else if (tool === 'asta') {
    h += `<p class="muted">Aste reticolari: lavorano solo a sforzo normale e sono incernierate agli estremi.</p>`;
    h += `<div class="row one">${field('EA [kN]', def.asta.EA, (v) => { if (v > 0) def.asta.EA = v; })}</div>`;
    h += coordEntry();
  } else if (tt === 'vinc') {
    const info = VINCOLI[tool];
    h += `<div class="desc">Blocca: <b>${info.blocca}</b>.<br>Reazioni: ${info.reazioni}.</div><p class="muted">Clicca un nodo per applicarlo. Clicca di nuovo sullo stesso nodo per ruotarlo di 90°.</p>`;
  } else if (tool === 'cerniera_interna') {
    h += `<p class="muted">Clicca un nodo per inserire (o togliere) una cerniera interna: le travi che vi arrivano restano collegate negli spostamenti ma ruotano liberamente. Per svincolare una sola trave, selezionala e spunta "cerniera all'estremo".</p>`;
  } else {
    h += `<p class="muted">Imposta il carico, poi clicca ${tt === 'nodeload' ? 'un nodo' : 'una trave'}.${tool === 'forza_campata' || tool === 'coppia_campata' ? ' La posizione si aggancia al passo della griglia.' : ''}</p>`;
    const d = def[tool];
    if (tool === 'forza' || tool === 'forza_campata') h += formForza(d, false);
    if (tool === 'coppia' || tool === 'coppia_campata') h += formCoppia(d, false);
    if (tool === 'distribuito') h += formDistribuito(d, false);
    if (tool === 'termico') h += formTermico(d, false);
  }
  return h;
}
function coordEntry() {
  setTimeout(() => {
    const b = $('#coord-add');
    if (!b) return;
    b.addEventListener('click', () => {
      const v = ['#c-x1', '#c-y1', '#c-x2', '#c-y2'].map((s) => num($(s).value, NaN));
      if (v.some((x) => !isFinite(x))) return flash('Inserisci quattro coordinate numeriche');
      if (Math.hypot(v[2] - v[0], v[3] - v[1]) < 1e-9) return flash('I due punti coincidono');
      pushUndo();
      const i = realize({ x: v[0], y: v[1] }), j = realize({ x: v[2], y: v[3] });
      addElem(i, j, tool === 'asta');
      $('#c-x1').value = v[2]; $('#c-y1').value = v[3];
      changed(); fit();
    });
  });
  return `<div class="sub">Oppure per coordinate</div><div class="row"><div class="f"><label>da x [m]</label><input id="c-x1" value="0"></div><div class="f"><label>da y [m]</label><input id="c-y1" value="0"></div></div>`
    + `<div class="row"><div class="f"><label>a x [m]</label><input id="c-x2" value="6"></div><div class="f"><label>a y [m]</label><input id="c-y2" value="0"></div></div><button id="coord-add">Aggiungi trave</button>`;
}

// ---------------------------------------------------------------- pannello risultati
function renderRisultati() {
  const box = $('#risultati');
  const R = risultato;
  let h = '<h3>Risultati</h3>';
  if (!R) { box.innerHTML = h + '<p class="muted">Calcolo…</p>'; return; }
  if (R.vuoto) { box.innerHTML = h + `<p class="muted">${R.msg}</p>`; return; }
  if (!R.ok) {
    if (R.labile) {
      h += `<span class="badge err">Struttura labile</span><p>${esc(R.errore)}</p>`;
      if (R.meccanismo) h += `<p class="muted">In rosso il cinematismo: il movimento che i vincoli non impediscono.</p>`;
    } else h += `<span class="badge err">Errore</span><p>${esc(R.errore)}</p>`;
    box.innerHTML = h; return;
  }
  const i = R.grado_iperstaticita;
  h += i === 0 ? '<span class="badge ok">Isostatica</span>' : `<span class="badge warn">Iperstatica · grado ${i}</span>`;
  h += R.equilibrio.ok ? ' <span class="badge ok" title="Somma di forze e momenti di carichi e reazioni">Equilibrio ✓</span>' : ' <span class="badge err">Equilibrio ✗</span>';
  h += `<div class="sub">Reazioni vincolari</div><table><tr><th>Nodo</th><th>Rx [kN]</th><th>Ry [kN]</th><th>M [kN·m]</th></tr>`;
  R.reazioni.forEach((r) => { h += `<tr><td>${r.nodo + 1} <small class="muted">${VINCOLI[r.tipo]?.nome || ''}</small></td><td>${fmt(r.Rx)}</td><td>${fmt(r.Ry)}</td><td>${fmt(r.M)}</td></tr>`; });
  h += '</table><p class="muted">Positive verso destra, verso l\'alto, antiorarie.</p>';
  h += `<div class="sub">Sollecitazioni massime</div><table><tr><th>Trave</th><th>M [kN·m]</th><th>V [kN]</th><th>N [kN]</th></tr>`;
  R.elementi.forEach((e, k) => {
    const s = e.estremi;
    const am = Math.abs(s.M.max) >= Math.abs(s.M.min) ? s.M.max : s.M.min;
    const av = Math.abs(s.V.max) >= Math.abs(s.V.min) ? s.V.max : s.V.min;
    const an = Math.abs(s.N.max) >= Math.abs(s.N.min) ? s.N.max : s.N.min;
    h += `<tr><td><a href="#" data-selelem="${k}">${k + 1}</a></td><td>${fmt(am)}</td><td>${fmt(av)}</td><td>${fmt(an)}</td></tr>`;
  });
  if (model.elementi.some((e) => e.tipo === 'asta')) h += '</table><p class="muted">Aste nella vista N: <b style="color:var(--N)">blu</b> = tiranti, <b style="color:var(--err)">rosso</b> = puntoni.</p><table>';
  h += '</table><p class="muted">Per ogni trave il valore con modulo massimo. N &gt; 0 trazione · M &gt; 0 tende le fibre inferiori (disegnato dal lato teso).</p>';
  h += `<details><summary class="sub" style="cursor:pointer">Spostamenti dei nodi</summary><table><tr><th>Nodo</th><th>ux [mm]</th><th>uy [mm]</th><th>φ [mrad]</th></tr>`;
  model.nodi.forEach((n, k) => { h += `<tr><td>${k + 1}</td><td>${fmt(R.U[3 * k] * 1e3)}</td><td>${fmt(R.U[3 * k + 1] * 1e3)}</td><td>${fmt(R.U[3 * k + 2] * 1e3)}</td></tr>`; });
  h += '</table></details>';
  if (vista === 'deformata') h += `<p class="muted">Deformata amplificata ×${fmt(scalaDeformata())}</p>`;
  box.innerHTML = h;
  box.querySelectorAll('[data-selelem]').forEach((a) => a.addEventListener('click', (ev) => { ev.preventDefault(); sel = { k: 'elem', i: +a.dataset.selelem }; renderAll(); }));
}

// ---------------------------------------------------------------- suggerimenti
function flash(msg) {
  flashMsg = msg; renderHint();
  clearTimeout(flashTimer);
  flashTimer = setTimeout(() => { flashMsg = ''; renderHint(); }, 2200);
}
function renderHint() {
  let m = flashMsg;
  if (!m) {
    const tt = TOOL_TIPO[tool];
    if (tt === 'draw') m = chainStart ? 'Clicca il punto successivo · Esc o tasto destro per finire' : 'Clicca il punto iniziale della trave';
    else if (tt === 'vinc') m = `${VINCOLI[tool].nome}: clicca un nodo (di nuovo per ruotarlo)`;
    else if (tt === 'hinge') m = 'Clicca un nodo per inserire o togliere la cerniera interna';
    else if (tt === 'nodeload') m = 'Clicca un nodo per applicare il carico';
    else if (tt === 'beamload') m = 'Clicca una trave per applicare il carico';
  }
  $('#hint').textContent = m;
}

// ---------------------------------------------------------------- vista
function fit() {
  if (!model.nodi.length) { view = { cx: 5, cy: 2, s: Math.min(W() / 14, H() / 9) }; renderAll(); return; }
  const xs = model.nodi.map((n) => n.x), ys = model.nodi.map((n) => n.y);
  const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
  const dx = Math.max(x1 - x0, 1), dy = Math.max(y1 - y0, 1);
  view.cx = (x0 + x1) / 2; view.cy = (y0 + y1) / 2;
  view.s = Math.max(5, Math.min((W() - 260) / dx, (H() - 250) / dy, 200));
  renderAll();
}

// ---------------------------------------------------------------- eventi canvas
function setTool(t) {
  tool = t; chainStart = null;
  if (t !== 'seleziona') sel = null;
  $$('.tool').forEach((b) => b.classList.toggle('on', b.dataset.tool === t));
  svg.style.cursor = t === 'seleziona' ? 'default' : 'crosshair';
  renderAll();
}
function targetOf(ev) {
  const g = ev.target.closest && ev.target.closest('[data-k]');
  return g ? { k: g.dataset.k, i: +g.dataset.i } : null;
}
function pos(ev) { const r = svg.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; }

svg.addEventListener('pointerdown', (ev) => {
  svg.focus();
  const [sx, sy] = pos(ev);
  if (ev.button === 1 || (ev.button === 0 && ev.getModifierState && ev.getModifierState(' '))) {
    drag = { tipo: 'pan', sx, sy, cx: view.cx, cy: view.cy }; svg.setPointerCapture(ev.pointerId); return;
  }
  if (ev.button === 2) { if (chainStart) { chainStart = null; renderAll(); } return; }
  const tt = TOOL_TIPO[tool];
  if (tt === 'sel') {
    // il nodo vince solo se si clicca proprio sul pallino; altrimenti conta il simbolo
    const vicino = nodeNear(sx, sy, 6);
    const dom = targetOf(ev);
    const k = vicino >= 0 ? vicino : (dom ? -1 : nodeNear(sx, sy, 10));
    const t = k >= 0 ? { k: 'nodo', i: k } : dom;
    if (t) {
      sel = t;
      if (t.k === 'nodo') drag = { tipo: 'nodo', i: k, sx, sy, moved: false };
      svg.setPointerCapture(ev.pointerId);
      renderAll();
    } else {
      drag = { tipo: 'pan', sx, sy, cx: view.cx, cy: view.cy, click: true };
      svg.setPointerCapture(ev.pointerId);
    }
    return;
  }
  if (tt === 'draw') {
    const p = snapPoint(sx, sy, ev);
    if (!chainStart) { chainStart = { x: p.x, y: p.y }; renderAll(); return; }
    if (Math.hypot(p.x - chainStart.x, p.y - chainStart.y) < 1e-9) return;
    pushUndo();
    const i = realize(chainStart), j = realize(p);
    if (addElem(i, j, tool === 'asta') < 0) flash('Questa trave esiste già');
    chainStart = { x: model.nodi[j].x, y: model.nodi[j].y };
    changed();
    return;
  }
  if (tt === 'vinc' || tt === 'hinge' || tt === 'nodeload') {
    const k = nodeNear(sx, sy, 16);
    if (k < 0) { flash('Clicca su un nodo (i punti alle estremità delle travi)'); return; }
    if (tt === 'vinc') setVincolo(k, tool);
    else if (tt === 'hinge') { pushUndo(); model.nodi[k].cerniera = !model.nodi[k].cerniera; flash(model.nodi[k].cerniera ? `Cerniera interna nel nodo ${k + 1}` : `Cerniera tolta dal nodo ${k + 1}`); changed(); }
    else addCarico(caricoDaTool(tool, k));
    return;
  }
  if (tt === 'beamload') {
    const b = beamTarget(sx, sy);
    if (!b) { flash('Clicca su una trave'); return; }
    const el = model.elementi[b.e];
    if (el.tipo === 'asta' && tool !== 'termico') { flash("Le aste accettano solo carichi nei nodi: usa una trave o aggiungi un nodo"); return; }
    if (el.tipo === 'asta' && tool === 'termico' && def.termico.dT_farfalla) { flash("Su un'asta vale solo il ΔT uniforme"); return; }
    addCarico(caricoDaTool(tool, b));
  }
});
svg.addEventListener('pointermove', (ev) => {
  const [sx, sy] = pos(ev);
  mouse = { sx, sy, x: WX(sx), y: WY(sy) };
  const tt = TOOL_TIPO[tool];
  if (drag) {
    if (drag.tipo === 'pan') {
      if (Math.hypot(sx - drag.sx, sy - drag.sy) > 3) drag.click = false;
      view.cx = drag.cx - (sx - drag.sx) / view.s; view.cy = drag.cy + (sy - drag.sy) / view.s;
      render(); return;
    }
    if (drag.tipo === 'nodo') {
      if (!drag.moved && Math.hypot(sx - drag.sx, sy - drag.sy) < 4) return;
      if (!drag.moved) { pushUndo(); drag.moved = true; }
      let x = WX(sx), y = WY(sy);
      if (!ev.altKey) { x = Math.round(x / griglia) * griglia; y = Math.round(y / griglia) * griglia; }
      const other = nodeAtWorld(x, y);
      if (other >= 0 && other !== drag.i) return;
      model.nodi[drag.i].x = r2(x); model.nodi[drag.i].y = r2(y);
      render(); schedule(); return;
    }
  }
  if (tt === 'draw') mouse.snap = snapPoint(sx, sy, ev);
  if (tt === 'beamload') mouse.beam = beamTarget(sx, sy);
  let hv = null;
  if (tt === 'vinc' || tt === 'hinge' || tt === 'nodeload' || tt === 'sel') { const k = nodeNear(sx, sy, tt === 'sel' ? 10 : 16); if (k >= 0) hv = { k: 'nodo', i: k }; }
  if (!hv && (tt === 'beamload' || tt === 'sel')) { const b = elemNear(sx, sy, 10); if (b) hv = { k: 'elem', i: b.e }; }
  hover = hv;
  $('#coords').textContent = `x = ${fmt(r2(mouse.snap ? mouse.snap.x : mouse.x))}  y = ${fmt(r2(mouse.snap ? mouse.snap.y : mouse.y))} m`;
  render();
});
svg.addEventListener('pointerup', (ev) => {
  if (drag && drag.tipo === 'pan' && drag.click) { sel = null; renderAll(); }
  if (drag && drag.tipo === 'nodo' && drag.moved) changed();
  drag = null;
});
svg.addEventListener('pointerleave', () => { mouse = null; hover = null; render(); });
svg.addEventListener('contextmenu', (ev) => ev.preventDefault());
svg.addEventListener('wheel', (ev) => {
  ev.preventDefault();
  const [sx, sy] = pos(ev);
  const wx = WX(sx), wy = WY(sy);
  const f = Math.exp(-ev.deltaY * 0.0015);
  view.s = Math.max(3, Math.min(2000, view.s * f));
  view.cx = wx - (sx - W() / 2) / view.s; view.cy = wy + (sy - H() / 2) / view.s;
  render();
}, { passive: false });
window.addEventListener('resize', () => render());

document.addEventListener('keydown', (ev) => {
  const typing = ['INPUT', 'SELECT', 'TEXTAREA'].includes(document.activeElement.tagName);
  if ((ev.ctrlKey || ev.metaKey) && ev.key.toLowerCase() === 'z') { if (typing) return; ev.preventDefault(); ev.shiftKey ? redo() : undo(); return; }
  if ((ev.ctrlKey || ev.metaKey) && ev.key.toLowerCase() === 'y') { if (typing) return; ev.preventDefault(); redo(); return; }
  if ((ev.ctrlKey || ev.metaKey) && ev.key.toLowerCase() === 's') { ev.preventDefault(); salva(); return; }
  if (typing || ev.ctrlKey || ev.metaKey || ev.altKey) return;
  if (ev.key === 'Escape') {
    if (chainStart) chainStart = null; else if (sel) sel = null; else setTool('seleziona');
    renderAll(); return;
  }
  if (ev.key === 'Delete' || ev.key === 'Backspace') { ev.preventDefault(); deleteSel(); return; }
  if (ev.key === '0') { fit(); return; }
  const b = $$('.tool').find((x) => x.dataset.key && x.dataset.key.toLowerCase() === ev.key.toLowerCase());
  if (b) setTool(b.dataset.tool);
});

// ---------------------------------------------------------------- barra superiore
$$('.tool').forEach((b) => b.addEventListener('click', () => setTool(b.dataset.tool)));
$$('#viste button').forEach((b) => b.addEventListener('click', () => {
  vista = b.dataset.vista;
  $$('#viste button').forEach((x) => x.classList.toggle('on', x === b));
  renderAll();
}));
$('#chk-reazioni').addEventListener('change', render);
$('#scala').addEventListener('input', () => { render(); renderRisultati(); });
$('#griglia').addEventListener('change', (e) => { griglia = num(e.target.value, 0.5); render(); });
$('#btn-fit').addEventListener('click', fit);
$('#btn-undo').addEventListener('click', undo);
$('#btn-redo').addEventListener('click', redo);
$('#btn-nuovo').addEventListener('click', () => {
  if (model.elementi.length && !confirm('Cancellare il modello attuale?')) return;
  pushUndo(); model = vuoto(); sel = null; chainStart = null; changed(); fit(); setTool('trave');
});
$('#btn-aiuto').addEventListener('click', () => $('#aiuto').showModal());
$('#empty-trave').addEventListener('click', () => setTool('trave'));
$('#empty-esempio').addEventListener('click', (ev) => { ev.stopPropagation(); $('#menu-esempi').classList.add('open'); });
$('#btn-esempi').addEventListener('click', (ev) => { ev.stopPropagation(); $('#menu-esempi').classList.toggle('open'); });
document.addEventListener('click', () => $('#menu-esempi').classList.remove('open'));

function download(nome, blob) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = nome;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}
function salva() {
  download('struttura.json', new Blob([JSON.stringify(model, null, 1)], { type: 'application/json' }));
}
$('#btn-salva').addEventListener('click', salva);
$('#btn-apri').addEventListener('click', () => $('#file-in').click());
$('#file-in').addEventListener('change', async (ev) => {
  const f = ev.target.files[0];
  if (!f) return;
  try {
    const m = JSON.parse(await f.text());
    if (!Array.isArray(m.nodi) || !Array.isArray(m.elementi)) throw new Error('formato non riconosciuto');
    pushUndo();
    model = { nodi: m.nodi, elementi: m.elementi, vincoli: m.vincoli || [], carichi: m.carichi || [] };
    normalizza(); sel = null; changed(); fit();
    flash(`Aperto ${f.name}`);
  } catch (e) { flash(`Impossibile aprire il file: ${e.message}`); }
  ev.target.value = '';
});
$('#btn-csv').addEventListener('click', async () => {
  if (!risultato || !risultato.ok) { flash('Prima serve una soluzione valida'); return; }
  const res = await fetch('api/csv', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(model) });
  download('risultati.csv', await res.blob());
});
$('#btn-png').addEventListener('click', () => {
  const w = W(), h = H();
  const copia = svg.cloneNode(true);
  const orig = svg.querySelectorAll('*'), dest = copia.querySelectorAll('*');
  const props = ['stroke', 'fill', 'stroke-width', 'stroke-dasharray', 'font-size', 'font-weight', 'font-family', 'opacity', 'paint-order', 'stroke-linecap'];
  orig.forEach((o, k) => {
    const cs = getComputedStyle(o);
    props.forEach((p) => dest[k].setAttribute(p, cs.getPropertyValue(p)));
  });
  copia.setAttribute('xmlns', NS);
  copia.setAttribute('width', w); copia.setAttribute('height', h);
  copia.insertAdjacentHTML('afterbegin', `<rect width="${w}" height="${h}" fill="${getComputedStyle(document.querySelector('.stage')).backgroundColor}"/>`);
  const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(copia)], { type: 'image/svg+xml' }));
  const img = new Image();
  img.onload = () => {
    const c = document.createElement('canvas'); c.width = 2 * w; c.height = 2 * h;
    const ctx = c.getContext('2d'); ctx.scale(2, 2); ctx.drawImage(img, 0, 0);
    c.toBlob((b) => download(`struttura-${vista}.png`, b));
    URL.revokeObjectURL(url);
  };
  img.src = url;
});

function normalizza() {
  model.nodi.forEach((n) => { n.x = +n.x; n.y = +n.y; n.cerniera = !!n.cerniera; });
  model.elementi.forEach((e) => { e.tipo = e.tipo || 'trave'; e.EI = e.EI || 1e4; e.EA = e.EA || 1e6; e.cerniera_i = !!e.cerniera_i; e.cerniera_j = !!e.cerniera_j; });
  model.vincoli.forEach((v) => { v.angolo = v.angolo || 0; v.kx = v.kx || 0; v.ky = v.ky || 0; v.kphi = v.kphi || 0; v.cedimento = v.cedimento || { x: 0, y: 0, phi: 0 }; });
  model.carichi.forEach((c) => {
    if (c.tipo === 'forza') { c.Fx = c.Fx || 0; c.Fy = c.Fy || 0; }
    if (c.tipo === 'distribuito') { c.qj = c.qj ?? c.qi; c.direzione = c.direzione || 'y'; c.proiezione = !!c.proiezione; }
    if (c.tipo === 'forza_campata') { c.Fx = c.Fx || 0; c.Fy = c.Fy || 0; }
  });
}

// ---------------------------------------------------------------- esempi
function costruisci(nodi, travi, vincoli, carichi, opt = {}) {
  const m = vuoto();
  m.nodi = nodi.map(([x, y, cer]) => ({ x, y, cerniera: !!cer }));
  m.elementi = travi.map(([i, j, tipo]) => ({ i, j, tipo: tipo || 'trave', EI: 1e4, EA: 1e6, inestensibile: tipo !== 'asta', cerniera_i: false, cerniera_j: false }));
  m.vincoli = vincoli.map(([nodo, tipo, angolo]) => ({ nodo, tipo, angolo: angolo || 0, kx: 0, ky: tipo === 'molla' ? 1e3 : 0, kphi: 0, cedimento: { x: 0, y: 0, phi: 0 } }));
  m.carichi = carichi;
  return m;
}
const ESEMPI = [
  ['Trave appoggiata', 'cerniera + carrello, carico uniforme', () => costruisci([[0, 0], [6, 0]], [[0, 1]], [[0, 'cerniera'], [1, 'carrello']], [{ tipo: 'distribuito', elem: 0, qi: -10, qj: -10, direzione: 'y', proiezione: false }])],
  ['Mensola', 'incastro, forza in punta e carico triangolare', () => costruisci([[0, 0], [4, 0]], [[0, 1]], [[0, 'incastro', -90]], [{ tipo: 'forza', nodo: 1, Fx: 0, Fy: -10 }, { tipo: 'distribuito', elem: 0, qi: -8, qj: 0, direzione: 'y', proiezione: false }])],
  ['Trave Gerber', 'incastro, cerniera interna, carrello', () => costruisci([[0, 0], [3, 0, true], [6, 0]], [[0, 1], [1, 2]], [[0, 'incastro', -90], [2, 'carrello']], [{ tipo: 'distribuito', elem: 1, qi: -10, qj: -10, direzione: 'y', proiezione: false }])],
  ['Trave continua', 'tre appoggi, iperstatica di grado 1', () => costruisci([[0, 0], [4, 0], [8, 0]], [[0, 1], [1, 2]], [[0, 'cerniera'], [1, 'carrello'], [2, 'carrello']], [{ tipo: 'distribuito', elem: 0, qi: -10, qj: -10, direzione: 'y', proiezione: false }, { tipo: 'distribuito', elem: 1, qi: -10, qj: -10, direzione: 'y', proiezione: false }])],
  ['Portale incastrato', 'telaio iperstatico con forza orizzontale', () => costruisci([[0, 0], [0, 4], [6, 4], [6, 0]], [[0, 1], [1, 2], [2, 3]], [[0, 'incastro'], [3, 'incastro']], [{ tipo: 'distribuito', elem: 1, qi: -20, qj: -20, direzione: 'y', proiezione: false }, { tipo: 'forza', nodo: 1, Fx: 15, Fy: 0 }])],
  ['Arco a tre cerniere', 'carico sulla proiezione orizzontale', () => costruisci([[0, 0], [4, 2, true], [8, 0]], [[0, 1], [1, 2]], [[0, 'cerniera'], [2, 'cerniera']], [{ tipo: 'distribuito', elem: 0, qi: -10, qj: -10, direzione: 'y', proiezione: true }, { tipo: 'distribuito', elem: 1, qi: -10, qj: -10, direzione: 'y', proiezione: true }])],
  ['Capriata Pratt', 'travatura reticolare', () => costruisci([[0, 0], [2, 0], [4, 0], [6, 0], [1, 2], [3, 2], [5, 2]], [[0, 1, 'asta'], [1, 2, 'asta'], [2, 3, 'asta'], [4, 5, 'asta'], [5, 6, 'asta'], [0, 4, 'asta'], [1, 4, 'asta'], [1, 5, 'asta'], [2, 5, 'asta'], [2, 6, 'asta'], [3, 6, 'asta']], [[0, 'cerniera'], [3, 'carrello']], [4, 5, 6].map((n) => ({ tipo: 'forza', nodo: n, Fx: 0, Fy: -20 })))],
  ['Carrello inclinato', 'forza in campata, appoggio su piano a 30°', () => costruisci([[0, 0], [6, 0]], [[0, 1]], [[0, 'cerniera'], [1, 'carrello', 30]], [{ tipo: 'forza_campata', elem: 0, a: 2, Fx: 0, Fy: -20 }])],
  ['Doppio pendolo', 'mensola con bipendolo in punta', () => costruisci([[0, 0], [4, 0]], [[0, 1]], [[0, 'incastro', -90], [1, 'doppio_pendolo', 90]], [{ tipo: 'forza', nodo: 1, Fx: 0, Fy: -10 }])],
];
$('#menu-esempi').innerHTML = ESEMPI.map(([n, d], k) => `<button data-es="${k}">${n}<small>${d}</small></button>`).join('');
$$('#menu-esempi [data-es]').forEach((b) => b.addEventListener('click', () => {
  pushUndo();
  model = ESEMPI[+b.dataset.es][2]();
  sel = null; setTool('seleziona'); changed(); fit();
}));

// ---------------------------------------------------------------- avvio
(function init() {
  try {
    const s = localStorage.getItem('fem-travi-modello');
    if (s) { const m = JSON.parse(s); if (m && Array.isArray(m.nodi)) { model = m; normalizza(); } }
  } catch (e) { /* storage non disponibile */ }
  setTool(model.elementi.length ? 'seleziona' : 'seleziona');
  fit();
  schedule();
})();
