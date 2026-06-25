#!/usr/bin/env python3
"""
generate_hgt_tree.py  —  Interactive HGT × phylogenetic-tree HTML visualiser

Reads:
  --tree          Newick tree file  (e.g. tree.txt)
  --ebr-overlaps  hgt_overlapping_EBRs.txt
  --block-overlaps hgt_overlapping_ancestoral_genes.txt
  --output        destination HTML file

Leaf nodes show HGT genes where EBR_species = "genus_<leaf_name>".
Internal '@' node (ancestor genome) shows genes from ancestoral_genes overlap.
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path


def parse_ebr_overlaps(path):
    """Return {short_species: [gene, ...]} from hgt_overlapping_EBRs.txt.
    EBR_species column is e.g. 'genus_sps2' → key becomes 'sps2'.
    """
    by_species = defaultdict(set)
    try:
        with open(path) as fh:
            next(fh)  # skip header
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                gene = f[3]
                ebr_species = f[7]
                short = ebr_species.split("_", 1)[-1] if "_" in ebr_species else ebr_species
                by_species[short].add(gene)
    except FileNotFoundError:
        print(f"[WARN] EBR overlaps file not found: {path}", file=sys.stderr)
    return {k: sorted(v) for k, v in by_species.items()}


def parse_block_overlaps(path):
    """Return sorted list of unique HGT genes from hgt_overlapping_ancestoral_genes.txt."""
    genes = set()
    try:
        with open(path) as fh:
            next(fh)  # skip header
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4:
                    genes.add(f[3])
    except FileNotFoundError:
        print(f"[WARN] Block-list overlaps file not found: {path}", file=sys.stderr)
    return sorted(genes)


def read_tree(path):
    return Path(path).read_text().strip()


# ── HTML template ─────────────────────────────────────────────────────────────
# __DATA_JSON__ is replaced at generation time with the embedded JSON payload.
HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>HGT Events</title>
<style>
/* ── CSS variables: dark (default) ──────────────────────────────────────────── */
:root {
  --bg:#0d1117; --bg-card:#161b27; --bg-panel:#111827; --bg-input:#0d1117;
  --bg-row-hover:#0f1f38; --bg-row-anc:#1c1000; --bg-gene:#0f2744;
  --bg-gene-anc:#2d1a00; --bg-stat:#0f2744; --bg-btn:#0f2744; --bg-btn-hover:#1e3a5f;
  --border:#21283a; --border-panel:#1e3a5f; --border-anc:#78350f;
  --border-gene:#1e3a5f; --border-input:#21283a;
  --text:#e2e8f0; --text-muted:#4a5568; --text-leaf:#94a3b8;
  --branch:#2d3e55; --root-dot:#2d3e55; --bl-label:#374151;
  --accent-blue:#3b82f6; --accent-amber:#f59e0b;
  --circle-fill:#3b82f6; --circle-empty:#1e293b;
  --circle-stroke:#1e3a5f; --circle-empty-stroke:#374151;
  --gene-color:#7dd3fc; --gene-anc:#fde68a;
  --panel-title:#60a5fa; --panel-title-anc:#fbbf24;
  --panel-sep:#1e2d45; --panel-sep2:#1a2535;
  --idx-color:#374151; --copy-color:#374151; --no-gene:#374151;
  --badge-color:#7dd3fc; --badge-anc:#fde68a;
  --shadow:rgba(0,0,0,0.7); --overlay:rgba(0,0,0,0.50);
  --bar-bg:#161b27; --bar-border:#21283a;
  --zoom-text:#64748b; --zoom-btn-bg:#1e293b; --zoom-btn-hover:#2d3e55;
  --theme-btn-bg:#1e293b; --theme-btn-border:#2d3e55;
  --layout-active-bg:#3b82f6; --layout-active-text:#fff;
}
[data-theme="light"] {
  --bg:#f1f5f9; --bg-card:#ffffff; --bg-panel:#ffffff; --bg-input:#f8fafc;
  --bg-row-hover:#eff6ff; --bg-row-anc:#fffbeb; --bg-gene:#dbeafe;
  --bg-gene-anc:#fef3c7; --bg-stat:#dbeafe; --bg-btn:#dbeafe; --bg-btn-hover:#bfdbfe;
  --border:#e2e8f0; --border-panel:#bfdbfe; --border-anc:#fde68a;
  --border-gene:#bfdbfe; --border-input:#cbd5e1;
  --text:#1e293b; --text-muted:#64748b; --text-leaf:#475569;
  --branch:#94a3b8; --root-dot:#94a3b8; --bl-label:#94a3b8;
  --accent-blue:#2563eb; --accent-amber:#d97706;
  --circle-fill:#2563eb; --circle-empty:#e2e8f0;
  --circle-stroke:#bfdbfe; --circle-empty-stroke:#94a3b8;
  --gene-color:#1d4ed8; --gene-anc:#92400e;
  --panel-title:#1d4ed8; --panel-title-anc:#92400e;
  --panel-sep:#e2e8f0; --panel-sep2:#f1f5f9;
  --idx-color:#94a3b8; --copy-color:#94a3b8; --no-gene:#94a3b8;
  --badge-color:#1d4ed8; --badge-anc:#92400e;
  --shadow:rgba(0,0,0,0.15); --overlay:rgba(0,0,0,0.25);
  --bar-bg:#f8fafc; --bar-border:#e2e8f0;
  --zoom-text:#475569; --zoom-btn-bg:#e2e8f0; --zoom-btn-hover:#cbd5e1;
  --theme-btn-bg:#e2e8f0; --theme-btn-border:#cbd5e1;
  --layout-active-bg:#2563eb; --layout-active-text:#fff;
}

* { box-sizing:border-box; margin:0; padding:0; }

body {
  font-family:'Segoe UI',system-ui,sans-serif;
  background:var(--bg); color:var(--text);
  min-height:100vh; display:flex; flex-direction:column;
  align-items:center; padding:28px 16px 60px;
  transition:background .25s,color .25s;
}

/* header */
#header {
  width:100%; max-width:1200px;
  display:flex; align-items:center; justify-content:space-between;
  margin-bottom:22px; flex-wrap:wrap; gap:10px;
}
#header-left h1 { font-size:1.45rem; font-weight:800; letter-spacing:.05em; color:var(--accent-blue); }
#header-left .subtitle { font-size:.76rem; color:var(--text-muted); margin-top:3px; }

#theme-btn {
  background:var(--theme-btn-bg); border:1px solid var(--theme-btn-border);
  border-radius:99px; color:var(--text-muted); cursor:pointer;
  padding:7px 14px; font-size:.78rem; display:flex; align-items:center;
  gap:7px; transition:background .15s,color .15s; white-space:nowrap;
}
#theme-btn:hover { color:var(--text); }

/* card */
#tree-card {
  width:100%; max-width:1200px;
  background:var(--bg-card); border:1px solid var(--border);
  border-radius:16px; overflow:hidden;
  transition:background .25s,border-color .25s;
}

/* toolbars share a row style */
.toolbar {
  display:flex; align-items:center; gap:6px;
  padding:9px 14px; border-bottom:1px solid var(--bar-border);
  background:var(--bar-bg); flex-wrap:wrap;
  transition:background .25s,border-color .25s;
}
.toolbar-label { font-size:.7rem; color:var(--text-muted); margin-right:4px; white-space:nowrap; }

/* toolbar buttons (zoom + layout share this) */
.tb-btn {
  background:var(--zoom-btn-bg); border:1px solid var(--bar-border);
  border-radius:7px; color:var(--zoom-text); cursor:pointer;
  font-size:.8rem; padding:4px 11px; line-height:1.4; font-weight:600;
  transition:background .1s,color .1s,border-color .1s;
}
.tb-btn:hover { background:var(--zoom-btn-hover); color:var(--text); }
.tb-btn.active {
  background:var(--layout-active-bg); border-color:var(--layout-active-bg);
  color:var(--layout-active-text);
}

#zoom-pct { font-size:.74rem; color:var(--zoom-text); min-width:44px; text-align:center; font-family:monospace; }
.tb-sep { width:1px; height:18px; background:var(--bar-border); margin:0 4px; }
#zoom-hint { margin-left:auto; font-size:.67rem; color:var(--text-muted); }

/* SVG area */
#tree-wrap {
  width:100%; height:580px; overflow:hidden; cursor:grab; position:relative;
}
#tree-wrap:active { cursor:grabbing; }
#tree-svg { display:block; width:100%; height:100%; user-select:none; }

/* SVG element styles */
.branch { fill:none; stroke-width:2; stroke-linecap:round; }
.node-leaf { cursor:pointer; }
.node-leaf circle { transition:filter .12s; }
.node-leaf:hover circle  { filter:drop-shadow(0 0 7px #3b82f6cc); }
.node-leaf.active circle { filter:drop-shadow(0 0 11px #60a5faff); }
.leaf-label { font-size:13px; cursor:pointer; user-select:none; transition:fill .1s; }
.node-leaf:hover .leaf-label, .node-leaf.active .leaf-label { fill:var(--text) !important; }
.node-ancestor { cursor:pointer; }
.node-ancestor polygon { transition:filter .12s; }
.node-ancestor:hover polygon  { filter:drop-shadow(0 0 8px #f59e0bcc); }
.node-ancestor.active polygon { filter:drop-shadow(0 0 13px #fbbf24ff); }
.ancestor-label { font-size:11px; font-weight:700; letter-spacing:.04em; user-select:none; }
.bl-label { font-size:9px; font-family:monospace; }

/* legend */
.legend {
  display:flex; gap:22px; padding:11px 16px;
  border-top:1px solid var(--bar-border); font-size:.73rem;
  color:var(--text-muted); flex-wrap:wrap;
  background:var(--bar-bg);
  transition:background .25s,border-color .25s;
}
.li { display:flex; align-items:center; gap:7px; }

/* overlay */
#overlay {
  display:none; position:fixed; inset:0;
  background:var(--overlay); z-index:300; backdrop-filter:blur(2px);
}
#overlay.show { display:block; }

/* side panel */
#panel {
  position:fixed; top:0; right:0; width:370px; max-width:93vw; height:100vh;
  background:var(--bg-panel); border-left:1px solid var(--border-panel);
  box-shadow:-8px 0 40px var(--shadow); z-index:400;
  display:flex; flex-direction:column;
  transform:translateX(100%);
  transition:transform .25s cubic-bezier(.4,0,.2,1),background .25s,border-color .25s;
}
#panel.show { transform:translateX(0); }
#panel-header {
  display:flex; align-items:flex-start; justify-content:space-between;
  padding:20px 20px 14px; border-bottom:1px solid var(--panel-sep); flex-shrink:0;
}
#panel-title { font-size:1rem; font-weight:700; color:var(--panel-title); line-height:1.3; }
#panel-title.anc { color:var(--panel-title-anc); }
#panel-meta { font-size:.7rem; color:var(--text-muted); margin-top:4px; }
#panel-close {
  background:none; border:none; color:var(--text-muted); font-size:1.3rem;
  cursor:pointer; padding:2px 6px; border-radius:6px; line-height:1;
  transition:color .1s,background .1s; flex-shrink:0; margin-left:10px;
}
#panel-close:hover { color:var(--text); background:var(--panel-sep); }
#panel-stats {
  padding:12px 20px; border-bottom:1px solid var(--panel-sep2);
  flex-shrink:0; display:flex; align-items:center; gap:10px;
}
.stat-badge {
  display:inline-flex; align-items:center; gap:5px;
  background:var(--bg-stat); border:1px solid var(--border-gene);
  border-radius:99px; padding:3px 10px; font-size:.75rem;
  color:var(--badge-color); font-weight:600;
}
.stat-badge.anc { background:var(--bg-gene-anc); border-color:var(--border-anc); color:var(--badge-anc); }
#panel-search { padding:10px 20px; flex-shrink:0; border-bottom:1px solid var(--panel-sep2); }
#panel-search input {
  width:100%; background:var(--bg-input); border:1px solid var(--border-input);
  border-radius:7px; color:var(--text); font-size:.8rem;
  padding:7px 10px; outline:none; transition:border-color .15s;
}
#panel-search input:focus { border-color:var(--accent-blue); }
#panel-search input::placeholder { color:var(--text-muted); }
#panel-body {
  flex:1; overflow-y:auto; padding:14px 20px 24px;
  scrollbar-width:thin; scrollbar-color:var(--border-gene) var(--bg-panel);
}
#panel-body::-webkit-scrollbar { width:5px; }
#panel-body::-webkit-scrollbar-thumb { background:var(--border-gene); border-radius:3px; }
.gene-row {
  display:flex; align-items:center; gap:10px;
  padding:9px 10px; border-radius:8px; margin-bottom:4px;
  border:1px solid transparent;
  transition:background .1s,border-color .1s; cursor:default;
}
.gene-row:hover { background:var(--bg-row-hover); border-color:var(--border-gene); }
.gene-row.anc:hover { background:var(--bg-row-anc); border-color:var(--border-anc); }
.gene-idx { font-size:.65rem; color:var(--idx-color); width:22px; text-align:right; flex-shrink:0; font-family:monospace; }
.gene-name { font-family:monospace; font-size:.82rem; color:var(--gene-color); flex:1; word-break:break-all; }
.gene-row.anc .gene-name { color:var(--gene-anc); }
.copy-btn {
  background:none; border:none; color:var(--copy-color); cursor:pointer;
  padding:2px 4px; border-radius:4px; font-size:.7rem;
  opacity:0; transition:color .1s,opacity .1s; flex-shrink:0;
}
.gene-row:hover .copy-btn { opacity:1; }
.copy-btn:hover { color:var(--gene-color); }
.copy-btn.anc:hover { color:var(--gene-anc); }
.no-genes-msg { text-align:center; color:var(--no-gene); font-style:italic; font-size:.85rem; margin-top:30px; }
#copy-all-wrap { padding:0 20px 14px; flex-shrink:0; border-top:1px solid var(--panel-sep2); }
#copy-all-btn {
  margin-top:10px; width:100%; padding:8px;
  background:var(--bg-btn); border:1px solid var(--border-gene);
  border-radius:8px; color:var(--badge-color); font-size:.78rem;
  cursor:pointer; transition:background .1s,border-color .1s;
}
#copy-all-btn:hover { background:var(--bg-btn-hover); border-color:var(--accent-blue); }
#copy-all-btn.anc { background:var(--bg-gene-anc); border-color:var(--border-anc); color:var(--badge-anc); }
#copy-all-btn.anc:hover { border-color:var(--accent-amber); }
</style>
</head>
<body>

<div id="header">
  <div id="header-left">
    <h1>HGT Events</h1>
    <div class="subtitle">Click a node to explore overlapping HGT genes</div>
  </div>
  <button id="theme-btn" onclick="toggleTheme()">
    <span id="theme-icon">☀️</span><span id="theme-label">Light mode</span>
  </button>
</div>

<div id="tree-card">

  <!-- Zoom toolbar -->
  <div class="toolbar">
    <button class="tb-btn" id="btn-zoom-in"    title="Zoom in">＋</button>
    <button class="tb-btn" id="btn-zoom-out"   title="Zoom out">－</button>
    <span   id="zoom-pct">100%</span>
    <div class="tb-sep"></div>
    <button class="tb-btn" id="btn-zoom-reset" title="Reset view">&#x2316; Reset</button>
    <span id="zoom-hint">Scroll to zoom &nbsp;·&nbsp; Drag to pan &nbsp;·&nbsp; Click node to explore</span>
  </div>

  <!-- Layout toolbar -->
  <div class="toolbar">
    <span class="toolbar-label">Tree layout:</span>
    <button class="tb-btn active" data-layout="rectangular" onclick="switchLayout('rectangular',this)">Rectangular</button>
    <button class="tb-btn"        data-layout="cladogram"   onclick="switchLayout('cladogram',this)">Cladogram</button>
    <button class="tb-btn"        data-layout="slanted"     onclick="switchLayout('slanted',this)">Slanted</button>
    <button class="tb-btn"        data-layout="circular"    onclick="switchLayout('circular',this)">Circular</button>
    <button class="tb-btn"        data-layout="unrooted"    onclick="switchLayout('unrooted',this)">Unrooted</button>
  </div>

  <!-- Tree canvas -->
  <div id="tree-wrap">
    <svg id="tree-svg"><g id="tree-g"></g></svg>
  </div>

  <!-- Legend -->
  <div class="legend">
    <div class="li">
      <svg width="14" height="14"><circle cx="7" cy="7" r="6" fill="#3b82f6" stroke="#1e3a5f" stroke-width="1.5"/></svg>
      Leaf — EBR overlaps
    </div>
    <div class="li">
      <svg width="14" height="14"><polygon points="7,0 14,14 0,14" fill="#f59e0b"/></svg>
      Ancestor node — ancestoral genes overlaps
    </div>
    <div class="li">
      <svg width="14" height="14"><circle cx="7" cy="7" r="6" fill="#1e293b" stroke="#374151" stroke-width="1.5"/></svg>
      No overlapping genes
    </div>
  </div>
</div>

<div id="overlay"></div>

<div id="panel">
  <div id="panel-header">
    <div><div id="panel-title"></div><div id="panel-meta"></div></div>
    <button id="panel-close" title="Close (Esc)">&#x2715;</button>
  </div>
  <div id="panel-stats"></div>
  <div id="panel-search">
    <input type="text" id="search-box" placeholder="Filter genes…" autocomplete="off" spellcheck="false">
  </div>
  <div id="panel-body"></div>
  <div id="copy-all-wrap">
    <button id="copy-all-btn">&#x2398; Copy all gene names</button>
  </div>
</div>

<script>
const DATA = __DATA_JSON__;

// ── Theme ─────────────────────────────────────────────────────────────────────
function toggleTheme() {
  const html = document.documentElement;
  const dark  = html.getAttribute("data-theme") === "dark";
  html.setAttribute("data-theme", dark ? "light" : "dark");
  document.getElementById("theme-icon").textContent  = dark ? "🌙" : "☀️";
  document.getElementById("theme-label").textContent = dark ? "Dark mode" : "Light mode";
}

// ── Newick parser ─────────────────────────────────────────────────────────────
function parseNewick(s) {
  s = s.trim().replace(/;$/, "").trim();
  let pos = 0;
  const eat = ch => { if (s[pos] === ch) pos++; };
  const readToken = () => { const i = pos; while (pos < s.length && !/[(),:;]/.test(s[pos])) pos++; return s.slice(i, pos).trim(); };
  const readFloat = () => { const i = pos; while (pos < s.length && /[0-9.eE+\-]/.test(s[pos])) pos++; return parseFloat(s.slice(i, pos)) || 0; };
  function node() {
    const n = { name:"", children:[], branchLength:0 };
    if (s[pos] === "(") { eat("("); n.children.push(node()); while (s[pos] === ",") { eat(","); n.children.push(node()); } eat(")"); n.name = readToken(); }
    else { n.name = readToken(); }
    if (s[pos] === ":") { eat(":"); n.branchLength = readFloat(); }
    return n;
  }
  return node();
}

// ── Tree measurements (done once) ────────────────────────────────────────────
let ROOT, N_LEAVES, TOTAL_DEPTH;

function measureTree(root) {
  let li = 0;
  function idx(n) { if (!n.children.length) n.leafIdx = li++; else n.children.forEach(idx); }
  idx(root);
  N_LEAVES = li;

  function depth(n, d) { n.depth = d; n.children.forEach(c => depth(c, d + c.branchLength)); }
  depth(root, 0);
  function maxD(n) { return n.children.length ? Math.max(...n.children.map(maxD)) : n.depth; }
  TOTAL_DEPTH = maxD(root);

  // fy: fractional leaf-y (same for all linear layouts)
  function fy(n) { if (!n.children.length) { n.fy = n.leafIdx; return; } n.children.forEach(fy); n.fy = n.children.reduce((s,c)=>s+c.fy,0)/n.children.length; }
  fy(root);

  // leafCount: for unrooted equal-angle
  function lc(n) { if (!n.children.length) { n.leafCount=1; return; } n.children.forEach(lc); n.leafCount=n.children.reduce((s,c)=>s+c.leafCount,0); }
  lc(root);

  // cladogram level
  function lvl(n,l){ n.level=l; n.children.forEach(c=>lvl(c,l+1)); }
  lvl(root,0);
  function maxLvl(n){ return n.children.length?Math.max(...n.children.map(maxLvl)):n.level; }
  n_maxLevel = maxLvl(root);
}

let n_maxLevel = 0;

// ── Zoom / pan ────────────────────────────────────────────────────────────────
let SVG_W=800, SVG_H=600, treeGroup, activeGroup=null;
const Z={k:1,tx:0,ty:0}, MIN_K=0.1, MAX_K=14;
let dragging=false, dragMoved=false, dragStart={x:0,y:0}, dragOrigin={tx:0,ty:0};

function applyZ(){
  if(treeGroup) treeGroup.setAttribute("transform",`translate(${Z.tx.toFixed(1)},${Z.ty.toFixed(1)}) scale(${Z.k.toFixed(5)})`);
  document.getElementById("zoom-pct").textContent = Math.round(Z.k*100)+"%";
}
function zoomTo(k,cx,cy){ k=Math.max(MIN_K,Math.min(MAX_K,k)); const d=k/Z.k; Z.tx=cx-d*(cx-Z.tx); Z.ty=cy-d*(cy-Z.ty); Z.k=k; applyZ(); }
function resetZoom(){
  const w=document.getElementById("tree-wrap"); const cw=w.clientWidth, ch=w.clientHeight;
  Z.k=Math.min(cw/SVG_W,ch/SVG_H)*0.88; Z.tx=(cw-SVG_W*Z.k)/2; Z.ty=(ch-SVG_H*Z.k)/2; applyZ();
}

const wrap = document.getElementById("tree-wrap");
wrap.addEventListener("wheel", e=>{
  e.preventDefault();
  const r=document.getElementById("tree-svg").getBoundingClientRect();
  zoomTo(Z.k*(e.deltaY<0?1.15:1/1.15), e.clientX-r.left, e.clientY-r.top);
},{passive:false});
wrap.addEventListener("mousedown", e=>{
  if(e.button!==0) return; dragging=true; dragMoved=false;
  dragStart={x:e.clientX,y:e.clientY}; dragOrigin={tx:Z.tx,ty:Z.ty}; e.preventDefault();
});
window.addEventListener("mousemove", e=>{
  if(!dragging) return;
  const dx=e.clientX-dragStart.x, dy=e.clientY-dragStart.y;
  if(Math.abs(dx)>3||Math.abs(dy)>3) dragMoved=true;
  Z.tx=dragOrigin.tx+dx; Z.ty=dragOrigin.ty+dy; applyZ();
});
window.addEventListener("mouseup", ()=>{ dragging=false; });

document.getElementById("btn-zoom-in").onclick   = ()=>{ const w=wrap; zoomTo(Z.k*1.3,w.clientWidth/2,w.clientHeight/2); };
document.getElementById("btn-zoom-out").onclick  = ()=>{ const w=wrap; zoomTo(Z.k/1.3,w.clientWidth/2,w.clientHeight/2); };
document.getElementById("btn-zoom-reset").onclick = resetZoom;

// ── Layout computation ────────────────────────────────────────────────────────
// All layout functions set n.x, n.y (and optionally n.angle, n.r) on every node.

const PAD  = {lx:60,rx:220,ty:70,by:70};   // linear layout padding
const ROW_H = 90, W_TREE = 520;

function computeRectangular(root) {
  SVG_W = PAD.lx+W_TREE+PAD.rx; SVG_H = PAD.ty+(N_LEAVES-1)*ROW_H+PAD.by;
  function a(n){ n.x=PAD.lx+(n.depth/TOTAL_DEPTH)*W_TREE; n.y=PAD.ty+n.fy*ROW_H; n.children.forEach(a); }
  a(root);
}
function computeCladogram(root) {
  SVG_W = PAD.lx+W_TREE+PAD.rx; SVG_H = PAD.ty+(N_LEAVES-1)*ROW_H+PAD.by;
  function a(n){
    n.x = PAD.lx + (n.children.length ? n.level/n_maxLevel : 1) * W_TREE;
    n.y = PAD.ty+n.fy*ROW_H; n.children.forEach(a);
  }
  a(root);
}
function computeSlanted(root) { computeRectangular(root); }  // same positions, different branch style

function computeCircular(root) {
  const S=620; SVG_W=S; SVG_H=S;
  const CX=S/2, CY=S/2, MAX_R=S*0.38;
  function a(n){
    if(!n.children.length){
      n.angle = (n.leafIdx/N_LEAVES)*2*Math.PI - Math.PI/2;
    } else {
      n.children.forEach(a);
      n.angle = n.children.reduce((s,c)=>s+c.angle,0)/n.children.length;
    }
    n.r = (n.depth/TOTAL_DEPTH)*MAX_R;
    n.x = CX+n.r*Math.cos(n.angle); n.y = CY+n.r*Math.sin(n.angle);
  }
  a(root); root.x=CX; root.y=CY; root.r=0;
}

function computeUnrooted(root) {
  const S=640; SVG_W=S; SVG_H=S;
  const CX=S/2, CY=S/2, MAX_R=S*0.40;
  root.x=CX; root.y=CY; root.angle=0;
  function a(n, startA, endA, px, py){
    let s=startA;
    n.children.forEach(c=>{
      const span=(c.leafCount/n.leafCount)*(endA-startA);
      c.angle=s+span/2;
      const len=(c.branchLength/TOTAL_DEPTH)*MAX_R;
      c.x=px+len*Math.cos(c.angle); c.y=py+len*Math.sin(c.angle);
      a(c, s, s+span, c.x, c.y);
      s+=span;
    });
  }
  a(root, -Math.PI/2, 3*Math.PI/2, CX, CY);
}

// ── Branch drawing ────────────────────────────────────────────────────────────
function mkEl(tag,attrs){
  const el=document.createElementNS("http://www.w3.org/2000/svg",tag);
  Object.entries(attrs).forEach(([k,v])=>el.setAttribute(k,v));
  return el;
}
function addBranch(x1,y1,x2,y2){
  const l=mkEl("line",{x1,y1,x2,y2,class:"branch"}); l.style.stroke="var(--branch)"; treeGroup.appendChild(l);
}
function addBL(x,y,txt){
  const t=mkEl("text",{x,y,"text-anchor":"middle",class:"bl-label"}); t.style.fill="var(--bl-label)"; t.textContent=txt; treeGroup.appendChild(t);
}

function drawBranch(par,ch,layout){
  switch(layout){
    case "rectangular":
    case "cladogram":
      addBranch(par.x,par.y, par.x,ch.y);
      addBranch(par.x,ch.y,  ch.x, ch.y);
      addBL((par.x+ch.x)/2, ch.y-7, ch.branchLength.toFixed(3));
      break;
    case "slanted":
      addBranch(par.x,par.y,ch.x,ch.y);
      addBL((par.x+ch.x)/2,(par.y+ch.y)/2-8,ch.branchLength.toFixed(3));
      break;
    case "circular":
      drawCircularBranch(par,ch);
      break;
    case "unrooted":
      addBranch(par.x,par.y,ch.x,ch.y);
      addBL((par.x+ch.x)/2,(par.y+ch.y)/2-8,ch.branchLength.toFixed(3));
      break;
  }
}

function drawCircularBranch(par,ch){
  if(par.r<0.5){ addBranch(par.x,par.y,ch.x,ch.y); return; }
  const CX=SVG_W/2, CY=SVG_H/2, r=par.r;
  const x1=CX+r*Math.cos(par.angle), y1=CY+r*Math.sin(par.angle);
  const x2=CX+r*Math.cos(ch.angle),  y2=CY+r*Math.sin(ch.angle);
  let da=ch.angle-par.angle;
  while(da> Math.PI) da-=2*Math.PI;
  while(da<-Math.PI) da+=2*Math.PI;
  const large=Math.abs(da)>Math.PI?1:0, sweep=da>0?1:0;
  const arc=mkEl("path",{d:`M${x1.toFixed(2)} ${y1.toFixed(2)} A${r.toFixed(2)} ${r.toFixed(2)} 0 ${large} ${sweep} ${x2.toFixed(2)} ${y2.toFixed(2)}`,class:"branch",fill:"none"});
  arc.style.stroke="var(--branch)"; arc.style.strokeWidth="2"; treeGroup.appendChild(arc);
  addBranch(x2,y2,ch.x,ch.y);
}

// ── Radial label helper ───────────────────────────────────────────────────────
function radialLabel(n, offset){
  const a=n.angle||0, lx=n.x+Math.cos(a)*offset, ly=n.y+Math.sin(a)*offset+4;
  const norm=((a%(2*Math.PI))+(2*Math.PI))%(2*Math.PI);
  return {x:lx, y:ly, anchor:(norm>Math.PI/2&&norm<3*Math.PI/2)?"end":"start"};
}

// ── Node drawing ──────────────────────────────────────────────────────────────
function drawAllNodes(root, layout){
  const radial = layout==="circular" || layout==="unrooted";

  function draw(n){
    const isAnc  = n.name==="@";
    const isLeaf = n.children.length===0;
    const isRoot = !n.name && !isLeaf;

    if(isRoot){
      const d=mkEl("circle",{cx:n.x,cy:n.y,r:5,class:"root-dot"});
      d.style.fill="var(--root-dot)"; treeGroup.appendChild(d);
      n.children.forEach(draw); return;
    }

    const genes    = isAnc ? DATA.ancestorGenes : (DATA.speciesGenes[n.name]||[]);
    const hasGenes = genes.length>0;
    const g = mkEl("g",{class:isAnc?"node-ancestor":"node-leaf"});

    if(isAnc){
      const s=11;
      const poly=mkEl("polygon",{points:`${n.x},${n.y-s} ${n.x+s},${n.y} ${n.x},${n.y+s} ${n.x-s},${n.y}`,stroke:"var(--border-anc)","stroke-width":1.5});
      poly.style.fill="var(--accent-amber)"; g.appendChild(poly);
      const lp = radial ? radialLabel(n,20) : {x:n.x, y:n.y-s-9, anchor:"middle"};
      const al=mkEl("text",{x:lp.x,y:lp.y,"text-anchor":lp.anchor,class:"ancestor-label"});
      al.style.fill="var(--accent-amber)"; al.textContent="Ancestor"; g.appendChild(al);

    } else {
      const circ=mkEl("circle",{cx:n.x,cy:n.y,r:8,"stroke-width":2});
      circ.style.fill   = hasGenes?"var(--circle-fill)":"var(--circle-empty)";
      circ.style.stroke = hasGenes?"var(--circle-stroke)":"var(--circle-empty-stroke)";
      g.appendChild(circ);
      const lp = radial ? radialLabel(n,16) : {x:n.x+16,y:n.y+5,anchor:"start"};
      const ll=mkEl("text",{x:lp.x,y:lp.y,"text-anchor":lp.anchor,class:"leaf-label"});
      ll.style.fill="var(--text-leaf)"; ll.textContent=n.name; g.appendChild(ll);
      if(hasGenes){
        const bx = lp.anchor==="end" ? lp.x-n.name.length*7.6-10 : lp.x+n.name.length*7.6+8;
        const badge=mkEl("text",{x:bx,y:lp.y-1,"text-anchor":lp.anchor==="end"?"end":"start"});
        badge.style.cssText="font-size:10px;font-weight:700;fill:var(--accent-blue);";
        badge.textContent=`(${genes.length})`; g.appendChild(badge);
      }
    }

    g.addEventListener("click",()=>{ if(!dragMoved) openPanel(n.name,genes,isAnc,g); });
    treeGroup.appendChild(g);
    n.children.forEach(draw);
  }
  draw(root);
}

// ── Layout switcher ───────────────────────────────────────────────────────────
let currentLayout = "rectangular";

function switchLayout(id, btnEl){
  currentLayout = id;
  document.querySelectorAll(".tb-btn[data-layout]").forEach(b=>b.classList.remove("active"));
  (btnEl || document.querySelector(`.tb-btn[data-layout="${id}"]`)).classList.add("active");

  treeGroup.innerHTML = "";
  if(activeGroup){ activeGroup=null; }

  switch(id){
    case "rectangular": computeRectangular(ROOT); break;
    case "cladogram":   computeCladogram(ROOT);   break;
    case "slanted":     computeSlanted(ROOT);      break;
    case "circular":    computeCircular(ROOT);     break;
    case "unrooted":    computeUnrooted(ROOT);     break;
  }

  // Draw branches first, then nodes on top
  function drawB(n){ n.children.forEach(c=>{ drawBranch(n,c,id); drawB(c); }); }
  drawB(ROOT);
  drawAllNodes(ROOT, id);

  requestAnimationFrame(resetZoom);
}

// ── Side panel ────────────────────────────────────────────────────────────────
const panel=document.getElementById("panel"), overlay=document.getElementById("overlay");
const ptitle=document.getElementById("panel-title"), pmeta=document.getElementById("panel-meta");
const pstats=document.getElementById("panel-stats"), pbody=document.getElementById("panel-body");
const searchBox=document.getElementById("search-box"), copyAllBtn=document.getElementById("copy-all-btn");
let currentGenes=[], currentIsAnc=false;

function openPanel(name,genes,isAnc,groupEl){
  currentGenes=genes; currentIsAnc=isAnc;
  if(activeGroup) activeGroup.classList.remove("active");
  activeGroup=groupEl; groupEl.classList.add("active");
  ptitle.textContent=isAnc?"◆ Ancestor genome node":`● ${name}`;
  ptitle.className=isAnc?"anc":"";
  pmeta.textContent=isAnc?"Source: ancestoral genes overlaps (DESCHRAMBLER)":`Source: EBR overlaps · species: genus_${name}`;
  const bc=isAnc?"stat-badge anc":"stat-badge";
  pstats.innerHTML=genes.length
    ?`<span class="${bc}">${genes.length} HGT gene${genes.length!==1?"s":""}</span>`
    :`<span class="${bc}" style="opacity:.5;">0 genes</span>`;
  copyAllBtn.className=isAnc?"anc":"";
  copyAllBtn.style.display=genes.length?"":"none";
  searchBox.value="";
  renderGeneList(genes,isAnc);
  panel.classList.add("show"); overlay.classList.add("show"); searchBox.focus();
}

function renderGeneList(genes,isAnc){
  if(!genes.length){ pbody.innerHTML=`<p class="no-genes-msg">No overlapping HGT genes found</p>`; return; }
  const rc=isAnc?"gene-row anc":"gene-row", cc=isAnc?"copy-btn anc":"copy-btn";
  pbody.innerHTML=genes.map((g,i)=>`
    <div class="${rc}">
      <span class="gene-idx">${i+1}</span>
      <span class="gene-name">${g}</span>
      <button class="${cc}" onclick="copyGene('${g}',this)" title="Copy">&#x2398;</button>
    </div>`).join("");
}

function filterGenes(q){
  q=q.trim().toLowerCase();
  renderGeneList(q?currentGenes.filter(g=>g.toLowerCase().includes(q)):currentGenes,currentIsAnc);
}
function closePanel(){
  panel.classList.remove("show"); overlay.classList.remove("show");
  if(activeGroup){ activeGroup.classList.remove("active"); activeGroup=null; }
}
function copyGene(gene,btn){
  navigator.clipboard.writeText(gene).then(()=>{ const o=btn.innerHTML; btn.innerHTML="✓"; setTimeout(()=>btn.innerHTML=o,1200); });
}
copyAllBtn.addEventListener("click",()=>{
  navigator.clipboard.writeText(currentGenes.join("\n")).then(()=>{ const o=copyAllBtn.textContent; copyAllBtn.textContent="✓ Copied!"; setTimeout(()=>copyAllBtn.textContent=o,1500); });
});
searchBox.addEventListener("input",e=>filterGenes(e.target.value));
document.getElementById("panel-close").addEventListener("click",closePanel);
overlay.addEventListener("click",closePanel);
document.addEventListener("keydown",e=>{ if(e.key==="Escape") closePanel(); });

// ── Boot ──────────────────────────────────────────────────────────────────────
treeGroup = document.getElementById("tree-g");
ROOT = parseNewick(DATA.newick);
measureTree(ROOT);
switchLayout("rectangular");
</script>
</body>
</html>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tree",           required=True, help="Newick tree file")
    ap.add_argument("--ebr-overlaps",   required=True, help="hgt_overlapping_EBRs.txt")
    ap.add_argument("--block-overlaps", required=True, help="hgt_overlapping_ancestoral_genes.txt")
    ap.add_argument("--output",         required=True, help="Output HTML file")
    args = ap.parse_args()

    newick         = read_tree(args.tree)
    species_genes  = parse_ebr_overlaps(args.ebr_overlaps)
    ancestor_genes = parse_block_overlaps(args.block_overlaps)

    data = {
        "newick":        newick,
        "speciesGenes":  species_genes,
        "ancestorGenes": ancestor_genes,
    }

    html = HTML_TEMPLATE.replace("__DATA_JSON__", json.dumps(data, indent=4))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html)

    print(
        f"\033[32;1m✔\033[0m HGT tree HTML generated → {args.output}",
        file=sys.stderr
    )
    print(f"  Leaf species with genes : {sorted(species_genes)}", file=sys.stderr)
    print(f"  Ancestor genes          : {len(ancestor_genes)}", file=sys.stderr)


if __name__ == "__main__":
    main()
