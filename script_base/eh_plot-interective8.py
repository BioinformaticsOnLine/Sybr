"""
eh_interactive.py  —  Interactive EH synteny plot  (zero external dependencies)
=================================================================================

Input files
-----------
Synteny (required, 9-col TSV, no header):
  chr  start  end  tarChr  tarSt  tarEnd  orient  ref  tar

Protein annotation (optional, 5-col TSV, no header):
  accession  start  end  chr_num  protein_id

Chromosome sizes (optional):
  --sizes 1:15000000,2:14000000,3:12000000
  --sizes chrom_sizes.tsv     (tab-delimited: chr<TAB>size)

Usage
-----
  python eh_interactive.py synteny.tsv --chr 1,2,3
  python eh_interactive.py synteny.tsv --chr 1,2,3 \\
      --annot protein_annotation.tsv --sizes 1:15000000,2:14000000
  python eh_interactive.py            # self-test
"""

import sys, json, argparse, os
import pandas as pd

SYN_COLS   = ["chr","start","end","tarChr","tarSt","tarEnd","orient","ref","tar"]
ANNOT_COLS = ["accession","start","end","chr","protein"]

def load_synteny(path):
    df = pd.read_csv(path, sep="\t", header=None, names=SYN_COLS,
                     dtype={"chr":str,"tarChr":str})
    for c in ("start","end","tarSt","tarEnd"):
        df[c] = pd.to_numeric(df[c])
    return df

def load_annot(path):
    df = pd.read_csv(path, sep="\t", header=None, names=ANNOT_COLS,
                     dtype={"chr":str})
    for c in ("start","end"):
        df[c] = pd.to_numeric(df[c])
    return df

def parse_sizes(sizes_arg):
    """Parse chromosome sizes from  'chr:size,chr:size' string or a TSV file."""
    sizes = {}
    if os.path.isfile(sizes_arg):
        with open(sizes_arg) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        sizes[str(parts[0])] = int(parts[1])
                    except ValueError:
                        pass
    else:
        for item in sizes_arg.split(","):
            item = item.strip()
            if ":" in item:
                ch, sz = item.rsplit(":", 1)
                try:
                    sizes[ch.strip()] = int(sz.strip())
                except ValueError:
                    pass
    return sizes

def build_html(syn_df, chr_range, annot_df, chrom_sizes, output):
    chr_range = [str(c) for c in chr_range]
    ref_label = str(syn_df["ref"].iloc[0])
    targets   = list(dict.fromkeys(syn_df["tar"]))

    # ── synteny data ──────────────────────────────────────────────────────
    syn_data = {}
    for chrom in chr_range:
        sub = syn_df[syn_df["chr"] == chrom]
        if sub.empty: continue
        # use provided size if given, else max block end
        data_end  = int(sub["end"].max())
        chr_size  = chrom_sizes.get(chrom, data_end)
        tar_chrs  = {}
        for tar in targets:
            tchrs = sorted(sub[sub["tar"]==tar]["tarChr"].unique().tolist(),
                           key=lambda x: (len(x), x))
            tar_chrs[tar] = tchrs
        syn_data[chrom] = {
            "max_end"  : chr_size,      # authoritative chromosome length
            "data_end" : data_end,      # last block position (for info)
            "tar_chrs" : tar_chrs,
            "blocks"   : [
                {"chr":r["chr"],"start":int(r["start"]),"end":int(r["end"]),
                 "tarChr":str(r["tarChr"]),"tarSt":int(r["tarSt"]),"tarEnd":int(r["tarEnd"]),
                 "orient":r["orient"],"ref":r["ref"],"tar":r["tar"]}
                for _,r in sub.iterrows()
            ]
        }

    # ── annotation data ───────────────────────────────────────────────────
    annot_data = {}
    if annot_df is not None:
        for chrom in chr_range:
            sub = annot_df[annot_df["chr"] == chrom]
            if sub.empty: continue
            annot_data[chrom] = [
                {"s":int(r["start"]),"e":int(r["end"]),"p":r["protein"]}
                for _,r in sub.iterrows()
            ]

    # ── dynamic layout constants (computed in Python) ─────────────────────
    # Vertical header height: based on longest species name (text becomes
    # height when rotated −90°). ~6.8 px per char at 11 px font.
    max_name_px = max((len(t.replace("_"," ")) * 6.8 for t in targets), default=80)
    V_HDR_H_val = max(100, min(220, int(max_name_px) + 24))

    SYN_J    = json.dumps(syn_data)
    ANN_J    = json.dumps(annot_data)
    TGTS_J   = json.dumps(targets)
    CHR_J    = json.dumps(chr_range)
    REF_LBL  = ref_label.replace("_"," ")
    HAS_ANN  = "true" if annot_data else "false"
    V_HDR_H_J = V_HDR_H_val

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>EH Plot — {REF_LBL}</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:'Segoe UI',Arial,sans-serif;background:#eef0f4;color:#2c3e50;overflow:hidden;height:100vh}}

#hdr{{background:linear-gradient(135deg,#1a252f,#2c3e50);color:#fff;
  padding:8px 14px;display:flex;align-items:center;gap:10px;flex-wrap:wrap;
  position:sticky;top:0;z-index:200;box-shadow:0 2px 8px rgba(0,0,0,.5)}}
#hdr h1{{font-size:13px;font-weight:700;white-space:nowrap}}
#hdr h1 em{{color:#5dade2;font-style:normal}}
.tabs{{display:flex;gap:4px;flex-wrap:wrap}}
.tab{{padding:3px 10px;border:1.5px solid rgba(255,255,255,.35);border-radius:20px;
  cursor:pointer;font-size:12px;color:#ccc;background:transparent;transition:all .15s}}
.tab:hover{{background:rgba(255,255,255,.15);color:#fff}}
.tab.active{{background:#2980b9;border-color:#2980b9;color:#fff;font-weight:600}}
.ctrl{{display:flex;gap:6px;align-items:center;margin-left:auto;flex-wrap:wrap}}
.btn{{padding:4px 11px;border:none;border-radius:4px;cursor:pointer;font-size:12px;
  font-weight:600;background:#2980b9;color:#fff;transition:background .15s}}
.btn:hover{{background:#1f6fa0}}
.btn.toggled{{background:#27ae60}}.btn.toggled:hover{{background:#1e8449}}
#zlbl{{font-size:11px;color:#aaa;min-width:72px;white-space:nowrap}}

#leg{{display:flex;gap:16px;align-items:center;flex-wrap:wrap;
  padding:5px 14px;background:#fff;border-bottom:1px solid #ddd;font-size:12px}}
.li{{display:flex;align-items:center;gap:5px}}
.sw{{width:18px;height:11px;border:1px solid #999;border-radius:2px}}
#hint{{margin-left:auto;color:#999;font-size:11px}}

#wrap{{display:flex;height:calc(100vh - 68px);overflow:hidden;position:relative}}

/* ── HORIZONTAL ── */
/* h-lcol hidden — labels are drawn inside hSppSvg for pixel-perfect alignment */
#h-lcol{{display:none}}
#h-lref{{background:#dde4ea;border-bottom:2px solid #aab5be;
  display:flex;align-items:center;justify-content:flex-end;
  padding-right:10px;font-size:11px;font-weight:700;color:#546e7a}}
#h-lspp{{overflow:hidden}}
.h-sl{{display:flex;align-items:center;justify-content:flex-end;padding-right:10px;
  font-size:11.5px;font-style:italic;color:#37474f;box-shadow:0 1px 0 0 #e0e3e8;
  cursor:default;transition:background .1s}}
.h-sl:hover{{background:#e8edf2}}
#h-pcol{{flex:1;display:flex;flex-direction:column;overflow:hidden}}
#h-refpanel{{flex:0 0 auto;overflow:hidden;background:#dde4ea;border-bottom:2px solid #aab5be;cursor:grab}}
#h-refpanel:active{{cursor:grabbing}}
#h-spppanel{{flex:1;overflow-y:auto;overflow-x:hidden;cursor:grab}}
#h-spppanel:active{{cursor:grabbing}}

/* ── VERTICAL ── */
#v-wrap{{display:flex;flex:1;overflow:hidden}}
#v-yaxis{{flex:0 0 58px;display:flex;flex-direction:column;overflow:hidden;
  background:#eceff1;border-right:1px solid #d5d8dc}}
#v-yaxis-top{{flex:0 0 auto;background:#eceff1;border-bottom:2px solid #aab5be}}
#v-yaxis-svg-wrap{{flex:1;overflow:hidden;position:relative}}
#v-right{{flex:1;display:flex;flex-direction:column;overflow:hidden}}
#v-col-header{{flex:0 0 auto;overflow:hidden;background:#fff}}
#v-resizer{{
  flex:0 0 7px;background:#b0bec5;cursor:row-resize;
  display:flex;align-items:center;justify-content:center;
  user-select:none;transition:background .15s;border-bottom:1px solid #90a4ae
}}
#v-resizer:hover,#v-resizer.dragging{{background:#78909c}}
#v-resizer::after{{content:'';width:44px;height:2px;background:rgba(255,255,255,.65);border-radius:2px}}
#v-scroll{{flex:1;overflow:auto;cursor:grab}}
#v-scroll:active{{cursor:grabbing}}

svg{{display:block;user-select:none}}
.hidden{{display:none!important}}

/* ── CONTROL PANEL ── */
#ctrl-panel{{position:absolute;right:0;top:0;bottom:0;width:250px;background:#fff;
  border-left:2px solid #b0bec5;display:flex;flex-direction:column;z-index:150;
  box-shadow:-4px 0 16px rgba(0,0,0,.18);transform:translateX(100%);transition:transform .25s ease}}
#ctrl-panel.open{{transform:translateX(0)}}
#ctrl-hdr{{background:#2c3e50;color:#fff;padding:10px 12px;
  display:flex;justify-content:space-between;align-items:center;font-size:13px;font-weight:700}}
#ctrl-hdr button{{background:transparent;border:1px solid rgba(255,255,255,.5);
  color:#fff;border-radius:4px;padding:2px 8px;cursor:pointer;font-size:12px}}
#ctrl-body{{flex:1;overflow-y:auto;padding:8px 0}}
.spp-item{{border-bottom:1px solid #eceff1}}
.spp-row{{display:flex;align-items:center;gap:8px;padding:7px 12px;
  background:#f8f9fb;cursor:pointer;user-select:none}}
.spp-row:hover{{background:#e8ecf0}}
.spp-row label{{font-size:12px;font-style:italic;cursor:pointer;flex:1}}
.spp-arrow{{font-size:10px;color:#888;transition:transform .2s}}
.spp-arrow.open{{transform:rotate(90deg)}}
.chr-list{{padding:4px 12px 4px 32px;display:none}}
.chr-list.open{{display:block}}
.chr-row{{display:flex;align-items:center;gap:6px;padding:3px 0}}
.chr-row label{{font-size:11.5px;color:#546e7a;cursor:pointer}}
input[type=checkbox]{{width:14px;height:14px;cursor:pointer;accent-color:#2980b9}}
#ctrl-footer{{padding:8px 12px;border-top:1px solid #eceff1;display:flex;gap:6px}}
#ctrl-footer button{{flex:1;padding:5px;border:none;border-radius:4px;cursor:pointer;font-size:11.5px;font-weight:600}}
#btn-all{{background:#27ae60;color:#fff}}#btn-all:hover{{background:#1e8449}}
#btn-none{{background:#c0392b;color:#fff}}#btn-none:hover{{background:#a93226}}

/* ── TOOLTIP ── */
#tip{{position:fixed;display:none;pointer-events:none;z-index:999;
  background:rgba(15,15,25,.93);color:#ecf0f1;padding:9px 13px;border-radius:7px;
  font-size:12px;line-height:1.85;min-width:220px;max-width:310px;
  box-shadow:0 4px 20px rgba(0,0,0,.55)}}
#tip .ttl{{font-weight:700;color:#5dade2;font-size:13px;margin-bottom:3px}}
#tip .row{{display:flex;justify-content:space-between;gap:8px}}
#tip .k{{color:#85929e}}
</style>
</head>
<body>

<div id="hdr">
  <h1>🧬 Reference: <em>{REF_LBL}</em></h1>
  <div class="tabs" id="tabs"></div>
  <div class="ctrl">
    <span id="zlbl">Zoom: 1×</span>
    <button class="btn" id="btn-view"  onclick="toggleView()">⇄ Vertical View</button>
    <button class="btn" id="btn-ctrl" onclick="toggleCtrl()">☰ Controls</button>
    <button class="btn" id="btn-reset" onclick="resetZoom()">⟳ Reset</button>
  </div>
</div>

<div id="leg">
  <div class="li"><div class="sw" style="background:lightblue"></div>Co-linear (+)</div>
  <div class="li"><div class="sw" style="background:lightpink"></div>Inverted (−)</div>
  <div class="li"><div class="sw" style="background:#dde1e7"></div>No synteny</div>
  <div class="li" id="ann-leg" style="display:none">
    <div class="sw" style="background:#e67e22;border-color:#b35c00"></div>Protein annotation
  </div>
  <div id="hint">🖱 Scroll to zoom &nbsp;·&nbsp; Drag to pan &nbsp;·&nbsp; Hover for details</div>
</div>

<div id="wrap">

  <!-- HORIZONTAL -->
  <div id="h-lcol">
    <div id="h-lref"></div>
    <div id="h-lspp"></div>
  </div>
  <div id="h-pcol">
    <div id="h-refpanel"><svg id="h-refsvg"></svg></div>
    <div id="h-spppanel"><svg id="h-sppsvg"></svg></div>
  </div>

  <!-- VERTICAL -->
  <div id="v-wrap" class="hidden">
    <div id="v-yaxis">
      <div id="v-yaxis-top"></div>
      <div id="v-yaxis-svg-wrap"><svg id="v-ysvg"></svg></div>
    </div>
    <div id="v-right">
      <div id="v-col-header"><svg id="v-hdrsvg"></svg></div>
      <div id="v-resizer" title="Drag to resize species name bar"></div>
      <div id="v-scroll"><svg id="v-mainsvg"></svg></div>
    </div>
  </div>

  <!-- CONTROL PANEL -->
  <div id="ctrl-panel">
    <div id="ctrl-hdr"><span>⚙ Show / Hide</span><button onclick="toggleCtrl()">✕ Close</button></div>
    <div id="ctrl-body"></div>
    <div id="ctrl-footer">
      <button id="btn-all"  onclick="setAll(true)">✓ Show All</button>
      <button id="btn-none" onclick="setAll(false)">✕ Hide All</button>
    </div>
  </div>

</div>
<div id="tip"></div>

<script>
// ── embedded data ─────────────────────────────────────────────────────────
const SYN    = {SYN_J};
const ANNOT  = {ANN_J};
const TGTS   = {TGTS_J};
const CHRS   = {CHR_J};
const HAS_ANN= {HAS_ANN};

// ── layout ────────────────────────────────────────────────────────────────
// horizontal
const H_ROW_H=54,H_REF_H=72,H_BPAD=9,H_BLK_H=H_ROW_H-H_BPAD*2,H_ML=10,H_ANN_H=14;
// Label area drawn inside the SVG (replaces the old h-lcol div).
// HLBL = label width, HGAP = gap before bar starts — total matches old h-lcol width.
const HLBL=175, HGAP=10;

// vertical — mirrors horizontal but rotated 90°
// V_HDR_H is computed from the longest species name (Python-side)
let V_HDR_H = {V_HDR_H_J};   // flexible: auto-computed from name length, user can drag-resize
const V_REF_W  = 52;
const V_SPP_W  = 54;   // one column per species  (≈ H_ROW_H)
const V_SPP_GAP= 12;
const V_BPAD   = 8;
const V_BLK_W  = V_SPP_W - V_BPAD*2;
const V_ML     = 10;
const V_LEFT   = 10;

// Annotation bar width — grows with zoom so labels become readable
function annBarW(z){{
  if(!HAS_ANN) return 0;
  // minimum 18 px; grows logarithmically with zoom up to 54 px
  return Math.min(54, Math.max(18, Math.round(18 * Math.log2(z + 1))));
}}

// ── state ─────────────────────────────────────────────────────────────────
let curChr  = CHRS[0];
let zoom    = 1, panX=0, panY=0;
let viewMode= "horizontal";
let isDrag  = false, dragX0=0, dragY0=0, pan0X=0, pan0Y=0;
let isResizing=false, resizeStartY=0, resizeStartH=0;
let visState= {{}};

const tip  = document.getElementById("tip");
const zlbl = document.getElementById("zlbl");
if(HAS_ANN) document.getElementById("ann-leg").style.display="flex";

// ── visibility ────────────────────────────────────────────────────────────
function initVis(chr){{
  const cd=SYN[chr];
  TGTS.forEach(tar=>{{
    if(!visState[tar]) visState[tar]={{visible:true,chrs:{{}}}};
    if(cd&&cd.tar_chrs[tar])
      cd.tar_chrs[tar].forEach(tc=>{{if(visState[tar].chrs[tc]===undefined) visState[tar].chrs[tc]=true;}});
  }});
}}
const isSppVis=(tar)=>visState[tar]?.visible!==false;
const isChrVis=(tar,tc)=>visState[tar]?.chrs[tc]!==false;
const visTgts=()=>TGTS.filter(isSppVis);

// ── SVG helpers ───────────────────────────────────────────────────────────
const NS="http://www.w3.org/2000/svg";
function el(tag,attrs,txt){{
  const e=document.createElementNS(NS,tag);
  for(const[k,v] of Object.entries(attrs)) e.setAttribute(k,v);
  if(txt!==undefined) e.textContent=txt; return e;
}}
function setSZ(s,w,h){{s.setAttribute("width",w);s.setAttribute("height",h);}}
function clear(s){{while(s.firstChild)s.removeChild(s.firstChild);}}
function addClip(s,id,x,y,w,h){{
  const d=el("defs",{{}}),cp=el("clipPath",{{id}}),r=el("rect",{{x,y,width:w,height:h}});
  cp.appendChild(r);d.appendChild(cp);s.appendChild(d);
}}
function fmtBp(v){{
  if(v>=1e6) return (v/1e6).toFixed(2).replace(/\\.?0+$/,"")+' Mb';
  if(v>=1e3) return (v/1e3).toFixed(1).replace(/\\.?0+$/,"")+' kb';
  return v+' bp';
}}
function niceTicks(lo,hi,n){{
  const r=hi-lo,rough=r/n,mag=Math.pow(10,Math.floor(Math.log10(rough)));
  const norm=rough/mag,nice=norm<1.5?1:norm<3.5?2:norm<7.5?5:10;
  const step=nice*mag,first=Math.ceil(lo/step)*step,t=[];
  for(let v=first;v<=hi;v+=step) t.push(v); return t;
}}

// ── species info tooltip (shared by both views) ───────────────────────────
function showSppTip(e,tar){{
  const cd=SYN[curChr]; if(!cd) return;
  const blks=cd.blocks.filter(b=>b.tar===tar&&isSppVis(tar));
  const visBlks=blks.filter(b=>isChrVis(tar,b.tarChr));
  const totalCov=visBlks.reduce((s,b)=>s+(b.end-b.start),0);
  const maxEnd=cd.max_end;
  const pct=((totalCov/maxEnd)*100).toFixed(1);
  const uChrs=[...new Set(visBlks.map(b=>b.tarChr))];
  const inverted=visBlks.filter(b=>b.orient==="-").length;
  tip.innerHTML=`<div class='ttl'>${{tar.replace(/_/g," ")}}</div>
    <div class='row'><span class='k'>Visible blocks</span><span>${{visBlks.length}} / ${{blks.length}}</span></div>
    <div class='row'><span class='k'>Co-linear</span><span>${{visBlks.length-inverted}}</span></div>
    <div class='row'><span class='k'>Inverted</span><span>${{inverted}}</span></div>
    <div class='row'><span class='k'>Coverage</span><span>${{totalCov.toLocaleString()}} bp (${{pct}}%)</span></div>
    <div class='row'><span class='k'>Target chrs (${{uChrs.length}})</span><span>${{uChrs.join(", ")}}</span></div>`;
  tip.style.display="block"; moveTip(e);
}}
function showSynTip(e,b){{
  const ori=b.orient==="+"?"<span style='color:#7ec8e3'>co-linear (+)</span>":"<span style='color:#f1948a'>inverted (−)</span>";
  tip.innerHTML=`<div class='ttl'>${{b.tar.replace(/_/g," ")}} &middot; chr ${{b.tarChr}}</div>
    <div class='row'><span class='k'>Ref chr</span><span>chr ${{b.chr}}</span></div>
    <div class='row'><span class='k'>Ref start</span><span>${{b.start.toLocaleString()}} bp</span></div>
    <div class='row'><span class='k'>Ref end</span><span>${{b.end.toLocaleString()}} bp</span></div>
    <div class='row'><span class='k'>Target chr</span><span>chr ${{b.tarChr}}</span></div>
    <div class='row'><span class='k'>Target start</span><span>${{b.tarSt.toLocaleString()}} bp</span></div>
    <div class='row'><span class='k'>Target end</span><span>${{b.tarEnd.toLocaleString()}} bp</span></div>
    <div class='row'><span class='k'>Block length</span><span>${{(b.end-b.start).toLocaleString()}} bp</span></div>
    <div class='row'><span class='k'>Orientation</span>${{ori}}</div>`;
  tip.style.display="block"; moveTip(e);
}}
function showAnnTip(e,a){{
  tip.innerHTML=`<div class='ttl'>Protein annotation</div>
    <div class='row'><span class='k'>Protein ID</span><span>${{a.p}}</span></div>
    <div class='row'><span class='k'>Start</span><span>${{a.s.toLocaleString()}} bp</span></div>
    <div class='row'><span class='k'>End</span><span>${{a.e.toLocaleString()}} bp</span></div>
    <div class='row'><span class='k'>Length</span><span>${{(a.e-a.s).toLocaleString()}} bp</span></div>`;
  tip.style.display="block"; moveTip(e);
}}
function moveTip(e){{
  let x=e.clientX+16,y=e.clientY+14;
  if(x+315>window.innerWidth)  x=e.clientX-325;
  if(y+280>window.innerHeight) y=e.clientY-290;
  tip.style.left=x+"px"; tip.style.top=y+"px";
}}
function hideTip(){{tip.style.display="none";}}

// ════════════════════════════════════════════════════════════════════════════
//  HORIZONTAL RENDER
// ════════════════════════════════════════════════════════════════════════════
const hRefSvg  =document.getElementById("h-refsvg");
const hSppSvg  =document.getElementById("h-sppsvg");
const hRefPanel=document.getElementById("h-refpanel");
const hSppPanel=document.getElementById("h-spppanel");
// hLref / hLspp divs still exist in DOM but h-lcol is display:none;
// species labels are drawn as SVG text inside hSppSvg below.

// hG2P: converts genomic coord → x pixel; barX is the bar-area left edge
// (defined as a module-level function so onWheel can use HLBL+HGAP)
function hG2P(g,maxEnd,iW){{return HLBL+HGAP+(g/maxEnd)*iW*zoom+panX;}}

function renderHorizontal(){{
  const cd=SYN[curChr]; if(!cd) return;
  const ann=ANNOT[curChr]||[];
  const maxEnd=cd.max_end;
  const pW=document.getElementById("h-pcol").clientWidth;
  // barX = HLBL label area + HGAP gap; bars start here (same x as old h-lcol width)
  const barX=HLBL+HGAP, iW=pW-barX-H_ML, sH=TGTS.length*H_ROW_H;
  const minP=-(zoom-1)*iW; panX=Math.min(0,Math.max(minP,panX));
  const visLo=Math.max(0,(-panX/(iW*zoom))*maxEnd);
  const visHi=Math.min(maxEnd,((iW-panX)/(iW*zoom))*maxEnd);
  zlbl.textContent="Zoom: "+(zoom<10?zoom.toFixed(1):Math.round(zoom))+"\xd7";

  const annH=HAS_ANN?H_ANN_H:0;
  const refH=H_REF_H-(HAS_ANN?0:H_ANN_H), barTop=annH+2, barH=refH-barTop-18;
  setSZ(hRefSvg,pW,refH); hRefPanel.style.height=refH+"px"; clear(hRefSvg);
  addClip(hRefSvg,"h-rclip",barX,0,iW,refH);

  // "Reference" label inside ref SVG — same coordinate system as the bar
  hRefSvg.appendChild(el("text",{{
    x:barX-8,y:barTop+barH/2,"text-anchor":"end","dominant-baseline":"middle",
    "font-size":"11px","font-weight":"700",fill:"#546e7a"
  }},"chr "+curChr));

  // axis ticks
  const axG=el("g",{{}});
  axG.appendChild(el("line",{{x1:barX,y1:refH-14,x2:barX+iW,y2:refH-14,stroke:"#90a4ae","stroke-width":"1"}}));
  niceTicks(visLo,visHi,Math.max(4,Math.floor(iW/110))).forEach(t=>{{
    const px=hG2P(t,maxEnd,iW); if(px<barX||px>barX+iW) return;
    axG.appendChild(el("line",{{x1:px,y1:refH-14,x2:px,y2:refH-8,stroke:"#90a4ae","stroke-width":"1"}}));
    axG.appendChild(el("text",{{x:px,y:refH-2,"text-anchor":"middle",fill:"#546e7a","font-size":"10px"}},fmtBp(t)));
  }}); hRefSvg.appendChild(axG);

  // annotation track
  if(HAS_ANN&&ann.length){{
    const aG=el("g",{{"clip-path":"url(#h-rclip)"}});
    aG.appendChild(el("rect",{{x:barX,y:2,width:iW,height:annH-2,fill:"#fef9e7",stroke:"#f0d9a0","stroke-width":"0.5"}}));
    hRefSvg.appendChild(el("text",{{x:barX+3,y:annH-3,fill:"#a04000","font-size":"9px","font-weight":"600"}},"proteins"));
    ann.forEach(a=>{{
      const ax=hG2P(a.s,maxEnd,iW),aw=hG2P(a.e,maxEnd,iW)-ax;
      if(aw<0.3||ax+aw<barX||ax>barX+iW) return;
      const r=el("rect",{{x:ax,y:3,width:Math.max(0.6,aw),height:annH-4,
        fill:"#e67e22",opacity:"0.8","clip-path":"url(#h-rclip)",cursor:"pointer"}});
      r.addEventListener("mouseenter",e=>showAnnTip(e,a));
      r.addEventListener("mousemove",moveTip); r.addEventListener("mouseleave",hideTip);
      aG.appendChild(r);
      if(aw>40&&zoom>8){{
        const maxChars=Math.floor(aw/5.5);
        const lbl=a.p.length>maxChars?a.p.slice(0,maxChars-1)+"\u2026":a.p;
        aG.appendChild(el("text",{{x:ax+aw/2,y:annH/2+3,"text-anchor":"middle",
          "font-size":"8px",fill:"#7d4e07","pointer-events":"none","clip-path":"url(#h-rclip)"}},lbl));
      }}
    }}); hRefSvg.appendChild(aG);
  }}

  // ref bar
  const rx=hG2P(0,maxEnd,iW),rw=hG2P(maxEnd,maxEnd,iW)-rx;
  hRefSvg.appendChild(el("rect",{{x:rx,y:barTop,width:Math.max(0,rw),height:barH,fill:"#455a64",rx:3,"clip-path":"url(#h-rclip)"}}));
  if(rw>50) hRefSvg.appendChild(el("text",{{
    x:Math.max(rx,barX)+Math.min(rw,iW)/2,y:barTop+barH/2,
    "text-anchor":"middle","dominant-baseline":"middle",
    fill:"#fff","font-size":"13px","font-weight":"bold","clip-path":"url(#h-rclip)"
  }},"chr "+curChr));

  // ── species rows: labels AND bars drawn in same SVG — pixel-perfect alignment ──
  setSZ(hSppSvg,pW,sH); clear(hSppSvg);
  addClip(hSppSvg,"h-sclip",barX,0,iW,sH);

  TGTS.forEach((tar,ri)=>{{
    const y=ri*H_ROW_H, blks=cd.blocks.filter(b=>b.tar===tar), vis=isSppVis(tar);
    // row background (full width)
    hSppSvg.appendChild(el("rect",{{x:0,y,width:pW,height:H_ROW_H,
      fill:ri%2===0?"#f8f9fb":"#fff",stroke:"#e4e7ea","stroke-width":".5",opacity:vis?1:0.35}}));
    // species name — drawn at exact same y as chromosome bar (same SVG, same coords)
    hSppSvg.appendChild(el("text",{{
      x:barX-8,y:y+H_ROW_H/2,"text-anchor":"end","dominant-baseline":"middle",
      "font-size":"11.5px","font-style":"italic",fill:vis?"#37474f":"#aaa"
    }},tar.replace(/_/g," ")));
    // invisible hover rect over label area
    const nameOver=el("rect",{{x:0,y,width:barX-4,height:H_ROW_H,fill:"transparent",cursor:"default"}});
    nameOver.addEventListener("mouseenter",e=>showSppTip(e,tar));
    nameOver.addEventListener("mousemove",moveTip); nameOver.addEventListener("mouseleave",hideTip);
    hSppSvg.appendChild(nameOver);

    if(!vis){{
      hSppSvg.appendChild(el("text",{{x:barX+4,y:y+H_ROW_H/2,"dominant-baseline":"middle",fill:"#aaa","font-size":"11px","font-style":"italic"}},"(hidden)"));
      return;
    }}
    const ex=hG2P(0,maxEnd,iW),ew=hG2P(maxEnd,maxEnd,iW)-ex;
    if(ew>0) hSppSvg.appendChild(el("rect",{{x:ex,y:y+H_BPAD,width:Math.max(0,ew),height:H_BLK_H,fill:"#dde1e7",rx:2,"clip-path":"url(#h-sclip)"}}));
    blks.forEach(b=>{{
      if(!isChrVis(tar,b.tarChr)) return;
      const bx=hG2P(b.start,maxEnd,iW),bw=hG2P(b.end,maxEnd,iW)-bx;
      if(bw<0.3||bx+bw<barX||bx>barX+iW) return;
      const fc=b.orient==="+"?"lightblue":"lightpink",ec=b.orient==="+"?"#4a90b8":"#c0607a";
      const r=el("rect",{{x:Math.max(bx,barX),y:y+H_BPAD,
        width:Math.min(bw,barX+iW-Math.max(bx,barX)),height:H_BLK_H,
        fill:fc,stroke:ec,"stroke-width":bw>4?.8:.3,rx:2,"clip-path":"url(#h-sclip)",cursor:"pointer"}});
      r.addEventListener("mouseenter",e=>showSynTip(e,b));
      r.addEventListener("mousemove",moveTip); r.addEventListener("mouseleave",hideTip);
      hSppSvg.appendChild(r);
      if(bw>22){{
        const vbx=Math.max(bx,barX),vbw=Math.min(bx+bw,barX+iW)-vbx;
        hSppSvg.appendChild(el("text",{{x:vbx+vbw/2,y:y+H_BPAD+H_BLK_H/2,
          "text-anchor":"middle","dominant-baseline":"middle","font-size":Math.min(11,bw/2.8)+"px",
          fill:"#1a252f","pointer-events":"none","clip-path":"url(#h-sclip)"}},b.tarChr));
      }}
    }});
  }});
}}

// ════════════════════════════════════════════════════════════════════════════
//  VERTICAL RENDER — horizontal layout rotated 90°
//  - Y-axis  = reference genomic coordinates
//  - Left column  = reference chromosome (vertical bar)
//  - Right columns = one vertical bar per target species
// ════════════════════════════════════════════════════════════════════════════
const vHdrSvg=document.getElementById("v-hdrsvg");
const vYSvg  =document.getElementById("v-ysvg");
const vMain  =document.getElementById("v-mainsvg");
const vScroll=document.getElementById("v-scroll");

function vG2Y(g,maxEnd,iH){{return V_ML+(g/maxEnd)*iH*zoom+panY;}}

function renderVertical(){{
  const cd=SYN[curChr]; if(!cd) return;
  const ann=ANNOT[curChr]||[];
  const maxEnd=cd.max_end;

  const wrapH=document.getElementById("wrap").clientHeight;
  const innerH=wrapH-V_ML*2;
  const totalH=V_ML+innerH*zoom+V_ML;
  panY=Math.min(0,Math.max(-(zoom-1)*innerH,panY));

  const visLo=Math.max(0,(-panY/(innerH*zoom))*maxEnd);
  const visHi=Math.min(maxEnd,((innerH-panY)/(innerH*zoom))*maxEnd);
  zlbl.textContent="Zoom: "+(zoom<10?zoom.toFixed(1):Math.round(zoom))+"×";

  const vt=visTgts();
  const annW=annBarW(zoom);  // FLEXIBLE annotation bar width
  const mainW=V_LEFT+V_REF_W+annW+vt.length*(V_SPP_W+V_SPP_GAP)+V_LEFT;

  // ── Y-axis (fixed left) ──────────────────────────────────────────────
  const yAxisW=58;
  document.getElementById("v-yaxis-top").style.height=V_HDR_H+"px";
  const yAreaH=wrapH-V_HDR_H;
  setSZ(vYSvg,yAxisW,yAreaH); clear(vYSvg);
  vYSvg.appendChild(el("line",{{x1:yAxisW-6,y1:0,x2:yAxisW-6,y2:yAreaH,stroke:"#90a4ae","stroke-width":"1"}}));
  const scrollTop=vScroll.scrollTop;
  niceTicks(visLo,visHi,Math.max(5,Math.floor(innerH/80))).forEach(t=>{{
    const py=vG2Y(t,maxEnd,innerH)-V_ML;
    const screenY=py-scrollTop;
    if(screenY<0||screenY>yAreaH) return;
    vYSvg.appendChild(el("line",{{x1:yAxisW-14,y1:screenY,x2:yAxisW-6,y2:screenY,stroke:"#90a4ae","stroke-width":"1"}}));
    vYSvg.appendChild(el("text",{{x:yAxisW-16,y:screenY,"text-anchor":"end","dominant-baseline":"middle",fill:"#546e7a","font-size":"9px"}},fmtBp(t)));
  }});
  vScroll.onscroll=()=>renderVertical();

  // ── Column header (fixed top) ────────────────────────────────────────
  setSZ(vHdrSvg,mainW,V_HDR_H); clear(vHdrSvg);
  addClip(vHdrSvg,"v-hclip",0,0,mainW,V_HDR_H);

  // reference column label
  const refLblX=V_LEFT+V_REF_W/2, refLblY=V_HDR_H/2;
  const refT=el("text",{{x:refLblX,y:refLblY,"text-anchor":"middle","dominant-baseline":"middle",
    "font-size":"11px","font-weight":"700",fill:"#455a64",
    transform:`rotate(-90,${{refLblX}},${{refLblY}})`}});
  refT.textContent="chr "+curChr+" (ref)";
  vHdrSvg.appendChild(refT);

  // annotation column label
  if(HAS_ANN){{
    const ax2=V_LEFT+V_REF_W+annW/2, ay2=V_HDR_H/2;
    vHdrSvg.appendChild(el("text",{{x:ax2,y:ay2,"text-anchor":"middle","dominant-baseline":"middle",
      "font-size":"8px",fill:"#a04000",transform:`rotate(-90,${{ax2}},${{ay2}})`}},"prot."));
  }}

  // species column labels — rotated −90°, with hover tooltip (no box, just separators)
  vt.forEach((tar,si)=>{{
    const bx=V_LEFT+V_REF_W+annW+si*(V_SPP_W+V_SPP_GAP)+V_SPP_GAP/2;
    const cx=bx+V_SPP_W/2, cy=V_HDR_H/2;
    // thin vertical separator between species (no box)
    if(si>0){{
      vHdrSvg.appendChild(el("line",{{
        x1:bx-V_SPP_GAP/2+1,y1:14,x2:bx-V_SPP_GAP/2+1,y2:V_HDR_H-14,
        stroke:"#b0bec5","stroke-width":"0.8"
      }}));
    }}
    // rotated name (truncate to fit header height)
    const full=tar.replace(/_/g," ");
    const maxChars=Math.floor((V_HDR_H-16)/6.5);
    const label=full.length>maxChars?full.slice(0,maxChars-1)+"…":full;
    const t=el("text",{{x:cx,y:cy,"text-anchor":"middle","dominant-baseline":"middle",
      "font-size":"11px","font-style":"italic",fill:"#2c3e50",
      transform:`rotate(-90,${{cx}},${{cy}})`}});
    t.textContent=label; vHdrSvg.appendChild(t);
    // invisible overlay for hover tooltip
    const over=el("rect",{{x:bx,y:0,width:V_SPP_W,height:V_HDR_H,fill:"transparent",cursor:"default"}});
    over.addEventListener("mouseenter",e=>showSppTip(e,tar));
    over.addEventListener("mousemove",moveTip);
    over.addEventListener("mouseleave",hideTip);
    vHdrSvg.appendChild(over);
  }});

  // ── Main scrollable area ─────────────────────────────────────────────
  setSZ(vMain,mainW,totalH); clear(vMain);
  addClip(vMain,"v-mclip",0,0,mainW,totalH);
  const mainG=el("g",{{"clip-path":"url(#v-mclip)"}});

  // reference bar (vertical, leftmost)
  const refTop=vG2Y(0,maxEnd,innerH), refBot=vG2Y(maxEnd,maxEnd,innerH);
  const refBarH=Math.max(0,refBot-refTop);
  mainG.appendChild(el("rect",{{x:V_LEFT,y:refTop,width:V_REF_W,height:refBarH,fill:"#455a64",rx:3}}));
  if(refBarH>30){{
    const midY=refTop+Math.min(refBarH/2,100);
    const rt=el("text",{{x:V_LEFT+V_REF_W/2,y:midY,"text-anchor":"middle","dominant-baseline":"middle",
      "font-size":"11px","font-weight":"bold",fill:"#fff",
      transform:`rotate(-90,${{V_LEFT+V_REF_W/2}},${{midY}})`}});
    rt.textContent="chr "+curChr; mainG.appendChild(rt);
  }}

  // annotation strip (right of ref bar, vertical)
  // FLEXIBLE: annW grows with zoom; labels appear when bar is tall enough
  if(HAS_ANN&&ann.length){{
    const ax=V_LEFT+V_REF_W;
    mainG.appendChild(el("rect",{{x:ax,y:refTop,width:annW,height:refBarH,
      fill:"#fef9e7",stroke:"#f0d9a0","stroke-width":"0.5"}}));
    ann.forEach(a=>{{
      const ay=vG2Y(a.s,maxEnd,innerH), ah=vG2Y(a.e,maxEnd,innerH)-ay;
      if(ah<0.3) return;
      const r=el("rect",{{x:ax,y:ay,width:annW,height:Math.max(0.6,ah),
        fill:"#e67e22",opacity:"0.8",cursor:"pointer"}});
      r.addEventListener("mouseenter",e=>showAnnTip(e,a));
      r.addEventListener("mousemove",moveTip); r.addEventListener("mouseleave",hideTip);
      mainG.appendChild(r);
      // protein ID label — rotated −90° — shown when block is tall enough
      if(ah>14&&zoom>6){{
        const px=ax+annW/2, py=ay+ah/2;
        const fs=Math.min(9,annW*0.45);
        const maxChars=Math.max(3,Math.floor(ah/(fs*0.75)));
        const lbl=a.p.length>maxChars?a.p.slice(0,maxChars-1)+"…":a.p;
        mainG.appendChild(el("text",{{x:px,y:py,"text-anchor":"middle","dominant-baseline":"middle",
          "font-size":fs+"px",fill:"#7d4e07","pointer-events":"none",
          transform:`rotate(-90,${{px}},${{py}})`}},lbl));
      }}
    }});
  }}

  // species bars — one vertical bar per species
  vt.forEach((tar,si)=>{{
    const bx=V_LEFT+V_REF_W+annW+si*(V_SPP_W+V_SPP_GAP)+V_SPP_GAP/2;
    const blks=cd.blocks.filter(b=>b.tar===tar);
    // alternating column background
    mainG.appendChild(el("rect",{{x:bx-V_SPP_GAP/2,y:0,
      width:V_SPP_W+V_SPP_GAP,height:totalH,
      fill:si%2===0?"rgba(248,249,251,0.55)":"rgba(255,255,255,0.55)",stroke:"none"}}));
    // empty bar (full chromosome height)
    mainG.appendChild(el("rect",{{x:bx+V_BPAD,y:refTop,width:V_BLK_W,height:refBarH,
      fill:"#dde1e7",rx:2,stroke:"#c5cad0","stroke-width":"0.5"}}));
    // synteny blocks at REFERENCE coordinates
    blks.forEach(b=>{{
      if(!isChrVis(tar,b.tarChr)) return;
      const by=vG2Y(b.start,maxEnd,innerH), bh=vG2Y(b.end,maxEnd,innerH)-by;
      if(bh<0.3) return;
      const fc=b.orient==="+"?"lightblue":"lightpink";
      const ec=b.orient==="+"?"#4a90b8":"#c0607a";
      const r=el("rect",{{x:bx+V_BPAD,y:by,width:V_BLK_W,height:Math.max(0.5,bh),
        fill:fc,stroke:ec,"stroke-width":bh>4?.8:.3,rx:2,cursor:"pointer"}});
      r.addEventListener("mouseenter",e=>showSynTip(e,b));
      r.addEventListener("mousemove",moveTip); r.addEventListener("mouseleave",hideTip);
      mainG.appendChild(r);
      // target chr label inside block when tall enough
      if(bh>16){{
        const fs=Math.min(10,bh/2.5,V_BLK_W*0.55);
        mainG.appendChild(el("text",{{x:bx+V_SPP_W/2,y:by+bh/2,
          "text-anchor":"middle","dominant-baseline":"middle",
          "font-size":fs+"px",fill:"#1a252f","pointer-events":"none"}},b.tarChr));
      }}
    }});
  }});

  vMain.appendChild(mainG);
}}

// ════════════════════════════════════════════════════════════════════════════
//  CONTROL PANEL
// ════════════════════════════════════════════════════════════════════════════
function buildCtrlPanel(){{
  const body=document.getElementById("ctrl-body"); body.innerHTML="";
  const cd=SYN[curChr];
  TGTS.forEach(tar=>{{
    const tchrs=cd&&cd.tar_chrs[tar]?cd.tar_chrs[tar]:[];
    const div=document.createElement("div"); div.className="spp-item";
    const srow=document.createElement("div"); srow.className="spp-row";
    const cb=document.createElement("input"); cb.type="checkbox";
    cb.checked=isSppVis(tar); cb.id="spp_"+tar;
    cb.addEventListener("change",()=>{{
      visState[tar].visible=cb.checked;
      div.querySelectorAll(".chr-cb").forEach(c=>{{c.disabled=!cb.checked;}});
      render();
    }});
    const lbl=document.createElement("label"); lbl.htmlFor="spp_"+tar;
    lbl.textContent=tar.replace(/_/g," ");
    const arrow=document.createElement("span"); arrow.className="spp-arrow"; arrow.textContent="▶";
    srow.appendChild(cb); srow.appendChild(lbl); srow.appendChild(arrow);
    const clist=document.createElement("div"); clist.className="chr-list";
    tchrs.forEach(tc=>{{
      const crow=document.createElement("div"); crow.className="chr-row";
      const ccb=document.createElement("input"); ccb.type="checkbox"; ccb.className="chr-cb";
      ccb.id="chr_"+tar+"_"+tc; ccb.disabled=!isSppVis(tar); ccb.checked=isChrVis(tar,tc);
      ccb.addEventListener("change",()=>{{visState[tar].chrs[tc]=ccb.checked;render();}});
      const clbl=document.createElement("label"); clbl.htmlFor="chr_"+tar+"_"+tc;
      clbl.textContent="chr "+tc;
      crow.appendChild(ccb); crow.appendChild(clbl); clist.appendChild(crow);
    }});
    srow.addEventListener("click",(e)=>{{
      if(e.target===cb||e.target===lbl) return;
      const open=clist.classList.toggle("open");
      arrow.classList.toggle("open",open);
    }});
    div.appendChild(srow); div.appendChild(clist); body.appendChild(div);
  }});
}}
function toggleCtrl(){{
  const p=document.getElementById("ctrl-panel");
  const open=p.classList.toggle("open");
  document.getElementById("btn-ctrl").classList.toggle("toggled",open);
  if(open) buildCtrlPanel();
}}
function setAll(val){{
  TGTS.forEach(tar=>{{
    visState[tar].visible=val;
    Object.keys(visState[tar].chrs).forEach(tc=>{{visState[tar].chrs[tc]=val;}});
  }});
  buildCtrlPanel(); render();
}}

// ════════════════════════════════════════════════════════════════════════════
//  ZOOM, PAN, TABS, VIEW TOGGLE
// ════════════════════════════════════════════════════════════════════════════
function render(){{ if(viewMode==="horizontal") renderHorizontal(); else renderVertical(); }}

function onWheel(e){{
  e.preventDefault();
  const delta=e.deltaY<0?1.15:1/1.15;
  const newZoom=Math.min(10000,Math.max(0.5,zoom*delta));
  if(viewMode==="horizontal"){{
    const rc=e.currentTarget.getBoundingClientRect();
    const mX=e.clientX-rc.left-(HLBL+HGAP); panX=mX-(mX-panX)*(newZoom/zoom);
  }} else {{
    const rc=e.currentTarget.getBoundingClientRect();
    const mY=e.clientY-rc.top+vScroll.scrollTop;
    panY=mY-(mY-panY)*(newZoom/zoom);
  }}
  zoom=newZoom; render();
}}
function onDown(e){{isDrag=true;dragX0=e.clientX;dragY0=e.clientY;pan0X=panX;pan0Y=panY;}}
function onMove(e){{
  if(isResizing){{
    const delta=e.clientY-resizeStartY;
    V_HDR_H=Math.max(40,Math.min(320,resizeStartH+delta));
    document.getElementById("v-yaxis-top").style.height=V_HDR_H+"px";
    if(viewMode==="vertical") render();
    return;
  }}
  if(!isDrag) return;
  if(viewMode==="horizontal") panX=pan0X+(e.clientX-dragX0);
  else                        panY=pan0Y+(e.clientY-dragY0);
  render();
}}
function onUp(){{isDrag=false; isResizing=false; document.getElementById("v-resizer").classList.remove("dragging");}}

[document.getElementById("h-refsvg"),
 document.getElementById("h-spppanel"),
 document.getElementById("v-scroll")].forEach(el2=>{{
  el2.addEventListener("wheel",onWheel,{{passive:false}});
  el2.addEventListener("mousedown",onDown);
}});
// drag-resize handler for species name bar
document.getElementById("v-resizer").addEventListener("mousedown",e=>{{
  isResizing=true;
  resizeStartY=e.clientY;
  resizeStartH=V_HDR_H;
  document.getElementById("v-resizer").classList.add("dragging");
  e.preventDefault(); e.stopPropagation();
}});
window.addEventListener("mousemove",onMove);
window.addEventListener("mouseup",onUp);

// chromosome tabs
const tabsDiv=document.getElementById("tabs");
CHRS.forEach(c=>{{
  const b=document.createElement("button");
  b.className="tab"+(c===curChr?" active":"");
  b.textContent="chr "+c; b.onclick=()=>selectChr(c); tabsDiv.appendChild(b);
}});

function selectChr(c){{
  curChr=c; zoom=1; panX=0; panY=0; initVis(c);
  document.querySelectorAll(".tab").forEach(b=>b.classList.toggle("active",b.textContent.trim()==="chr "+c));
  if(document.getElementById("ctrl-panel").classList.contains("open")) buildCtrlPanel();
  render();
}}
function resetZoom(){{zoom=1;panX=0;panY=0;render();}}

function toggleView(){{
  viewMode=viewMode==="horizontal"?"vertical":"horizontal";
  zoom=1;panX=0;panY=0;
  const btn=document.getElementById("btn-view");
  const hL=document.getElementById("h-lcol"),hP=document.getElementById("h-pcol");
  const vW=document.getElementById("v-wrap");
  if(viewMode==="vertical"){{
    btn.textContent="⇄ Horizontal View"; btn.classList.add("toggled");
    hL.classList.add("hidden"); hP.classList.add("hidden"); vW.classList.remove("hidden");
  }} else {{
    btn.textContent="⇄ Vertical View"; btn.classList.remove("toggled");
    hL.classList.remove("hidden"); hP.classList.remove("hidden"); vW.classList.add("hidden");
  }}
  render();
}}

// ── init ─────────────────────────────────────────────────────────────────
initVis(curChr);
render();
window.addEventListener("resize",render);
</script>
</body>
</html>"""

    out = f"{output}.html"
    with open(out, "w") as f:
        f.write(html)
    print(f"Saved \u2192 {out}")

# ── CLI ───────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description="Interactive EH synteny plot (HTML, zero dependencies)")
    p.add_argument("infile")
    p.add_argument("--chr",   dest="chr_range", required=True,
                   help='Comma-separated reference chromosomes, e.g. "1,2,3"')
    p.add_argument("--annot", default=None,
                   help="5-col protein annotation TSV (accession start end chr protein)")
    p.add_argument("--sizes", default=None,
                   help='Chromosome sizes: "1:15000000,2:14000000" or a chr<TAB>size TSV file')
    p.add_argument("--output", default="eh_interactive")
    a = p.parse_args()

    syn_df = load_synteny(a.infile)
    ann_df = load_annot(a.annot) if a.annot else None
    sizes  = parse_sizes(a.sizes) if a.sizes else {}

    build_html(syn_df, [c.strip() for c in a.chr_range.split(",")],
               ann_df, sizes, a.output)

# ── self-test ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) == 1:
        rows = [
            ["1",3304587,3309784,"2",984714,989338,"+","Adineta_vaga:100k","Brachionus_angularis"],
            ["1",6570163,15147997,"6",1037884,2557947,"+","Adineta_vaga:100k","Brachionus_angularis"],
            ["2",10520628,10531937,"8",2691072,2698511,"+","Adineta_vaga:100k","Brachionus_angularis"],
            ["2",12393409,12500425,"3",1962858,2035300,"-","Adineta_vaga:100k","Brachionus_angularis"],
            ["2",14384211,14431493,"3",3568836,3653844,"+","Adineta_vaga:100k","Brachionus_angularis"],
            ["3",6262071,11168946,"2",1106029,4395677,"+","Adineta_vaga:100k","Brachionus_angularis"],
            ["3",12313479,12347362,"1",2706760,2794912,"+","Adineta_vaga:100k","Brachionus_angularis"],
            ["4",508335,2897387,"1",1157471,2075909,"+","Adineta_vaga:100k","Brachionus_angularis"],
            ["4",7296630,7629422,"5",204286,1699989,"+","Adineta_vaga:100k","Brachionus_angularis"],
            ["5",1957266,2863271,"8",339413,2713622,"+","Adineta_vaga:100k","Brachionus_angularis"],
            ["6",5699009,5762272,"1",5134849,5175869,"+","Adineta_vaga:100k","Brachionus_angularis"],
            ["1",2116609,2144663,"36",198711,291091,"+","Adineta_vaga:100k","Brachionus_koreanus1"],
            ["2",14133754,14799343,"30",115486,458293,"+","Adineta_vaga:100k","Brachionus_koreanus1"],
            ["3",6922680,7543645,"3",2396660,2531648,"+","Adineta_vaga:100k","Brachionus_koreanus1"],
            ["3",8149650,9625642,"23",837241,1157246,"+","Adineta_vaga:100k","Brachionus_koreanus1"],
            ["4",7629060,7662784,"6",1384987,1391245,"+","Adineta_vaga:100k","Brachionus_koreanus1"],
            ["1",4594703,4597541,"19",1308074,1310442,"+","Adineta_vaga:100k","Brachionus_koreanus2"],
            ["1",7693115,7710556,"6",683511,692393,"+","Adineta_vaga:100k","Brachionus_koreanus2"],
            ["1",8601918,8844268,"8",2035490,3693261,"+","Adineta_vaga:100k","Brachionus_koreanus2"],
        ]
        # example sizes: chr1 = 16 Mb
        sizes = {"1":16000000,"2":15000000,"3":13000000,"4":8000000,"5":10000000,"6":11000000}
        build_html(pd.DataFrame(rows, columns=SYN_COLS),
                   ["1","2","3","4","5","6"], None, sizes, "eh_interactive_output")
    else:
        main()
