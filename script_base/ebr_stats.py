#!/usr/bin/env python3
"""
ebr_stats.py — EBR statistics reporter for the SYBR pipeline.

Reads EBRs.txt (chr  start  end  species) and produces:
  • EBRs_stats.txt  — plain-text summary tables
  • EBRs_stats.html — interactive HTML report with charts (Chart.js)

Usage:
    python3 ebr_stats.py <EBRs.txt> [output_dir]

If output_dir is omitted the files are written next to the input file.
"""

import sys
import os
import math
from collections import defaultdict


# ─────────────────────────────────────────────────────────────────────────────
#  Parse
# ─────────────────────────────────────────────────────────────────────────────
def parse_ebrs(path):
    records = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 4:
                print(f"  [WARN] line {lineno}: expected 4 cols, got {len(parts)} — skipped",
                      file=sys.stderr)
                continue
            chr_id, start, end, species = parts
            try:
                start = int(start)
                end   = int(end)
            except ValueError:
                print(f"  [WARN] line {lineno}: non-integer coords — skipped", file=sys.stderr)
                continue
            size = end - start
            if size < 0:
                # swap coords if inverted
                start, end = end, start
                size = end - start
            records.append({
                "chr":     chr_id,
                "start":   start,
                "end":     end,
                "size":    size,
                "species": species.strip(),
            })
    return records


# ─────────────────────────────────────────────────────────────────────────────
#  Compute stats
# ─────────────────────────────────────────────────────────────────────────────
def compute_stats(records):
    total = len(records)
    sizes = [r["size"] for r in records]

    def size_stats(sz_list):
        if not sz_list:
            return dict(n=0, total_bp=0, mean=0, median=0, min=0, max=0, std=0)
        n = len(sz_list)
        s = sorted(sz_list)
        total_bp = sum(s)
        mean = total_bp / n
        median = (s[n // 2] if n % 2 == 1
                  else (s[n // 2 - 1] + s[n // 2]) / 2)
        mn = s[0]
        mx = s[-1]
        variance = sum((x - mean) ** 2 for x in s) / n
        std = math.sqrt(variance)
        return dict(n=n, total_bp=total_bp, mean=mean,
                    median=median, min=mn, max=mx, std=std)

    # Global
    global_stats = size_stats(sizes)

    # Per species
    by_species = defaultdict(list)
    for r in records:
        by_species[r["species"]].append(r["size"])
    species_stats = {sp: size_stats(szs) for sp, szs in sorted(by_species.items())}

    # Per chromosome
    by_chr = defaultdict(list)
    for r in records:
        by_chr[r["chr"]].append(r["size"])

    def chr_sort_key(c):
        try:
            return (0, int(c))
        except ValueError:
            return (1, c)

    chr_stats = {ch: size_stats(szs)
                 for ch, szs in sorted(by_chr.items(), key=lambda x: chr_sort_key(x[0]))}

    # Per chromosome per species
    by_chr_sp = defaultdict(lambda: defaultdict(list))
    for r in records:
        by_chr_sp[r["chr"]][r["species"]].append(r["size"])

    return dict(
        total=total,
        global_stats=global_stats,
        species_stats=species_stats,
        chr_stats=chr_stats,
        by_chr_sp=by_chr_sp,
        by_species=by_species,
        by_chr=by_chr,
    )


# ─────────────────────────────────────────────────────────────────────────────
#  Text report
# ─────────────────────────────────────────────────────────────────────────────
def write_text(stats, out_path, input_file):
    total = stats["total"]
    g     = stats["global_stats"]

    def fmt_bp(n):
        return f"{n:,.0f} bp"

    def pct(n):
        return f"{100 * n / total:.1f}%" if total else "0.0%"

    lines = []
    w = 70
    lines.append("=" * w)
    lines.append("  SYBR Pipeline — EBR Statistics Report")
    lines.append(f"  Input: {os.path.basename(input_file)}")
    lines.append("=" * w)
    lines.append("")

    # ── Global summary ────────────────────────────────────────────────────
    lines.append("─" * w)
    lines.append("  GLOBAL SUMMARY")
    lines.append("─" * w)
    lines.append(f"  Total EBRs               : {total:,}")
    lines.append(f"  Total EBR sequence (bp)  : {g['total_bp']:,}")
    lines.append(f"  Minimum EBR size         : {g['min']:,} bp")
    lines.append(f"  Maximum EBR size         : {g['max']:,} bp")
    lines.append(f"  Mean EBR size            : {g['mean']:,.1f} bp")
    lines.append(f"  Median EBR size          : {g['median']:,.1f} bp")
    lines.append(f"  Std deviation            : {g['std']:,.1f} bp")
    lines.append(f"  Chromosomes represented  : {len(stats['chr_stats'])}")
    lines.append(f"  Species represented      : {len(stats['species_stats'])}")
    lines.append("")

    # ── Per species ───────────────────────────────────────────────────────
    lines.append("─" * w)
    lines.append("  EBR COUNTS AND SIZES BY SPECIES")
    lines.append("─" * w)
    hdr = f"  {'Species':<28} {'Count':>6} {'%Total':>7}  {'Total bp':>12}  {'Min bp':>10}  {'Mean bp':>10}  {'Max bp':>10}"
    lines.append(hdr)
    lines.append("  " + "-" * (w - 2))
    for sp, s in stats["species_stats"].items():
        lines.append(
            f"  {sp:<28} {s['n']:>6} {pct(s['n']):>7}  {s['total_bp']:>12,}  "
            f"{s['min']:>10,}  {s['mean']:>10,.1f}  {s['max']:>10,}"
        )
    lines.append("")

    # ── Per chromosome ────────────────────────────────────────────────────
    lines.append("─" * w)
    lines.append("  EBR COUNTS AND SIZES BY CHROMOSOME")
    lines.append("─" * w)
    hdr = f"  {'Chr':>5} {'Count':>6} {'%Total':>7}  {'Total bp':>12}  {'Min bp':>10}  {'Mean bp':>10}  {'Max bp':>10}"
    lines.append(hdr)
    lines.append("  " + "-" * (w - 2))
    for ch, s in stats["chr_stats"].items():
        lines.append(
            f"  {ch:>5} {s['n']:>6} {pct(s['n']):>7}  {s['total_bp']:>12,}  "
            f"{s['min']:>10,}  {s['mean']:>10,.1f}  {s['max']:>10,}"
        )
    lines.append("")

    # ── Per chromosome × species ──────────────────────────────────────────
    lines.append("─" * w)
    lines.append("  EBR COUNTS PER CHROMOSOME PER SPECIES")
    lines.append("─" * w)
    all_species = sorted(stats["species_stats"].keys())
    sp_width = max(len(s) for s in all_species)
    col_w = max(sp_width, 6)

    header_parts = ["  " + f"{'Chr':>5}"]
    for sp in all_species:
        header_parts.append(f"  {sp:>{col_w}}")
    header_parts.append(f"  {'Total':>6}")
    lines.append("".join(header_parts))
    lines.append("  " + "-" * (w - 2))

    for ch in stats["chr_stats"]:
        row = [f"  {ch:>5}"]
        ch_total = 0
        for sp in all_species:
            n = len(stats["by_chr_sp"][ch].get(sp, []))
            ch_total += n
            row.append(f"  {n:>{col_w}}")
        row.append(f"  {ch_total:>6}")
        lines.append("".join(row))
    lines.append("")
    lines.append("=" * w)

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  Text report   → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
#  HTML report
# ─────────────────────────────────────────────────────────────────────────────
def write_html(stats, records, out_path, input_file):
    import json

    total   = stats["total"]
    g       = stats["global_stats"]
    sp_data = stats["species_stats"]
    ch_data = stats["chr_stats"]
    all_sp  = sorted(sp_data.keys())
    all_ch  = list(ch_data.keys())

    def pct(n):
        return round(100 * n / total, 2) if total else 0

    # Colour palette per species (up to 8 species)
    PALETTE = [
        "#378ADD", "#1D9E75", "#D85A30", "#D4537E",
        "#7F77DD", "#BA7517", "#639922", "#E24B4A",
    ]
    sp_colours = {sp: PALETTE[i % len(PALETTE)] for i, sp in enumerate(all_sp)}

    # ── Data for charts ───────────────────────────────────────────────────
    # 1. Species count bar
    sp_counts   = [sp_data[sp]["n"] for sp in all_sp]
    sp_pcts     = [pct(sp_data[sp]["n"]) for sp in all_sp]
    sp_colours_list = [sp_colours[sp] for sp in all_sp]

    # 2. Chr count bar (stacked by species)
    chr_sp_counts = {sp: [len(stats["by_chr_sp"][ch].get(sp, [])) for ch in all_ch]
                     for sp in all_sp}

    # 3. Mean size per species
    sp_means = [round(sp_data[sp]["mean"]) for sp in all_sp]
    sp_mins  = [sp_data[sp]["min"] for sp in all_sp]
    sp_maxs  = [sp_data[sp]["max"] for sp in all_sp]

    # 4. Size distribution histogram (all EBRs, 20 bins)
    sizes = [r["size"] for r in records]
    if sizes:
        hist_min = min(sizes)
        hist_max = max(sizes)
        n_bins   = 20
        bin_w    = max(1, (hist_max - hist_min) // n_bins)
        hist_bins   = []
        hist_counts = []
        for i in range(n_bins):
            lo = hist_min + i * bin_w
            hi = lo + bin_w
            hist_bins.append(f"{lo//1000}k")
            hist_counts.append(sum(1 for s in sizes if lo <= s < hi))
        # last bin is inclusive
        hist_counts[-1] += sum(1 for s in sizes if s == hist_max)
    else:
        hist_bins = hist_counts = []

    # 5. Chr total bp bar
    chr_total_bp = [ch_data[ch]["total_bp"] for ch in all_ch]

    # ── Tables ────────────────────────────────────────────────────────────
    def stat_row(label, s):
        return f"""
        <tr>
          <td>{label}</td>
          <td>{s['n']:,}</td>
          <td>{pct(s['n']):.1f}%</td>
          <td>{s['total_bp']:,}</td>
          <td>{s['min']:,}</td>
          <td>{s['median']:,.0f}</td>
          <td>{s['mean']:,.1f}</td>
          <td>{s['max']:,}</td>
          <td>{s['std']:,.1f}</td>
        </tr>"""

    sp_rows = "".join(stat_row(sp, sp_data[sp]) for sp in all_sp)

    chr_rows = "".join(stat_row(f"Chr {ch}", ch_data[ch]) for ch in all_ch)

    # Chr × species table
    sp_header = "".join(f"<th>{sp}</th>" for sp in all_sp)
    cross_rows = []
    for ch in all_ch:
        cells = "".join(
            f"<td>{len(stats['by_chr_sp'][ch].get(sp, [])):,}</td>"
            for sp in all_sp
        )
        row_total = ch_data[ch]["n"]
        cross_rows.append(f"<tr><td><b>Chr {ch}</b></td>{cells}<td><b>{row_total:,}</b></td></tr>")
    cross_rows_html = "\n".join(cross_rows)

    # Totals row
    sp_totals = "".join(f"<td><b>{sp_data[sp]['n']:,}</b></td>" for sp in all_sp)

    # ── Inline JS data ────────────────────────────────────────────────────
    chart_data = dict(
        allSp     = all_sp,
        allCh     = [f"Chr {c}" for c in all_ch],
        spCounts  = sp_counts,
        spPcts    = sp_pcts,
        spColours = sp_colours_list,
        chrSpCounts = chr_sp_counts,
        spMeans   = sp_means,
        spMins    = sp_mins,
        spMaxs    = sp_maxs,
        histBins  = hist_bins,
        histCounts= hist_counts,
        chrTotalBp= chr_total_bp,
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>EBR Statistics — {os.path.basename(input_file)}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js"></script>
<style>
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         font-size: 13px; color: #1a1a1a; background: #f4f5f7; line-height: 1.5; }}
  .page {{ max-width: 1200px; margin: 0 auto; padding: 24px 20px 60px; }}
  h1 {{ font-size: 20px; font-weight: 600; color: #111; margin-bottom: 4px; }}
  .subtitle {{ font-size: 12px; color: #666; margin-bottom: 24px; }}
  h2 {{ font-size: 14px; font-weight: 600; color: #222; margin-bottom: 12px; letter-spacing: .02em; }}

  /* cards */
  .card {{ background: #fff; border: 0.5px solid #ddd; border-radius: 10px;
           padding: 20px 24px; margin-bottom: 20px; }}

  /* summary boxes */
  .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(160px,1fr)); gap: 12px; }}
  .metric {{ background: #f8f9fb; border: 0.5px solid #e2e4e8; border-radius: 8px;
             padding: 14px 16px; }}
  .metric .val {{ font-size: 22px; font-weight: 600; color: #185FA5; }}
  .metric .lbl {{ font-size: 11px; color: #666; margin-top: 2px; }}

  /* tables */
  .tbl-wrap {{ overflow-x: auto; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
  thead tr {{ background: #185FA5; color: #fff; }}
  thead th {{ padding: 8px 10px; text-align: right; font-weight: 500; white-space: nowrap; }}
  thead th:first-child {{ text-align: left; }}
  tbody tr:nth-child(even) {{ background: #f7f8fb; }}
  tbody tr:hover {{ background: #eaf2ff; }}
  td, th {{ padding: 7px 10px; border-bottom: 0.5px solid #e8e8e8; }}
  td {{ text-align: right; }}
  td:first-child {{ text-align: left; font-weight: 500; }}
  tfoot td {{ background: #e8f0fb; font-weight: 600; border-top: 1px solid #b0c8f0; }}

  /* charts grid */
  .chart-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 20px; }}
  .chart-grid .card {{ margin-bottom: 0; }}
  @media (max-width: 760px) {{ .chart-grid {{ grid-template-columns: 1fr; }} }}
  .chart-wrap {{ position: relative; height: 260px; }}
  .chart-wrap-tall {{ position: relative; height: 320px; }}
</style>
</head>
<body>
<div class="page">

  <h1>EBR Statistics Report</h1>
  <p class="subtitle">Input: <code>{os.path.basename(input_file)}</code></p>

  <!-- ── Global summary metrics ── -->
  <div class="card">
    <h2>Global summary</h2>
    <div class="metric-grid">
      <div class="metric"><div class="val">{total:,}</div><div class="lbl">Total EBRs</div></div>
      <div class="metric"><div class="val">{len(all_sp)}</div><div class="lbl">Species</div></div>
      <div class="metric"><div class="val">{len(all_ch)}</div><div class="lbl">Chromosomes</div></div>
      <div class="metric"><div class="val">{g['min']:,}</div><div class="lbl">Min EBR size (bp)</div></div>
      <div class="metric"><div class="val">{g['median']:,.0f}</div><div class="lbl">Median EBR size (bp)</div></div>
      <div class="metric"><div class="val">{g['mean']:,.0f}</div><div class="lbl">Mean EBR size (bp)</div></div>
      <div class="metric"><div class="val">{g['max']:,}</div><div class="lbl">Max EBR size (bp)</div></div>
      <div class="metric"><div class="val">{g['std']:,.0f}</div><div class="lbl">Std deviation (bp)</div></div>
      <div class="metric"><div class="val">{g['total_bp']:,}</div><div class="lbl">Total EBR sequence (bp)</div></div>
    </div>
  </div>

  <!-- ── Charts row 1 ── -->
  <div class="chart-grid">
    <div class="card">
      <h2>EBR count per species</h2>
      <div class="chart-wrap"><canvas id="spCountChart"></canvas></div>
    </div>
    <div class="card">
      <h2>EBR size distribution (all EBRs)</h2>
      <div class="chart-wrap"><canvas id="histChart"></canvas></div>
    </div>
  </div>

  <!-- ── Charts row 2 ── -->
  <div class="chart-grid">
    <div class="card">
      <h2>EBRs per chromosome (stacked by species)</h2>
      <div class="chart-wrap-tall"><canvas id="chrStackChart"></canvas></div>
    </div>
    <div class="card">
      <h2>Min / Mean / Max EBR size per species (bp)</h2>
      <div class="chart-wrap-tall"><canvas id="spSizeChart"></canvas></div>
    </div>
  </div>

  <!-- ── Chart row 3: chr total bp ── -->
  <div class="card">
    <h2>Total EBR sequence per chromosome (bp)</h2>
    <div class="chart-wrap"><canvas id="chrBpChart"></canvas></div>
  </div>

  <!-- ── Table: per species ── -->
  <div class="card">
    <h2>Statistics by species</h2>
    <div class="tbl-wrap">
      <table>
        <thead><tr>
          <th>Species</th><th>Count</th><th>% Total</th>
          <th>Total bp</th><th>Min bp</th><th>Median bp</th>
          <th>Mean bp</th><th>Max bp</th><th>Std dev bp</th>
        </tr></thead>
        <tbody>{sp_rows}</tbody>
        <tfoot><tr>
          <td>All species</td>
          <td>{total:,}</td><td>100.0%</td>
          <td>{g['total_bp']:,}</td>
          <td>{g['min']:,}</td><td>{g['median']:,.0f}</td>
          <td>{g['mean']:,.1f}</td><td>{g['max']:,}</td>
          <td>{g['std']:,.1f}</td>
        </tr></tfoot>
      </table>
    </div>
  </div>

  <!-- ── Table: per chromosome ── -->
  <div class="card">
    <h2>Statistics by chromosome</h2>
    <div class="tbl-wrap">
      <table>
        <thead><tr>
          <th>Chromosome</th><th>Count</th><th>% Total</th>
          <th>Total bp</th><th>Min bp</th><th>Median bp</th>
          <th>Mean bp</th><th>Max bp</th><th>Std dev bp</th>
        </tr></thead>
        <tbody>{chr_rows}</tbody>
      </table>
    </div>
  </div>

  <!-- ── Table: cross chromosome × species ── -->
  <div class="card">
    <h2>EBR counts per chromosome × species</h2>
    <div class="tbl-wrap">
      <table>
        <thead><tr>
          <th>Chromosome</th>{sp_header}<th>Row total</th>
        </tr></thead>
        <tbody>{cross_rows_html}</tbody>
        <tfoot><tr>
          <td>Column total</td>{sp_totals}<td><b>{total:,}</b></td>
        </tr></tfoot>
      </table>
    </div>
  </div>

</div><!-- /page -->

<script>
const D = {json.dumps(chart_data)};

const sp   = D.allSp;
const ch   = D.allCh;
const pal  = D.spColours;

function mkAlpha(hex, a) {{
  const r=parseInt(hex.slice(1,3),16), g=parseInt(hex.slice(3,5),16), b=parseInt(hex.slice(5,7),16);
  return `rgba(${{r}},${{g}},${{b}},${{a}})`;
}}

// 1. Species count bar
new Chart(document.getElementById('spCountChart'), {{
  type: 'bar',
  data: {{
    labels: sp,
    datasets: [{{
      label: 'EBR count',
      data: D.spCounts,
      backgroundColor: pal.map(c => mkAlpha(c, 0.8)),
      borderColor: pal,
      borderWidth: 1,
      borderRadius: 4,
    }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ display: false }},
               tooltip: {{ callbacks: {{ afterLabel: (ctx) => `${{D.spPcts[ctx.dataIndex]}}% of total` }} }} }},
    scales: {{
      x: {{ ticks: {{ font: {{ size: 11 }}, maxRotation: 30 }} }},
      y: {{ beginAtZero: true, ticks: {{ font: {{ size: 11 }} }}, title: {{ display: true, text: 'EBR count', font: {{ size: 11 }} }} }}
    }}
  }}
}});

// 2. Histogram
new Chart(document.getElementById('histChart'), {{
  type: 'bar',
  data: {{
    labels: D.histBins,
    datasets: [{{
      label: 'EBRs',
      data: D.histCounts,
      backgroundColor: mkAlpha('#378ADD', 0.75),
      borderColor: '#185FA5',
      borderWidth: 1,
      borderRadius: 3,
    }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ display: false }} }},
    scales: {{
      x: {{ ticks: {{ font: {{ size: 10 }}, maxRotation: 45 }},
           title: {{ display: true, text: 'EBR size', font: {{ size: 11 }} }} }},
      y: {{ beginAtZero: true, ticks: {{ font: {{ size: 11 }} }},
           title: {{ display: true, text: 'Count', font: {{ size: 11 }} }} }}
    }}
  }}
}});

// 3. Stacked bar chr × species
new Chart(document.getElementById('chrStackChart'), {{
  type: 'bar',
  data: {{
    labels: ch,
    datasets: sp.map((s, i) => ({{
      label: s,
      data: D.chrSpCounts[s],
      backgroundColor: mkAlpha(pal[i], 0.8),
      borderColor: pal[i],
      borderWidth: 0.5,
      borderRadius: 2,
    }}))
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ position: 'bottom', labels: {{ font: {{ size: 11 }}, boxWidth: 12 }} }} }},
    scales: {{
      x: {{ stacked: true, ticks: {{ font: {{ size: 11 }} }} }},
      y: {{ stacked: true, beginAtZero: true,
           ticks: {{ font: {{ size: 11 }} }},
           title: {{ display: true, text: 'EBR count', font: {{ size: 11 }} }} }}
    }}
  }}
}});

// 4. Min/Mean/Max per species (floating bar trick via min + [min→mean] + [mean→max])
new Chart(document.getElementById('spSizeChart'), {{
  type: 'bar',
  data: {{
    labels: sp,
    datasets: [
      {{
        label: 'Min',
        data: D.spMins,
        backgroundColor: sp.map((_, i) => mkAlpha(pal[i], 0.35)),
        borderColor: 'transparent', borderWidth: 0, borderRadius: 3,
      }},
      {{
        label: 'Mean',
        data: D.spMeans,
        backgroundColor: sp.map((_, i) => mkAlpha(pal[i], 0.75)),
        borderColor: sp.map((_, i) => pal[i]), borderWidth: 1, borderRadius: 3,
      }},
      {{
        label: 'Max',
        data: D.spMaxs,
        backgroundColor: 'transparent',
        borderColor: sp.map((_, i) => pal[i]),
        borderWidth: 2, type: 'bar', borderRadius: 3,
      }}
    ]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ position: 'bottom', labels: {{ font: {{ size: 11 }}, boxWidth: 12 }} }} }},
    scales: {{
      x: {{ ticks: {{ font: {{ size: 11 }}, maxRotation: 30 }} }},
      y: {{ beginAtZero: true, ticks: {{ font: {{ size: 11 }} }},
           title: {{ display: true, text: 'Size (bp)', font: {{ size: 11 }} }} }}
    }}
  }}
}});

// 5. Chr total bp
new Chart(document.getElementById('chrBpChart'), {{
  type: 'bar',
  data: {{
    labels: ch,
    datasets: [{{
      label: 'Total EBR bp',
      data: D.chrTotalBp,
      backgroundColor: ch.map((_, i) => mkAlpha('#1D9E75', 0.65 + (i % 3) * 0.1)),
      borderColor: '#0F6E56',
      borderWidth: 1,
      borderRadius: 4,
    }}]
  }},
  options: {{
    responsive: true, maintainAspectRatio: false,
    plugins: {{ legend: {{ display: false }} }},
    scales: {{
      x: {{ ticks: {{ font: {{ size: 11 }} }} }},
      y: {{ beginAtZero: true, ticks: {{ font: {{ size: 11 }} }},
           title: {{ display: true, text: 'Total bp', font: {{ size: 11 }} }} }}
    }}
  }}
}});
</script>
</body>
</html>"""

    with open(out_path, "w") as fh:
        fh.write(html)
    print(f"  HTML report   → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    input_path = sys.argv[1]
    if not os.path.isfile(input_path):
        print(f"[ERROR] File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(input_path))
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n  Reading: {input_path}")
    records = parse_ebrs(input_path)
    if not records:
        print("[ERROR] No valid records found.", file=sys.stderr)
        sys.exit(1)
    print(f"  Parsed {len(records)} EBRs")

    stats = compute_stats(records)

    base = os.path.splitext(os.path.basename(input_path))[0]
    write_text(stats, os.path.join(out_dir, f"{base}_stats.txt"),  input_path)
    write_html(stats, records, os.path.join(out_dir, f"{base}_stats.html"), input_path)

    print(f"\n  Done.\n")


if __name__ == "__main__":
    main()
