#!/usr/bin/env python3
"""PCoA ordination of L4.11 strains — 3 sub-groups: L4.11.2, core L4.11.1, Proto-L4.11.1."""

import csv, re
import numpy as np
from ete3 import Tree
import plotly.graph_objects as go
import plotly.io as pio
from collections import Counter
from scipy.linalg import eigh

# ── Paths ────────────────────────────────────────────────────────────────────
TREE = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/T3.raxml.bestTree"
CSV  = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/strains.csv"
OUT_HTML = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/pcoa_l411_intra.html"
OUT_PNG  = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/article/figures/pcoa_l411.png"

# ── L4.11.1 strain lists ─────────────────────────────────────────────────────
# Proto-L4.11.1 (basal sub-clade: 5 South Africa + 2 UK)
PROTO = set("""SRR3675589 ERR13289532 SRR35281598 SRR35281599
SRR35281596 SRR35281602 SRR35281600""".split())

# All L4.11.1 (91 strains total)
ALL_L4111 = set("""ERR13289532 ERR2513221 ERR2514395 ERR4553410 ERR4553419 ERR4553448
ERR4553470 ERR4553478 ERR4553511 ERR4553515 ERR4553546 ERR4553566 ERR4553568
ERR4553613 ERR4553633 ERR4553665 ERR4553670 ERR4553716 ERR4553720 ERR4553721
ERR4553747 ERR4553770 ERR4553815 ERR4553821 ERR4553824 ERR4553830 ERR4553834
ERR4553841 ERR4553856 ERR4553857 ERR4553887 ERR4553896 ERR4553923 ERR4553942
ERR4553947 ERR4553968 ERR4553972 SRR1049729 SRR1049730 SRR1062930 SRR1140739
SRR21661641 SRR29016766 SRR29016810 SRR29016812 SRR29016829 SRR29016865
SRR29016881 SRR29017028 SRR29017092 SRR29017094 SRR29017107 SRR29017138
SRR29017144 SRR29017147 SRR29017178 SRR29017188 SRR29017213 SRR29017234
SRR29017242 SRR29017384 SRR29017410 SRR29017445 SRR29017548 SRR29017583
SRR29017611 SRR29017620 SRR29017622 SRR29017648 SRR29017661 SRR29017674
SRR29017683 SRR29017690 SRR29017695 SRR29017701 SRR29055490 SRR29341528
SRR29440700 SRR30443755 SRR35281596 SRR35281598 SRR35281599 SRR35281600
SRR35281602 SRR3675589 SRR4423155 SRR4423179 SRR4423181 SRR6650186
SRR6797722 SRR6797801""".split())

# Core L4.11.1 = all L4.11.1 minus Proto
CORE_L4111 = ALL_L4111 - PROTO

def extract_accession(leaf_name):
    m = re.search(r'(ERR\d+|SRR\d+)', leaf_name)
    return m.group(1) if m else leaf_name

def assign_subgroup(acc):
    if acc in PROTO:
        return "Proto-L4.11.1"
    elif acc in CORE_L4111:
        return "L4.11.1"
    else:
        return "L4.11.2"

# ── Load tree and prune to L4.11 only ────────────────────────────────────────
print("Loading tree...")
t = Tree(TREE)
l411_leaves = [l for l in t.get_leaves() if l.name.startswith("L4.11_")]
print(f"  Total: {len(t.get_leaves())} tips, L4.11: {len(l411_leaves)}")

t.prune(l411_leaves, preserve_branch_length=True)
print(f"  Pruned: {len(t.get_leaves())} tips")

# ── Cophenetic distance matrix ───────────────────────────────────────────────
print("Computing cophenetic distance matrix...")
coph_data, ordered_names = t.cophenetic_matrix()
dist_np = np.array(coph_data, dtype=float)
n = len(ordered_names)
dist_np = (dist_np + dist_np.T) / 2
np.fill_diagonal(dist_np, 0)
print(f"  Matrix: {n}x{n}")

# ── Load metadata ────────────────────────────────────────────────────────────
meta = {}
with open(CSV) as f:
    for row in csv.DictReader(f):
        country = row.get("Country", "").strip()
        meta[row["Id"]] = {"country": country if country else "Unknown"}

leaf_meta = {}
counts = Counter()
for tip in ordered_names:
    acc = extract_accession(tip)
    sub = assign_subgroup(acc)
    country = meta.get(acc, {}).get("country", "Unknown")
    leaf_meta[tip] = {"subgroup": sub, "country": country, "accession": acc}
    counts[sub] += 1

for k, v in sorted(counts.items()):
    print(f"  {k}: {v}")

# ── PCoA ─────────────────────────────────────────────────────────────────────
print("Running PCoA...")
D2 = dist_np ** 2
H = np.eye(n) - np.ones((n, n)) / n
B = -0.5 * H @ D2 @ H
eigenvalues, eigenvectors = eigh(B)
idx = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

n_pos = (eigenvalues > 0).sum()
print(f"  {n_pos} positive eigenvalues out of {n}")

top_k = 3
evals = eigenvalues[:top_k]
evecs = eigenvectors[:, :top_k]
coords = evecs * np.sqrt(np.abs(evals))
sum_pos = eigenvalues[eigenvalues > 0].sum()
prop_explained = evals / sum_pos * 100

pc1, pc2, pc3 = coords[:, 0], coords[:, 1], coords[:, 2]
pct1, pct2, pct3 = prop_explained[0], prop_explained[1], prop_explained[2]
print(f"  PC1: {pct1:.1f}%, PC2: {pct2:.1f}%, PC3: {pct3:.1f}%")

# ── Metadata arrays ──────────────────────────────────────────────────────────
subgroups = [leaf_meta[nm]["subgroup"] for nm in ordered_names]
countries = [leaf_meta[nm]["country"] for nm in ordered_names]
accessions = [leaf_meta[nm]["accession"] for nm in ordered_names]

# ── Colors and symbols ───────────────────────────────────────────────────────
GROUP_STYLE = {
    "L4.11.2":       {"color": "#dc2626", "symbol": "diamond",  "size": 6},
    "L4.11.1":       {"color": "#2563eb", "symbol": "circle",   "size": 7},
    "Proto-L4.11.1": {"color": "#16a34a", "symbol": "triangle-up", "size": 9},
}

# ── Figure 1: 2D by sub-group ───────────────────────────────────────────────
fig2d = go.Figure()
for grp in ["L4.11.2", "L4.11.1", "Proto-L4.11.1"]:
    sty = GROUP_STYLE[grp]
    mask = [i for i, s in enumerate(subgroups) if s == grp]
    fig2d.add_trace(go.Scatter(
        x=[pc1[i] for i in mask], y=[pc2[i] for i in mask],
        mode="markers", name=f"{grp} (n={len(mask)})",
        marker=dict(size=sty["size"], color=sty["color"], opacity=0.8,
                    symbol=sty["symbol"], line=dict(width=0.5, color="white")),
        text=[f"{accessions[i]}<br>{countries[i]}" for i in mask],
        hovertemplate="<b>%{text}</b><br>PC1=%{x:.5f}<br>PC2=%{y:.5f}<extra></extra>"
    ))
fig2d.update_layout(
    title=dict(text="PCoA of L4.11 intra-lineage patristic distances", font=dict(size=15)),
    xaxis_title=f"PC1 ({pct1:.1f}%)", yaxis_title=f"PC2 ({pct2:.1f}%)",
    template="plotly_white", width=950, height=680,
    legend=dict(x=0.02, y=0.98, bgcolor="rgba(255,255,255,0.85)",
                font=dict(size=13), bordercolor="#ccc", borderwidth=1)
)

# ── Figure 2: 2D by country ─────────────────────────────────────────────────
fig_country = go.Figure()
country_counts = Counter(countries)
top_countries = [c for c, _ in country_counts.most_common(10) if c != "Unknown"]
country_colors = {
    "Peru": "#dc2626", "Bangladesh": "#2563eb", "South Africa": "#16a34a",
    "United States of America": "#f59e0b", "USA": "#f59e0b",
    "Uganda": "#8b5cf6", "Argentina": "#78716c",
    "China": "#ec4899", "Brazil": "#0891b2", "India": "#65a30d",
    "Vietnam": "#0d9488", "Italy": "#d97706", "Sweden": "#6366f1",
    "United Kingdom of Great Britain and Northern Ireland": "#0ea5e9",
}
for country in top_countries:
    mask = [i for i, c in enumerate(countries) if c == country]
    col = country_colors.get(country, "#6b7280")
    fig_country.add_trace(go.Scatter(
        x=[pc1[i] for i in mask], y=[pc2[i] for i in mask],
        mode="markers", name=f"{country} ({len(mask)})",
        marker=dict(size=7, color=col, opacity=0.75,
                    line=dict(width=0.5, color="white")),
        text=[f"{accessions[i]}<br>{subgroups[i]}" for i in mask],
        hovertemplate="<b>%{text}</b><br>PC1=%{x:.5f}<br>PC2=%{y:.5f}<extra></extra>"
    ))
mask_other = [i for i, c in enumerate(countries) if c not in top_countries]
if mask_other:
    fig_country.add_trace(go.Scatter(
        x=[pc1[i] for i in mask_other], y=[pc2[i] for i in mask_other],
        mode="markers", name=f"Other/Unknown ({len(mask_other)})",
        marker=dict(size=4, color="#d1d5db", opacity=0.5,
                    line=dict(width=0.3, color="white")),
        text=[f"{accessions[i]}<br>{countries[i]}<br>{subgroups[i]}" for i in mask_other],
        hovertemplate="<b>%{text}</b><br>PC1=%{x:.5f}<br>PC2=%{y:.5f}<extra></extra>"
    ))
fig_country.update_layout(
    title=dict(text="PCoA of L4.11 — colored by country", font=dict(size=15)),
    xaxis_title=f"PC1 ({pct1:.1f}%)", yaxis_title=f"PC2 ({pct2:.1f}%)",
    template="plotly_white", width=950, height=680,
    legend=dict(x=0.02, y=0.98, bgcolor="rgba(255,255,255,0.85)",
                font=dict(size=12), bordercolor="#ccc", borderwidth=1)
)

# ── Figure 3: 3D by sub-group ───────────────────────────────────────────────
fig3d = go.Figure()
for grp in ["L4.11.2", "L4.11.1", "Proto-L4.11.1"]:
    sty = GROUP_STYLE[grp]
    mask = [i for i, s in enumerate(subgroups) if s == grp]
    fig3d.add_trace(go.Scatter3d(
        x=[pc1[i] for i in mask], y=[pc2[i] for i in mask], z=[pc3[i] for i in mask],
        mode="markers", name=f"{grp} (n={len(mask)})",
        marker=dict(size=3, color=sty["color"], opacity=0.7,
                    line=dict(width=0.3, color="white")),
        text=[f"{accessions[i]}<br>{countries[i]}" for i in mask],
        hovertemplate="<b>%{text}</b><br>PC1=%{x:.5f}<br>PC2=%{y:.5f}<br>PC3=%{z:.5f}<extra></extra>"
    ))
fig3d.update_layout(
    title=dict(text="PCoA 3D — L4.11 intra-lineage distances", font=dict(size=15)),
    scene=dict(xaxis_title=f"PC1 ({pct1:.1f}%)", yaxis_title=f"PC2 ({pct2:.1f}%)",
               zaxis_title=f"PC3 ({pct3:.1f}%)"),
    template="plotly_white", width=950, height=720, legend=dict(x=0.02, y=0.98)
)

# ── Save HTML ────────────────────────────────────────────────────────────────
with open(OUT_HTML, "w") as f:
    f.write("<html><head><title>PCoA L4.11 (intra-lineage)</title></head><body>\n")
    f.write(f"<h1>PCoA — L4.11 Intra-lineage Distances ({n} strains)</h1>\n")
    f.write(f"<p>Variance explained: PC1={pct1:.1f}%, PC2={pct2:.1f}%, PC3={pct3:.1f}%</p>\n")
    f.write("<h2>1. By Sub-group (2D)</h2>\n")
    f.write(pio.to_html(fig2d, full_html=False, include_plotlyjs="cdn"))
    f.write("<h2>2. By Country (2D)</h2>\n")
    f.write(pio.to_html(fig_country, full_html=False, include_plotlyjs=False))
    f.write("<h2>3. By Sub-group (3D)</h2>\n")
    f.write(pio.to_html(fig3d, full_html=False, include_plotlyjs=False))
    f.write("</body></html>\n")
print(f"\nHTML: {OUT_HTML}")

# ── Publication-quality static PNG ───────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=200)

# Panel A: by sub-group (draw L4.11.2 first, then core, then proto on top)
for grp, col, mk, ms, zord, lab in [
    ("L4.11.2",       "#dc2626", "D", 10, 2, "L4.11.2"),
    ("L4.11.1",       "#2563eb", "o", 14, 3, "L4.11.1 (core)"),
    ("Proto-L4.11.1", "#16a34a", "^", 40, 4, "Proto-L4.11.1"),
]:
    mask = [i for i, s in enumerate(subgroups) if s == grp]
    ax1.scatter([pc1[i] for i in mask], [pc2[i] for i in mask],
                c=col, marker=mk, s=ms, alpha=0.8, edgecolors="white",
                linewidths=0.3, label=f"{lab} (n={len(mask)})", zorder=zord)

ax1.set_xlabel(f"PC1 ({pct1:.1f}%)", fontsize=11)
ax1.set_ylabel(f"PC2 ({pct2:.1f}%)", fontsize=11)
ax1.set_title("(a) By sub-lineage", fontsize=12, fontweight="bold")
ax1.legend(fontsize=8, loc="upper left", framealpha=0.9)
ax1.grid(True, alpha=0.15)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Panel B: by country
for country in top_countries:
    mask = [i for i, c in enumerate(countries) if c == country]
    col = country_colors.get(country, "#6b7280")
    lab = country.replace("United States of America", "USA").replace(
        "United Kingdom of Great Britain and Northern Ireland", "UK")
    ax2.scatter([pc1[i] for i in mask], [pc2[i] for i in mask],
                c=col, s=12, alpha=0.7, edgecolors="white", linewidths=0.3,
                label=f"{lab} ({len(mask)})")

if mask_other:
    ax2.scatter([pc1[i] for i in mask_other], [pc2[i] for i in mask_other],
                c="#d1d5db", s=6, alpha=0.4, edgecolors="white", linewidths=0.2,
                label=f"Other ({len(mask_other)})")

ax2.set_xlabel(f"PC1 ({pct1:.1f}%)", fontsize=11)
ax2.set_ylabel(f"PC2 ({pct2:.1f}%)", fontsize=11)
ax2.set_title("(b) By country of origin", fontsize=12, fontweight="bold")
ax2.legend(fontsize=7, loc="upper left", framealpha=0.9, ncol=2)
ax2.grid(True, alpha=0.15)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(OUT_PNG, bbox_inches="tight", facecolor="white")
print(f"PNG: {OUT_PNG}")

# Also copy script
import shutil
shutil.copy(__file__, "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/pcoa_l411_intra.py")

print("\nDone!")
