#!/usr/bin/env python3
"""PCoA ordination of L4.11 strains — fast version using ete3 + correction."""

import csv, sys
import numpy as np
from ete3 import Tree
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
import plotly.graph_objects as go
import plotly.io as pio
from collections import Counter

# ── Paths ────────────────────────────────────────────────────────────────────
TREE = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/T3.raxml.bestTree"
CSV  = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/strains.csv"
OUT  = "/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/pcoa_l411.html"

# ── L4.11.1 strain list ──────────────────────────────────────────────────────
L4111 = set("""ERR13289532 ERR2513221 ERR2514395 ERR4553410 ERR4553419 ERR4553448
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

# ── Load tree ────────────────────────────────────────────────────────────────
print("Loading tree with ete3...")
t = Tree(TREE)
leaves = t.get_leaves()
tip_names = [l.name for l in leaves]
n = len(tip_names)
print(f"  {n} tips")

import re

def extract_accession(leaf_name):
    """Extract SRR/ERR accession from prefixed leaf name."""
    m = re.search(r'(ERR\d+|SRR\d+)', leaf_name)
    return m.group(1) if m else leaf_name

# ── Cophenetic distance matrix (BEFORE metadata, to get ordered_names) ───────
print("Computing cophenetic distance matrix...")
coph_data, ordered_names = t.cophenetic_matrix()
dist_np = np.array(coph_data, dtype=float)

# Make sure it's symmetric and zero-diagonal
dist_np = (dist_np + dist_np.T) / 2
np.fill_diagonal(dist_np, 0)

print(f"  Matrix shape: {dist_np.shape}")
print(f"  Range: [{dist_np[dist_np > 0].min():.6f}, {dist_np.max():.6f}]")

# ── Load metadata ────────────────────────────────────────────────────────────
meta = {}
with open(CSV) as f:
    for row in csv.DictReader(f):
        country = row.get("Country", "").strip()
        meta[row["Id"]] = {"country": country if country else "Unknown"}

# Build accession->sublineage mapping
acc_to_sub = {acc: "L4.11.1" for acc in L4111}

# Map ordered leaf names to metadata using extracted accessions
leaf_meta = {}
n_l4111 = 0
for tip in ordered_names:
    acc = extract_accession(tip)
    sub = acc_to_sub.get(acc, "L4.11.2")
    country = meta.get(acc, {}).get("country", "Unknown")
    leaf_meta[tip] = {"sublineage": sub, "country": country}
    if sub == "L4.11.1":
        n_l4111 += 1

print(f"  L4.11.1: {n_l4111}, L4.11.2: {n - n_l4111}")


# ── PCoA with manual eigendecomposition ──────────────────────────────────────
print("Running PCoA...")
dm = DistanceMatrix(dist_np, ids=ordered_names)

# Use a classical PCoA with manual eigenvalue handling
from scipy.spatial.distance import squareform
from scipy.linalg import eigh

# Classical MDS / PCoA
D2 = dist_np ** 2
H = np.eye(n) - np.ones((n, n)) / n
B = -0.5 * H @ D2 @ H

# Eigendecomposition
eigenvalues, eigenvectors = eigh(B)
# Sort descending
idx = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# Keep only positive eigenvalues
pos_mask = eigenvalues > 0
n_pos = pos_mask.sum()
print(f"  {n_pos} positive eigenvalues out of {n}")

# Take top 3
top_k = 3
evals = eigenvalues[:top_k]
evecs = eigenvectors[:, :top_k]

# Coordinates
coords = evecs * np.sqrt(np.abs(evals))

# Proportion explained (only positive eigenvalues)
sum_pos = eigenvalues[eigenvalues > 0].sum()
prop_explained = evals / sum_pos * 100

pc1, pc2, pc3 = coords[:, 0], coords[:, 1], coords[:, 2]
pct1, pct2, pct3 = prop_explained[0], prop_explained[1], prop_explained[2]
print(f"  PC1: {pct1:.1f}%, PC2: {pct2:.1f}%, PC3: {pct3:.1f}%")

# ── Metadata arrays ──────────────────────────────────────────────────────────
sublineages = [leaf_meta.get(nm, {}).get("sublineage", "L4.11.2") for nm in ordered_names]
countries   = [leaf_meta.get(nm, {}).get("country", "Unknown") for nm in ordered_names]

# ── 2D by sub-lineage ────────────────────────────────────────────────────────
fig2d = go.Figure()
for sub, col, sym in [("L4.11.1", "#2563eb", "circle"), ("L4.11.2", "#dc2626", "diamond")]:
    mask = [i for i, s in enumerate(sublineages) if s == sub]
    fig2d.add_trace(go.Scatter(
        x=[pc1[i] for i in mask], y=[pc2[i] for i in mask],
        mode="markers", name=f"{sub} (n={len(mask)})",
        marker=dict(size=7, color=col, opacity=0.7, symbol=sym,
                    line=dict(width=0.5, color="white")),
        text=[f"{ordered_names[i]}<br>{countries[i]}" for i in mask],
        hovertemplate="<b>%{text}</b><br>PC1=%{x:.4f}<br>PC2=%{y:.4f}<extra></extra>"
    ))
fig2d.update_layout(
    title=dict(text="PCoA of L4.11 patristic distances — by sub-lineage", font=dict(size=16)),
    xaxis_title=f"PC1 ({pct1:.1f}%)", yaxis_title=f"PC2 ({pct2:.1f}%)",
    template="plotly_white", width=950, height=680,
    legend=dict(x=0.02, y=0.98, bgcolor="rgba(255,255,255,0.85)",
                font=dict(size=13), bordercolor="#ccc", borderwidth=1)
)

# ── 2D by country ────────────────────────────────────────────────────────────
fig_country = go.Figure()
country_counts = Counter(countries)
top_countries = [c for c, _ in country_counts.most_common(10) if c != "Unknown"]
country_colors = {
    "Peru": "#dc2626", "Bangladesh": "#2563eb", "South Africa": "#16a34a",
    "USA": "#f59e0b", "Uganda": "#8b5cf6", "Argentina": "#78716c",
    "China": "#ec4899", "Brazil": "#0891b2", "India": "#65a30d",
    "Vietnam": "#0d9488",
}
for country in top_countries:
    mask = [i for i, c in enumerate(countries) if c == country]
    col = country_colors.get(country, "#6b7280")
    fig_country.add_trace(go.Scatter(
        x=[pc1[i] for i in mask], y=[pc2[i] for i in mask],
        mode="markers", name=f"{country} ({len(mask)})",
        marker=dict(size=7, color=col, opacity=0.75,
                    line=dict(width=0.5, color="white")),
        text=[f"{ordered_names[i]}<br>{sublineages[i]}" for i in mask],
        hovertemplate="<b>%{text}</b><br>PC1=%{x:.4f}<br>PC2=%{y:.4f}<extra></extra>"
    ))
mask_other = [i for i, c in enumerate(countries) if c not in top_countries]
if mask_other:
    fig_country.add_trace(go.Scatter(
        x=[pc1[i] for i in mask_other], y=[pc2[i] for i in mask_other],
        mode="markers", name=f"Other/Unknown ({len(mask_other)})",
        marker=dict(size=4, color="#d1d5db", opacity=0.5,
                    line=dict(width=0.3, color="white")),
        text=[f"{ordered_names[i]}<br>{countries[i]}<br>{sublineages[i]}" for i in mask_other],
        hovertemplate="<b>%{text}</b><br>PC1=%{x:.4f}<br>PC2=%{y:.4f}<extra></extra>"
    ))
fig_country.update_layout(
    title=dict(text="PCoA of L4.11 — colored by country of sampling", font=dict(size=16)),
    xaxis_title=f"PC1 ({pct1:.1f}%)", yaxis_title=f"PC2 ({pct2:.1f}%)",
    template="plotly_white", width=950, height=680,
    legend=dict(x=0.02, y=0.98, bgcolor="rgba(255,255,255,0.85)",
                font=dict(size=12), bordercolor="#ccc", borderwidth=1)
)

# ── 3D by sub-lineage ────────────────────────────────────────────────────────
fig3d = go.Figure()
for sub, col in [("L4.11.1", "#2563eb"), ("L4.11.2", "#dc2626")]:
    mask = [i for i, s in enumerate(sublineages) if s == sub]
    fig3d.add_trace(go.Scatter3d(
        x=[pc1[i] for i in mask], y=[pc2[i] for i in mask], z=[pc3[i] for i in mask],
        mode="markers", name=f"{sub} (n={len(mask)})",
        marker=dict(size=3, color=col, opacity=0.7,
                    line=dict(width=0.3, color="white")),
        text=[f"{ordered_names[i]}<br>{countries[i]}" for i in mask],
        hovertemplate="<b>%{text}</b><br>PC1=%{x:.4f}<br>PC2=%{y:.4f}<br>PC3=%{z:.4f}<extra></extra>"
    ))
fig3d.update_layout(
    title=dict(text="PCoA 3D — L4.11 genomic distances", font=dict(size=16)),
    scene=dict(
        xaxis_title=f"PC1 ({pct1:.1f}%)",
        yaxis_title=f"PC2 ({pct2:.1f}%)",
        zaxis_title=f"PC3 ({pct3:.1f}%)",
    ),
    template="plotly_white", width=950, height=720,
    legend=dict(x=0.02, y=0.98)
)

# ── Save combined HTML ───────────────────────────────────────────────────────
with open(OUT, "w") as f:
    f.write("<html><head><title>PCoA L4.11</title></head><body>\n")
    f.write("<h1>PCoA Ordination — L4.11 Patristic Distances (1291 strains)</h1>\n")
    f.write(f"<p>Variance explained: PC1={pct1:.1f}%, PC2={pct2:.1f}%, PC3={pct3:.1f}%</p>\n")
    f.write("<h2>1. By Sub-lineage (2D)</h2>\n")
    f.write(pio.to_html(fig2d, full_html=False, include_plotlyjs="cdn"))
    f.write("<h2>2. By Country (2D)</h2>\n")
    f.write(pio.to_html(fig_country, full_html=False, include_plotlyjs=False))
    f.write("<h2>3. By Sub-lineage (3D)</h2>\n")
    f.write(pio.to_html(fig3d, full_html=False, include_plotlyjs=False))
    f.write("</body></html>\n")

print(f"\nDone! Output: {OUT}")
