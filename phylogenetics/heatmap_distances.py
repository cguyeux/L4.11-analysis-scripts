#!/usr/bin/env python3
"""
Pairwise phylogenetic distance heatmap for L4.11 strains.
Uses the RAxML best tree to compute patristic distances,
then displays a seaborn clustermap colored by sub-lineage.

Strategy: with 828+ tips, a full heatmap is unreadable.
We subsample ~150 strains (proportional to sub-lineage sizes)
to keep the figure legible, while showing the full structure.
"""
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from Bio import Phylo

# ── Configuration ──────────────────────────────────────────────────────────
TREE_FILE = '/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/résultats/T3.raxml.bestTree'
OUTPUT_PNG = '/home/christophe/Documents/codes/MTBC/TBannotator/articles/L4.11/article/figures/heatmap_snp_distances.png'

MAX_STRAINS = 150  # subsample for readability

# L4.11.1 strain IDs
L4111_IDS = set("""
ERR13289532 ERR2513221 ERR2514395 ERR4553410 ERR4553419 ERR4553448 ERR4553470
ERR4553478 ERR4553511 ERR4553515 ERR4553546 ERR4553566 ERR4553568 ERR4553613
ERR4553633 ERR4553665 ERR4553670 ERR4553716 ERR4553720 ERR4553721 ERR4553747
ERR4553770 ERR4553815 ERR4553821 ERR4553824 ERR4553830 ERR4553834 ERR4553841
ERR4553856 ERR4553857 ERR4553887 ERR4553896 ERR4553923 ERR4553942 ERR4553947
ERR4553968 ERR4553972 SRR1049729 SRR1049730 SRR1062930 SRR1140739 SRR21661641
SRR29016766 SRR29016810 SRR29016812 SRR29016829 SRR29016865 SRR29016881
SRR29017028 SRR29017092 SRR29017094 SRR29017107 SRR29017138 SRR29017144
SRR29017147 SRR29017178 SRR29017188 SRR29017213 SRR29017234 SRR29017242
SRR29017384 SRR29017410 SRR29017445 SRR29017548 SRR29017583 SRR29017611
SRR29017620 SRR29017622 SRR29017648 SRR29017661 SRR29017674 SRR29017683
SRR29017690 SRR29017695 SRR29017701 SRR29055490 SRR29341528 SRR29440700
SRR30443755 SRR35281596 SRR35281598 SRR35281599 SRR35281600 SRR35281602
SRR3675589 SRR4423155 SRR4423179 SRR4423181 SRR6650186 SRR6797722 SRR6797801
""".split())

# Colors
C_L4111 = '#1976D2'  # Blue
C_L4112 = '#E53935'  # Red
C_OTHER = '#9E9E9E'  # Grey for outgroup/basal

random.seed(42)

# ── Load tree ──────────────────────────────────────────────────────────────
print("Loading tree...")
tree = Phylo.read(TREE_FILE, 'newick')
all_tips = list(tree.get_terminals())
print(f"  Total tips: {len(all_tips)}")

# Separate L4.11 tips from outgroup
l411_tips = [t for t in all_tips if t.name and t.name.startswith('L4.11_')]
other_tips = [t for t in all_tips if t not in l411_tips]
print(f"  L4.11 tips: {len(l411_tips)}, outgroup: {len(other_tips)}")

# Classify into sub-lineages
def get_raw_id(tip_name):
    """Extract SRR/ERR ID from tip name like 'L4.11_SRR28393599'."""
    if '_' in tip_name:
        return tip_name.split('_', 1)[1]
    return tip_name

l4111_tips = [t for t in l411_tips if get_raw_id(t.name) in L4111_IDS]
l4112_tips = [t for t in l411_tips if get_raw_id(t.name) not in L4111_IDS]
print(f"  L4.11.1: {len(l4111_tips)}, L4.11.2: {len(l4112_tips)}")

# ── Subsample proportionally ───────────────────────────────────────────────
n1 = max(10, int(MAX_STRAINS * len(l4111_tips) / len(l411_tips)))
n2 = MAX_STRAINS - n1
# Include a few outgroup for contrast
n_out = min(5, len(other_tips))
n2 = n2 - n_out

sample_1 = random.sample(l4111_tips, min(n1, len(l4111_tips)))
sample_2 = random.sample(l4112_tips, min(n2, len(l4112_tips)))
sample_other = random.sample(other_tips, n_out) if other_tips else []

sampled = sample_1 + sample_2 + sample_other
print(f"  Subsampled: {len(sample_1)} L4.11.1 + {len(sample_2)} L4.11.2 + {len(sample_other)} outgroup = {len(sampled)}")

# ── Compute pairwise distance matrix ──────────────────────────────────────
print("Computing pairwise distances (this may take a minute)...")
n = len(sampled)
dist_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i + 1, n):
        d = tree.distance(sampled[i], sampled[j])
        dist_matrix[i, j] = d
        dist_matrix[j, i] = d
    if (i + 1) % 20 == 0:
        print(f"  Row {i+1}/{n}")

print(f"  Distance range: {dist_matrix[dist_matrix > 0].min():.6f} – {dist_matrix.max():.6f}")

# ── Create labels and colors ──────────────────────────────────────────────
labels = []
row_colors = []
for t in sampled:
    raw_id = get_raw_id(t.name)
    if raw_id in L4111_IDS:
        labels.append(raw_id)
        row_colors.append(C_L4111)
    elif t in other_tips:
        labels.append(raw_id)
        row_colors.append(C_OTHER)
    else:
        labels.append(raw_id)
        row_colors.append(C_L4112)

# ── Create clustermap ─────────────────────────────────────────────────────
print("Creating clustermap...")

# Convert to pandas for seaborn
import pandas as pd
df = pd.DataFrame(dist_matrix, index=labels, columns=labels)
row_color_series = pd.Series(row_colors, index=labels, name='Sub-lineage')

# Use a sequential colormap — lower distances = darker
g = sns.clustermap(
    df,
    method='ward',
    metric='euclidean',
    cmap='YlOrRd_r',
    figsize=(12, 11),
    row_colors=row_color_series,
    col_colors=row_color_series,
    xticklabels=False,  # too many for individual labels
    yticklabels=False,
    linewidths=0,
    dendrogram_ratio=(0.15, 0.15),
    cbar_pos=(0.02, 0.03, 0.03, 0.15),
    cbar_kws={'label': 'Patristic distance (subst./site)'},
)

# Legend — place in upper-right area of the figure, clear of dendrogram
legend_patches = [
    mpatches.Patch(facecolor=C_L4111, label=f'L4.11.1 (n={len(sample_1)})'),
    mpatches.Patch(facecolor=C_L4112, label=f'L4.11.2 (n={len(sample_2)})'),
    mpatches.Patch(facecolor=C_OTHER, label=f'Outgroup (n={len(sample_other)})'),
]
g.fig.legend(handles=legend_patches, loc='upper right',
             fontsize=9, frameon=True, facecolor='white',
             edgecolor='#cccccc', bbox_to_anchor=(0.98, 0.98),
             title='Sub-lineage', title_fontsize=10)

plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n✓ Heatmap saved to: {OUTPUT_PNG}")

# ── Print summary statistics ──────────────────────────────────────────────
# Intra-group distances
idx_1 = [i for i, t in enumerate(sampled) if get_raw_id(t.name) in L4111_IDS]
idx_2 = [i for i, t in enumerate(sampled) if t not in other_tips and get_raw_id(t.name) not in L4111_IDS]

intra_1 = [dist_matrix[i, j] for i in idx_1 for j in idx_1 if i < j]
intra_2 = [dist_matrix[i, j] for i in idx_2 for j in idx_2 if i < j]
inter = [dist_matrix[i, j] for i in idx_1 for j in idx_2]

print(f"\nDistance summary:")
print(f"  Intra L4.11.1: mean={np.mean(intra_1):.6f}, std={np.std(intra_1):.6f}")
print(f"  Intra L4.11.2: mean={np.mean(intra_2):.6f}, std={np.std(intra_2):.6f}")
print(f"  Inter (1 vs 2): mean={np.mean(inter):.6f}, std={np.std(inter):.6f}")
print(f"  Ratio inter/intra_1: {np.mean(inter)/np.mean(intra_1):.2f}x")
print(f"  Ratio inter/intra_2: {np.mean(inter)/np.mean(intra_2):.2f}x")
